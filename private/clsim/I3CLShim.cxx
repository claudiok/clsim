
#include <icetray/I3SingleServiceFactory.h>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include "clsim/I3CLShim.h"

I3CLShim::~I3CLShim() {}

I3CLShim::I3CLShim(const I3Context &ctx) : I3ServiceBase(ctx)
{
	AddParameter("GeometryFactory", "", geometryFactory_);
	AddParameter("StepConverter", "", stepConverter_);
	AddParameter("PhotonConverter", "", photonConverter_);
	AddParameter("WavelengthAcceptance", "",wavelengthAcceptance_);
	AddParameter("AngularAcceptance", "", angularAcceptance_);
	AddParameter("OversampleFactor", "", oversampleFactor_);
	AddParameter("CacheDepth", "", cacheDepth_);
}

void
I3CLShim::Configure()
{
	GetParameter("GeometryFactory", geometryFactory_);
	GetParameter("StepConverter", stepConverter_);
	GetParameter("PhotonConverter", photonConverter_);
	GetParameter("WavelengthAcceptance", wavelengthAcceptance_);
	GetParameter("AngularAcceptance", angularAcceptance_);
	GetParameter("OversampleFactor", oversampleFactor_);
	GetParameter("CacheDepth", cacheDepth_);
}

void
I3CLShim::GetMeanAmplitudes(std::vector<LightSource> &sources, const std::vector<Receiver> &receivers)
{
	// Create a potemkin geometry from the provided vector of positions, and compile
	// it into an OpenCL kernel.
	if (!photonConverter_->IsInitialized()) {
		boost::python::list positions;
		BOOST_FOREACH(const I3PhotonicsService::Receiver &r, receivers)
			positions.append(r.pos);
		I3CLSimSimpleGeometryPtr geo = boost::python::extract<I3CLSimSimpleGeometryPtr>(geometryFactory_(positions));
		photonConverter_->SetGeometry(geo);
		for (size_t i=0; i < geo->size(); i++)
			receiverMap_[std::make_pair(geo->GetStringID(i), geo->GetDomID(i))] = i;
		photonConverter_->Initialize();
		
		// Ensure that the step generator and propagator work in batches of the same size
		granularity_ = photonConverter_->GetWorkgroupSize();
		uint32_t bunchSize = photonConverter_->GetMaxNumWorkitems();
		uint32_t maxBunchSize = bunchSize - bunchSize%granularity_;
		
		stepConverter_->SetMaxBunchSize(maxBunchSize);
		stepConverter_->SetBunchSizeGranularity(granularity_);
		stepConverter_->Initialize();
	}
	
#if 0
	for (std::list<CacheEntry>::iterator entry = resultCache_.begin(); entry != resultCache_.end(); entry++) {
		bool equal = (entry->first.size() == sources.size());
		for (unsigned i=0; i < sources.size(); i++)
			equal &= entry->first[i].first == sources[i].first;
		if (equal) {
			receivers = entry->second;
			// Move to front
			if (entry != resultCache_.begin())
				resultCache_.splice(resultCache_.begin(), resultCache_, entry);
			return;
		}
	}
#endif

	double energyScale = 1e4;
	// Enqueue N copies of each light source
	typedef std::pair<PhotonicsSource, double> source_pair;
	uint32_t index = 0;
	
	BOOST_FOREACH(const LightSource &lsource, sources) {
		I3Particle particle(*lsource.source.particle);
		particle.SetEnergy(lsource.source.E*energyScale);
		I3CLSimLightSource source(particle);
		for (int i=0; i < oversampleFactor_; i++, index++)
			stepConverter_->EnqueueLightSource(source, index);
		log_error("Enqueued %.1e GeV source", particle.GetEnergy());
	}
	stepConverter_->EnqueueBarrier();
	
	// Harvest steps and pass them on to the propagator, using
	// the same ID as before.
	uint32_t pindex = 0;
	for (;; pindex++) {
		bool barrierReset = false;
		I3CLSimStepSeriesConstPtr steps = stepConverter_->GetConversionResultWithBarrierInfo(barrierReset);
#if 0
		while (steps->size % granularity_ != 0) {
			steps->push_back(I3CLSimStep());
			steps->back().SetNumPhotons(0);
		}
#endif
		log_info("Got %zu steps", steps->size());
		if (steps->size() > 0)
			photonConverter_->EnqueueSteps(steps, pindex);
		if (barrierReset)
			break;
	}
	
	// Compute the offset of each receiver 
	std::vector<uint32_t> receiver_offsets;
	{
		uint32_t offset = 0;
		BOOST_FOREACH(const Receiver &receiver, receivers) {
			receiver_offsets.push_back(offset);
			offset += std::max(std::distance(receiver.time_bin_edges.first, receiver.time_bin_edges.second)-1, ptrdiff_t(0));
			assert(std::distance(receiver.time_bin_edges.first, receiver.time_bin_edges.second) >= 2);
		}
	}
	
	// Harvest photons, adding their hit probability weight to the appropriate bin
	// in each receiver.
	for (uint32_t j=0; j < pindex+1; j++) {
		I3CLSimStepToPhotonConverter::ConversionResult_t result = photonConverter_->GetConversionResult();
		log_info("Got %zu photons", result.photons->size());
		BOOST_FOREACH(const I3CLSimPhoton &p, *result.photons) {
			// look up the index in the simple geometry from the [potemkin] string/om ids
			std::map<std::pair<int16_t, uint16_t>, size_t>::const_iterator
			    receiver_index = receiverMap_.find(std::make_pair(p.GetStringID(), p.GetOMID()));
			if (receiver_index == receiverMap_.end())
				continue;
			const I3PhotonicsService::Receiver &receiver = receivers[receiver_index->second];
			// ensure that we actually have bins
			assert(receiver.time_bin_edges.second >= receiver.time_bin_edges.first+1);
			// bail if the photon arrives before or after readout
			if (p.GetTime() < *receiver.time_bin_edges.first || p.GetTime() >= *(receiver.time_bin_edges.second-1))
				continue;
			// find the time bin this photon belongs in
			int binidx = std::max(std::distance(receiver.time_bin_edges.first,
			    std::lower_bound(receiver.time_bin_edges.first, receiver.time_bin_edges.second-1, p.GetTime()))-1, ptrdiff_t(0));
			// each real source is repeated oversampleFactor_ times
			int sourceidx = p.GetID() / oversampleFactor_;
			// ensure that we haven't stepped out of bounds
			if (receiver_index == boost::prior(receiverMap_.end()))
				assert(receiver_offsets[receiver_index->second] + binidx < sources[sourceidx].amplitudes.size());
			else {
				assert(receiver_offsets[receiver_index->second] + binidx < receiver_offsets[receiver_index->second+1]);
			}
			// compute detection weight for wavelength and impact angle
			double weight = p.GetWeight()*wavelengthAcceptance_->GetValue(p.GetWavelength())*angularAcceptance_->GetValue(std::cos(p.GetDirTheta()));
			// add to the bin
			sources[sourceidx].amplitudes(receiver_offsets[receiver_index->second] + binidx)
			    += weight/(oversampleFactor_*energyScale);
		}
	}

#if 0
	// Cache the result if we have room.
	if (cacheDepth_ > 0) {
		if (resultCache_.size() < cacheDepth_)
			resultCache_.push_back(CacheEntry());
		resultCache_.back() = CacheEntry(sources, receivers);
	}
#endif
}

typedef I3SingleServiceFactory<I3CLShim, I3BlockPhotonicsService>
    I3CLShimFactory;
I3_SERVICE_FACTORY(I3CLShimFactory);

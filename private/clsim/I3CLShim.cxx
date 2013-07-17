
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
I3CLShim::GetMeanAmplitudes(const std::vector<std::pair<PhotonicsSource, double> > &sources, std::vector<I3PhotonicsService::Receiver> &receivers)
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
	
	double energyScale = 1e4;
	// Enqueue N copies of each light source
	typedef std::pair<PhotonicsSource, double> source_pair;
	uint32_t index = 0;
	
	BOOST_FOREACH(const source_pair &spair, sources) {
		I3Particle particle(*spair.first.particle);
		particle.SetEnergy(spair.first.E*energyScale);
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
	
	// Harvest photons, adding their hit probability weight to the appropriate bin
	// in each receiver.
	for (uint32_t j=0; j < pindex; j++) {
		I3CLSimStepToPhotonConverter::ConversionResult_t result = photonConverter_->GetConversionResult();
		log_info("Got %zu photons", result.photons->size());
		BOOST_FOREACH(const I3CLSimPhoton &p, *result.photons) {
			I3PhotonicsService::Receiver &receiver = receivers[receiverMap_[std::make_pair(p.GetStringID(), p.GetOMID())]];
			if (p.GetTime() < receiver.time_edges.front() || p.GetTime() >= receiver.time_edges.back())
				continue;
			int binidx = std::distance(receiver.time_edges.begin(), std::lower_bound(receiver.time_edges.begin(), receiver.time_edges.end(), p.GetTime()));
			int sourceidx = p.GetID() / oversampleFactor_;
			double weight = p.GetWeight()*wavelengthAcceptance_->GetValue(p.GetWavelength())*angularAcceptance_->GetValue(std::cos(p.GetDirTheta()));
			receiver.amplitudes(sourceidx, binidx == 0 ? 0 : binidx-1) += weight/(oversampleFactor_*energyScale);
		}
	}
	
	// Cache the result if we have room.
	if (cacheDepth_ > 0) {
		if (resultCache_.size() < cacheDepth_)
			resultCache_.push_back(CacheEntry());
		resultCache_.back() = CacheEntry(sources, receivers);
	}
}

typedef I3SingleServiceFactory<I3CLShim, I3PhotonicsService>
    I3CLShimFactory;
I3_SERVICE_FACTORY(I3CLShimFactory);

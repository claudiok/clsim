
#include <icetray/I3Module.h>
#include <phys-services/I3RandomService.h>
#include <dataclasses/physics/I3MCTree.h>


#include "clsim/I3CLSimOpenCLDevice.h"
#include "clsim/I3CLSimLightSourceParameterization.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimModuleHelper.h"
#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
#include "clsim/tabulator/Axes.h"

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>

#include <fstream>

namespace fs = boost::filesystem;

class I3CLSimTabulatorModule : public I3Module {
public:
	I3CLSimTabulatorModule(const I3Context&);
	virtual ~I3CLSimTabulatorModule();
	void Configure();
	
	void DAQ(I3FramePtr);
	void Finish();

private:
	void HarvestSteps();
	
	I3CLSimLightSourceParameterizationSeries parameterizationList_;
	I3RandomServicePtr randomService_;
	I3CLSimFunctionConstPtr wavelengthGenerationBias_;
	I3CLSimMediumPropertiesConstPtr mediumProperties_;
	I3CLSimFunctionConstPtr angularAcceptance_;
	// I3CLSimSpectrumTableConstPtr spectrumTable_;
	I3CLSimOpenCLDeviceSeries openCLDeviceList_;
	
	I3CLSimLightSourceToStepConverterPtr particleToStepsConverter_;
	I3CLSimStepToTableConverterPtr tabulator_;
	
	std::string tablePath_;
	boost::python::dict tableHeader_;
	clsim::tabulator::AxesPtr axes_;
	
	boost::thread stepHarvester_;
	struct {
		boost::mutex mutex;
		boost::condition_variable cv;
	} semaphore;
	I3CLSimQueue<I3ParticleConstPtr> sourceQueue_;
	bool run_;
	
	SET_LOGGER("I3CLSimTabulatorModule");
};

I3_MODULE(I3CLSimTabulatorModule);

I3CLSimTabulatorModule::I3CLSimTabulatorModule(const I3Context &ctx)
    : I3Module(ctx), run_(true)
{
	AddOutBox("OutBox");
	
	AddParameter("RandomService", "", "I3RandomService");
	AddParameter("WavelengthAcceptance", "", wavelengthGenerationBias_);
	AddParameter("AngularAcceptance", "", angularAcceptance_);
	AddParameter("MediumProperties", "", mediumProperties_);
	AddParameter("ParameterizationList","", parameterizationList_);
	AddParameter("OpenCLDeviceList", "", openCLDeviceList_);
	AddParameter("Filename", "", "");
	AddParameter("TableHeader", "", boost::python::dict());
	AddParameter("Axes", "", axes_);
}

void I3CLSimTabulatorModule::Configure()
{
	GetParameter("RandomService", randomService_);
	GetParameter("WavelengthAcceptance", wavelengthGenerationBias_);
	GetParameter("AngularAcceptance", angularAcceptance_);
	GetParameter("MediumProperties", mediumProperties_);
	GetParameter("ParameterizationList",parameterizationList_);
	GetParameter("OpenCLDeviceList",openCLDeviceList_);
	GetParameter("Filename", tablePath_);
	GetParameter("TableHeader", tableHeader_);
	GetParameter("Axes", axes_);
	
	if (tablePath_.empty())
		log_fatal("You must specify an output filename!");
	else if (fs::exists(tablePath_))
		log_fatal_stream(tablePath_ << " already exists!");
	try {
		std::ofstream f(tablePath_.c_str());
	} catch (...) {
		log_fatal_stream("Could not open " << tablePath_ << " for writing");
	}
	fs::remove(tablePath_);
	
	tabulator_ = boost::make_shared<I3CLSimStepToTableConverter>(
	    openCLDeviceList_[0], axes_, mediumProperties_,
	    wavelengthGenerationBias_, angularAcceptance_, randomService_);
	
	particleToStepsConverter_ =
	    I3CLSimModuleHelper::initializeGeant4(randomService_,
	                                          mediumProperties_,
	                                          wavelengthGenerationBias_,
	                                          tabulator_->GetBunchSize() /*granularity*/,
	                                          tabulator_->GetBunchSize() /*maxBunchSize*/,
	                                          parameterizationList_,
	                                          "QGSP_BERT_EMV" /*geant4PhysicsListName_*/,
	                                          10.*I3Units::perCent /*geant4MaxBetaChangePerStep_*/,
	                                          200 /*geant4MaxNumPhotonsPerStep_*/,
	                                          false); // the multiprocessor version is not yet safe to use

	stepHarvester_ = boost::thread(boost::bind(&I3CLSimTabulatorModule::HarvestSteps, this));
	
	
}

void I3CLSimTabulatorModule::Finish()
{
	log_trace("finish called");
	run_ = false;
	if (!particleToStepsConverter_->BarrierActive())
		particleToStepsConverter_->EnqueueBarrier();
	stepHarvester_.join();
	tabulator_->Finish();
	
	tabulator_->WriteFITSFile(tablePath_, tableHeader_);
}

I3CLSimTabulatorModule::~I3CLSimTabulatorModule()
{
	log_trace("dtor called");
}

/// Harvest steps and feed them to the tabulator
void I3CLSimTabulatorModule::HarvestSteps()
{
	// Get the first reference source. This will block until something is
	// added to the queue.
	I3ParticleConstPtr reference = sourceQueue_.Get();
	for (;;) {
		I3CLSimStepSeriesConstPtr steps;
		bool barrierWasJustReset=false;
		steps = particleToStepsConverter_->GetConversionResultWithBarrierInfo(barrierWasJustReset);
		
		if (steps && !steps->empty()) {
			// Do stuff
			tabulator_->EnqueueSteps(I3CLSimStepSeriesPtr(new I3CLSimStepSeries(*steps)), reference);
			log_debug_stream("enqueued " << steps->size() << " steps");
		}
		
		if (barrierWasJustReset) {
			if (!run_) {
				log_trace("Exiting on barrier");
				return;
			} else {
				// Signal main thread to continue
				semaphore.cv.notify_one();
				// Get the reference source for the next frame
				reference = sourceQueue_.Get();
			}
		}
	}
	
}

void I3CLSimTabulatorModule::DAQ(I3FramePtr frame)
{
	I3ParticleConstPtr reference = frame->Get<I3ParticleConstPtr>("ReferenceParticle");
	if (!reference)
		log_fatal("Frame does not contain an I3Particle 'ReferenceParticle'!");
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	
	{
		// Wait for particles from the previous event to be consumed
		boost::unique_lock<boost::mutex> lock(semaphore.mutex);
		while (particleToStepsConverter_->BarrierActive())
			semaphore.cv.wait(lock);
	}
	
	sourceQueue_.Put(reference);
	BOOST_FOREACH(const I3Particle &p, *mctree) {
		if (p.GetShape() != I3Particle::Dark && p.GetLocationType() == I3Particle::InIce)
			particleToStepsConverter_->EnqueueLightSource(I3CLSimLightSource(p), 0);
	}
	particleToStepsConverter_->EnqueueBarrier();
}


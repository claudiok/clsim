
#include <icetray/I3Module.h>
#include <phys-services/I3RandomService.h>
#include <dataclasses/physics/I3MCTree.h>


#include "clsim/I3CLSimOpenCLDevice.h"
#include "clsim/I3CLSimLightSourceParameterization.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimModuleHelper.h"

#include <boost/foreach.hpp>
#include <boost/thread.hpp>

class I3CLSimTabulatorModule : public I3Module {
public:
	I3CLSimTabulatorModule(const I3Context&);
	void Configure();
	
	void DAQ(I3FramePtr);
	void Finish();

private:
	void HarvestSteps();
	
	I3CLSimLightSourceParameterizationSeries parameterizationList_;
	I3RandomServicePtr randomService_;
	I3CLSimFunctionConstPtr wavelengthGenerationBias_;
	I3CLSimMediumPropertiesConstPtr mediumProperties_;
	// I3CLSimSpectrumTableConstPtr spectrumTable_;
	I3CLSimOpenCLDeviceSeries openCLDeviceList_;
	
	I3CLSimLightSourceToStepConverterPtr particleToStepsConverter_;
	
	uint32_t sourceCounter_;
	boost::thread stepHarvester_;
};

I3_MODULE(I3CLSimTabulatorModule);

I3CLSimTabulatorModule::I3CLSimTabulatorModule(const I3Context &ctx)
    : I3Module(ctx)
{
	AddOutBox("OutBox");
	
	AddParameter("RandomService", "", "I3RandomService");
	AddParameter("WavelengthGenerationBias", "", wavelengthGenerationBias_);
	AddParameter("MediumProperties", "", mediumProperties_);
	AddParameter("ParameterizationList","", parameterizationList_);
	AddParameter("OpenCLDeviceList", "", openCLDeviceList_);
}

void I3CLSimTabulatorModule::Configure()
{
	GetParameter("RandomService", randomService_);
	GetParameter("WavelengthGenerationBias", wavelengthGenerationBias_);
	GetParameter("MediumProperties", mediumProperties_);
	GetParameter("ParameterizationList",parameterizationList_);
	GetParameter("OpenCLDeviceList",openCLDeviceList_);
	
	particleToStepsConverter_ =
	    I3CLSimModuleHelper::initializeGeant4(randomService_,
	                                          mediumProperties_,
	                                          wavelengthGenerationBias_,
	                                          1 /*granularity*/,
	                                          1 /*maxBunchSize*/,
	                                          parameterizationList_,
	                                          "QGSP_BERT_EMV" /*geant4PhysicsListName_*/,
	                                          10.*I3Units::perCent /*geant4MaxBetaChangePerStep_*/,
	                                          200 /*geant4MaxNumPhotonsPerStep_*/,
	                                          false); // the multiprocessor version is not yet safe to use
	sourceCounter_ = 0;
	
	stepHarvester_ = boost::thread(boost::bind(&I3CLSimTabulatorModule::HarvestSteps, this));
	
}

void I3CLSimTabulatorModule::Finish()
{
	stepHarvester_.interrupt();
	stepHarvester_.join();
}

/// Harvest steps and feed them to OpenCL
void I3CLSimTabulatorModule::HarvestSteps()
{
	for (;;) {
		I3CLSimStepSeriesConstPtr steps;
		bool barrierWasJustReset=false;
		try {
			steps = particleToStepsConverter_->GetConversionResultWithBarrierInfo(barrierWasJustReset);
		} catch (boost::thread_interrupted &interrupt) {
			break;
		}
		if (steps && ! steps->empty()) {
			// Do stuff
		}
	}
	
}

void I3CLSimTabulatorModule::DAQ(I3FramePtr frame)
{
	
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	
	BOOST_FOREACH(const I3Particle &p, *mctree) {
		if (p.GetShape() != I3Particle::Dark && p.GetLocationType() == I3Particle::InIce)
			particleToStepsConverter_->EnqueueLightSource(I3CLSimLightSource(p), sourceCounter_++);
	}

}


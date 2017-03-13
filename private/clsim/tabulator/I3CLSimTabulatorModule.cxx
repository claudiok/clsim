/**
 * Copyright (c) 2015
 * Jakob van Santen <jvansanten@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimStepToTableConverter.cxx
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#include <icetray/I3Module.h>
#include <phys-services/I3RandomService.h>
#include <dataclasses/physics/I3MCTree.h>

#include "clsim/I3CLSimOpenCLDevice.h"
#include "clsim/I3CLSimLightSourceParameterization.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSpectrumTable.h"
#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimModuleHelper.h"
#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
#include "clsim/tabulator/Axes.h"

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/make_shared.hpp>
#include <boost/scoped_ptr.hpp>
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
	
	std::string mctreeName_, flasherPulseSeriesName_;
	I3CLSimLightSourceParameterizationSeries parameterizationList_;
	I3RandomServicePtr randomService_;
	I3CLSimFunctionConstPtr wavelengthGenerationBias_;
	I3CLSimMediumPropertiesConstPtr mediumProperties_;
	I3CLSimFunctionConstPtr angularAcceptance_;
	I3CLSimSpectrumTableConstPtr spectrumTable_;
	I3CLSimOpenCLDeviceSeries openCLDeviceList_;
	double referenceArea_;
	size_t photonsPerBunch_, entriesPerPhoton_;
	bool recordErrors_;
	
	I3CLSimLightSourceToStepConverterPtr particleToStepsConverter_;
	boost::scoped_ptr<I3CLSimStepToTableConverter> tabulator_;
	
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
	
	AddParameter("MCTreeName", "", "I3MCTree");
	AddParameter("FlasherPulseSeriesName", "", "");
	AddParameter("RandomService", "", "I3RandomService");
	AddParameter("Area", "Geometric area of the sensor", M_PI*std::pow(0.16510*I3Units::m, 2));
	AddParameter("WavelengthAcceptance", "Wavelength acceptance (relative to geometric area)", wavelengthGenerationBias_);
	AddParameter("AngularAcceptance", "Angular acceptance (relative to head-on geometric area)", angularAcceptance_);
	AddParameter("MediumProperties", "", mediumProperties_);
	AddParameter("ParameterizationList","", parameterizationList_);
	AddParameter("SpectrumTable", "", spectrumTable_);
	AddParameter("OpenCLDeviceList", "", openCLDeviceList_);
	AddParameter("PhotonsPerBunch", "", 200);
	AddParameter("EntriesPerPhoton", "", 3000);
	AddParameter("Filename", "", "");
	AddParameter("RecordErrors", "", false);
	AddParameter("TableHeader", "", boost::python::dict());
	AddParameter("Axes", "", axes_);
}

void I3CLSimTabulatorModule::Configure()
{
	GetParameter("MCTreeName", mctreeName_);
	GetParameter("FlasherPulseSeriesName", flasherPulseSeriesName_);
	GetParameter("RandomService", randomService_);
	GetParameter("WavelengthAcceptance", wavelengthGenerationBias_);
	GetParameter("Area", referenceArea_);
	GetParameter("AngularAcceptance", angularAcceptance_);
	GetParameter("MediumProperties", mediumProperties_);
	GetParameter("ParameterizationList",parameterizationList_);
	GetParameter("SpectrumTable", spectrumTable_);
	GetParameter("OpenCLDeviceList",openCLDeviceList_);
	GetParameter("PhotonsPerBunch", photonsPerBunch_);
	GetParameter("EntriesPerPhoton", entriesPerPhoton_);
	GetParameter("Filename", tablePath_);
	GetParameter("RecordErrors", recordErrors_);
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
	
	if (openCLDeviceList_.size() == 0)
		log_fatal_stream("No OpenCL devices provided. Does your OpenCL runtime support CPUs?");
	
	tabulator_.reset(
	    new I3CLSimStepToTableConverter(
	    openCLDeviceList_[0], axes_, entriesPerPhoton_*photonsPerBunch_,
	    recordErrors_,
	    mediumProperties_, spectrumTable_, referenceArea_,
	    wavelengthGenerationBias_, angularAcceptance_, randomService_));
	
	particleToStepsConverter_ =
	    I3CLSimModuleHelper::initializeGeant4(randomService_,
	                                          mediumProperties_,
	                                          wavelengthGenerationBias_,
	                                          tabulator_->GetBunchSize() /*granularity*/,
	                                          tabulator_->GetBunchSize() /*maxBunchSize*/,
	                                          parameterizationList_,
	                                          "QGSP_BERT_EMV" /*geant4PhysicsListName_*/,
	                                          10.*I3Units::perCent /*geant4MaxBetaChangePerStep_*/,
	                                          photonsPerBunch_ /*geant4MaxNumPhotonsPerStep_*/,
	                                          false); // the multiprocessor version is not yet safe to use

	stepHarvester_ = boost::thread(boost::bind(&I3CLSimTabulatorModule::HarvestSteps, this));
	
	
}

namespace {
    class ScopedGILRelease
    {
    public:
        inline ScopedGILRelease()
        {
            m_thread_state = PyEval_SaveThread();
        }
        
        inline ~ScopedGILRelease()
        {
            PyEval_RestoreThread(m_thread_state);
            m_thread_state = NULL;
        }
        
    private:
        PyThreadState *m_thread_state;
    };    
}


void I3CLSimTabulatorModule::Finish()
{
	log_trace("finish called");
	run_ = false;
	if (!particleToStepsConverter_->BarrierActive())
		particleToStepsConverter_->EnqueueBarrier();
	{
		// Release the GIL while we wait
		ScopedGILRelease release;
		stepHarvester_.join();
		tabulator_->Finish();
	}
	
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
			log_trace_stream("enqueued " << steps->size() << " steps");
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
	I3MCTreeConstPtr mctree;
	I3CLSimFlasherPulseSeriesConstPtr flashers;
	I3ParticleConstPtr reference = frame->Get<I3ParticleConstPtr>("ReferenceParticle");
	if (!reference)
		log_fatal("Frame does not contain an I3Particle 'ReferenceParticle'!");
	
	{
		// Release the Python interpreter lock
		ScopedGILRelease release;
		// Wait for steps from the previous event to be consumed
		boost::unique_lock<boost::mutex> lock(semaphore.mutex);
		while (particleToStepsConverter_->BarrierActive())
			semaphore.cv.wait(lock);
	}
	
	// Enqueue a copy to ensure that the deleter does not invoke the Python interpreter
	sourceQueue_.Put(I3ParticlePtr(new I3Particle(*reference)));
	if ((mctree = frame->Get<I3MCTreeConstPtr>(mctreeName_))) {
		BOOST_FOREACH(const I3Particle &p, *mctree) {
			if (p.GetShape() != I3Particle::Dark && p.GetLocationType() == I3Particle::InIce)
				particleToStepsConverter_->EnqueueLightSource(I3CLSimLightSource(p), 0);
		}
	}
	if ((flashers = frame->Get<I3CLSimFlasherPulseSeriesConstPtr>(flasherPulseSeriesName_))) {
		BOOST_FOREACH(const I3CLSimFlasherPulse &p, *flashers) {
			particleToStepsConverter_->EnqueueLightSource(I3CLSimLightSource(p), 0);
			log_trace_stream("enqueued a pulse of "<<p.GetNumberOfPhotonsNoBias()<<" photons");
		}
	}

	particleToStepsConverter_->EnqueueBarrier();
}


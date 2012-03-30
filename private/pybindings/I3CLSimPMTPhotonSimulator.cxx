/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
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
 * @file I3CLSimPMTPhotonSimulator.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimPMTPhotonSimulator.h>

#ifdef USE_HIT_MAKER
#include <clsim/I3CLSimPMTPhotonSimulatorIceCube.h>
#endif

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimPMTPhotonSimulatorWrapper : I3CLSimPMTPhotonSimulator, bp::wrapper<I3CLSimPMTPhotonSimulator>
{
    // pure virtual
    virtual void ApplyAfterPulseLatePulseAndJitterSim(const OMKey &key, const I3MCHit &input_hit, std::vector<I3MCHit> &output_vector) const
    {
        this->get_override("ApplyAfterPulseLatePulseAndJitterSim")(key, input_hit, output_vector);
    }
    
    virtual void SetCalibration(I3CalibrationConstPtr calibration)
    {
        this->get_override("SetCalibration")(calibration);
    }
    
    virtual void SetDetectorStatus(I3DetectorStatusConstPtr status)
    {
        this->get_override("SetDetectorStatus")(status);
    }

    virtual void SetRandomService(I3RandomServicePtr random)
    {
        this->get_override("SetRandomService")(random);
    }


};

void register_I3CLSimPMTPhotonSimulator()
{
    {
        bp::scope I3CLSimPMTPhotonSimulator_scope = 
        bp::class_<I3CLSimPMTPhotonSimulatorWrapper, boost::shared_ptr<I3CLSimPMTPhotonSimulatorWrapper>, boost::noncopyable>("I3CLSimPMTPhotonSimulator", no_init)
        .def("ApplyAfterPulseLatePulseAndJitterSim", bp::pure_virtual(&I3CLSimPMTPhotonSimulator::ApplyAfterPulseLatePulseAndJitterSim))
        .def("SetCalibration", bp::pure_virtual(&I3CLSimPMTPhotonSimulator::SetCalibration))
        .def("SetDetectorStatus", bp::pure_virtual(&I3CLSimPMTPhotonSimulator::SetDetectorStatus))
        .def("SetRandomService", bp::pure_virtual(&I3CLSimPMTPhotonSimulator::SetRandomService))
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorWrapper>, shared_ptr<const I3CLSimPMTPhotonSimulator> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorWrapper>, shared_ptr<I3CLSimPMTPhotonSimulator> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorWrapper>, shared_ptr<const I3CLSimPMTPhotonSimulatorWrapper> >();
    utils::register_const_ptr<I3CLSimPMTPhotonSimulator>();
    
#ifdef USE_HIT_MAKER
    // IceCube code (calls hit-maker)
    {
        bp::class_<
        I3CLSimPMTPhotonSimulatorIceCube, 
        boost::shared_ptr<I3CLSimPMTPhotonSimulatorIceCube>, 
        bases<I3CLSimPMTPhotonSimulator>,
        boost::noncopyable
        >
        (
         "I3CLSimPMTPhotonSimulatorIceCube",
         bp::init<
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("jitter")=I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_jitter,
            bp::arg("pre_pulse_probability")=I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_pre_pulse_probability,
            bp::arg("late_pulse_probability")=I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_late_pulse_probability,
            bp::arg("after_pulse_probability")=I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_after_pulse_probability
           )
          )
        )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorIceCube>, shared_ptr<const I3CLSimPMTPhotonSimulatorIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorIceCube>, shared_ptr<I3CLSimPMTPhotonSimulator> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimPMTPhotonSimulatorIceCube>, shared_ptr<const I3CLSimPMTPhotonSimulator> >();
    utils::register_const_ptr<I3CLSimPMTPhotonSimulatorIceCube>();
#endif
    
}

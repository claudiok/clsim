//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module g4sim-interface.
//
//   g4sim-intrface is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <sstream>

#include <clsim/I3CLSimPMTPhotonSimulator.h>

#ifdef USE_HIT_MAKER
#include <clsim/I3CLSimPMTPhotonSimulatorIceCube.h>
#endif

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>
#include "const_ptr_helpers.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimPMTPhotonSimulatorWrapper : I3CLSimPMTPhotonSimulator, bp::wrapper<I3CLSimPMTPhotonSimulator>
{
    // pure virtual
    virtual void ApplyAfterPulseLatePulseAndJitterSim(const I3MCHit &input_hit, std::vector<I3MCHit> &output_vector) const
    {
        this->get_override("ApplyAfterPulseLatePulseAndJitterSim")(input_hit, output_vector);
    }
    
    virtual void SetCalibration(I3CalibrationConstPtr calibration)
    {
        this->get_override("SetCalibration")(calibration);
    }
    
    virtual void SetDetectorStatus(I3DetectorStatusConstPtr status)
    {
        this->get_override("SetDetectorStatus")(status);
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
         I3RandomServicePtr,
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("randomService"),
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

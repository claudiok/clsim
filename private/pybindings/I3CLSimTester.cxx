//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module g4sim-interface.
//
//   this file is free software; you can redistribute it and/or modify
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

#include <test/I3CLSimTesterBase.h>
#include <test/I3CLSimRandomDistributionTester.h>
#include <test/I3CLSimWlenDependentValueTester.h>
#include <test/I3CLSimMediumPropertiesTester.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

#include <boost/foreach.hpp>

using namespace boost::python;
namespace bp = boost::python;


struct I3CLSimTesterBaseWrapper : I3CLSimTesterBase, bp::wrapper<I3CLSimTesterBase>
{

};

void register_I3CLSimTester()
{
    {
        bp::scope I3CLSimTesterBase_scope = 
        bp::class_<I3CLSimTesterBaseWrapper, boost::shared_ptr<I3CLSimTesterBaseWrapper>, boost::noncopyable>
        ("I3CLSimTesterBase", bp::no_init)

        .def("GetMaxWorkgroupSize", &I3CLSimTesterBase::GetMaxWorkgroupSize)
        .add_property("maxWorkgroupSize", &I3CLSimTesterBase::GetMaxWorkgroupSize)
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimTesterBaseWrapper>, shared_ptr<const I3CLSimTesterBase> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimTesterBaseWrapper>, shared_ptr<I3CLSimTesterBase> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimTesterBaseWrapper>, shared_ptr<const I3CLSimTesterBaseWrapper> >();

    // I3CLSimRandomDistributionTester
    {
        bp::scope I3CLSimRandomDistributionTester_scope = 
        bp::class_<I3CLSimRandomDistributionTester, 
                   boost::shared_ptr<I3CLSimRandomDistributionTester>,
                   bases<I3CLSimTesterBase>,
                   boost::noncopyable>
        ("I3CLSimRandomDistributionTester",
         bp::init<const I3CLSimOpenCLDevice &, uint64_t, uint64_t, I3RandomServicePtr, I3CLSimRandomValueConstPtr>
         (
          (
           bp::arg("device"),
           bp::arg("workgroupSize"),
           bp::arg("workItemsPerIteration"),
           bp::arg("randomService"),
           bp::arg("randomDistribution")
          )
         )
        )
        .def("GenerateRandomNumbers", &I3CLSimRandomDistributionTester::GenerateRandomNumbers, bp::arg("iterations"))
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomDistributionTester>, shared_ptr<const I3CLSimRandomDistributionTester> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomDistributionTester>, shared_ptr<I3CLSimTesterBase> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomDistributionTester>, shared_ptr<const I3CLSimTesterBase> >();

    
    // I3CLSimWlenDependentValueTester
    {
        bp::scope I3CLSimWlenDependentValueTester_scope = 
        bp::class_<I3CLSimWlenDependentValueTester, 
        boost::shared_ptr<I3CLSimWlenDependentValueTester>,
        bases<I3CLSimTesterBase>,
        boost::noncopyable>
        ("I3CLSimWlenDependentValueTester",
         bp::init<const I3CLSimOpenCLDevice &, uint64_t, uint64_t, I3CLSimWlenDependentValueConstPtr>
         (
          (
           bp::arg("device"),
           bp::arg("workgroupSize"),
           bp::arg("workItemsPerIteration"),
           bp::arg("wlenDependentValue")
           )
          )
         )
        .def("EvaluateFunction", &I3CLSimWlenDependentValueTester::EvaluateFunction, bp::arg("xValues"))
        .def("EvaluateDerivative", &I3CLSimWlenDependentValueTester::EvaluateDerivative, bp::arg("xValues"))

        .def("EvaluateReferenceFunction", &I3CLSimWlenDependentValueTester::EvaluateReferenceFunction, bp::arg("xValues"))
        .def("EvaluateReferenceDerivative", &I3CLSimWlenDependentValueTester::EvaluateReferenceDerivative, bp::arg("xValues"))
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueTester>, shared_ptr<const I3CLSimWlenDependentValueTester> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueTester>, shared_ptr<I3CLSimTesterBase> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueTester>, shared_ptr<const I3CLSimTesterBase> >();

    
    // I3CLSimMediumPropertiesTester
    {
        bp::scope I3CLSimMediumPropertiesTester_scope = 
        bp::class_<I3CLSimMediumPropertiesTester, 
        boost::shared_ptr<I3CLSimMediumPropertiesTester>,
        bases<I3CLSimTesterBase>,
        boost::noncopyable>
        ("I3CLSimMediumPropertiesTester",
         bp::init<const I3CLSimOpenCLDevice &, uint64_t, uint64_t, I3CLSimMediumPropertiesConstPtr, I3RandomServicePtr>
         (
          (
           bp::arg("platformAndDeviceName"),
           bp::arg("workgroupSize"),
           bp::arg("workItemsPerIteration"),
           bp::arg("mediumProperties"),
           bp::arg("randomService") = I3RandomServicePtr()
           )
          )
         )
        .def("EvaluatePhaseRefIndex", &I3CLSimMediumPropertiesTester::EvaluatePhaseRefIndex, bp::arg("xValues"), bp::arg("layer"))
        .def("EvaluateDispersion", &I3CLSimMediumPropertiesTester::EvaluateDispersion, bp::arg("xValues"), bp::arg("layer"))
        .def("EvaluateGroupVelocity", &I3CLSimMediumPropertiesTester::EvaluateGroupVelocity, bp::arg("xValues"), bp::arg("layer"))
        .def("EvaluateAbsorptionLength", &I3CLSimMediumPropertiesTester::EvaluateAbsorptionLength, bp::arg("xValues"), bp::arg("layer"))
        .def("EvaluateScatteringLength", &I3CLSimMediumPropertiesTester::EvaluateScatteringLength, bp::arg("xValues"), bp::arg("layer"))
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimMediumPropertiesTester>, shared_ptr<const I3CLSimMediumPropertiesTester> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimMediumPropertiesTester>, shared_ptr<I3CLSimTesterBase> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimMediumPropertiesTester>, shared_ptr<const I3CLSimTesterBase> >();

}

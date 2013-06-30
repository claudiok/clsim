
#include "clsim/I3CLShim.h"

namespace bp = boost::python;

void register_I3CLShim()
{
	bp::import("icecube.photonics_service");
	
	bp::class_<I3CLShim, boost::shared_ptr<I3CLShim>,
	    bp::bases<I3PhotonicsService>, boost::noncopyable>
	    ("I3CLShim", bp::init<const I3Context &>())
	    .def("SetParameter", &I3CLShim::SetParameter)
	    .def("Configure", &I3CLShim::Configure)
	;
}

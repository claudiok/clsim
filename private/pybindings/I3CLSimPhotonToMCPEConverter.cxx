
#include "clsim/dom/I3CLSimPhotonToMCPEConverter.h"
#include "clsim/dom/I3PhotonToMCPEConverter.h"

namespace bp = boost::python;

void register_I3CLSimPhotonToMCPEConverter()
{
    bp::class_<I3CLSimPhotonToMCPEConverter, boost::shared_ptr<I3CLSimPhotonToMCPEConverter>, boost::noncopyable>("I3CLSimPhotonToMCPEConverter", bp::no_init)
    ;
}

void register_I3CLSimPhotonToMCPEConverterForDOMs()
{
    bp::class_<I3CLSimPhotonToMCPEConverterForDOMs, boost::shared_ptr<I3CLSimPhotonToMCPEConverterForDOMs>,
        bp::bases<I3CLSimPhotonToMCPEConverter>, boost::noncopyable>
        ("I3CLSimPhotonToMCPEConverterForDOMs",
          bp::init<I3RandomServicePtr,boost::shared_ptr<const std::map<OMKey, I3CLSimFunctionConstPtr>>,I3CLSimFunctionConstPtr>(
                   bp::args("randomService", "wavelengthAcceptance", "angularAcceptance")))
    ;
}
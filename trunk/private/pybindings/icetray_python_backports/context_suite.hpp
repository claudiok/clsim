// Construct python objects from shared pointers in the context

#ifndef ICETRAY_PYTHON_CONTEXT_SUITE_HPP_INCLUDED_SUITE_HPP_INCLUDED
#define ICETRAY_PYTHON_CONTEXT_SUITE_HPP_INCLUDED_SUITE_HPP_INCLUDED

#include <boost/python/def_visitor.hpp>
#include <icetray/I3Context.h>
#include <icetray/name_of.h>

namespace icetray { namespace python {

template <typename ServiceType>
class context_suite : public boost::python::def_visitor<context_suite<ServiceType > > {
public:

	template <typename Class>
	static void
	visit (Class& cl)
	{
		cl
			.def("from_context", &context_suite<ServiceType>::GetFromContext)
			.staticmethod("from_context")
			.def("put_in_context", &context_suite<ServiceType>::PutInContext)
			;
	}

private:

	typedef boost::shared_ptr<ServiceType> ServicePtr;

	static ServicePtr
	GetFromContext(const I3Context &ctx, const std::string &name)
	{
		ServicePtr service = ctx.Get<ServicePtr>(name);
		if (!service) {
			std::string mesg("The context contains no object of type '");
			mesg += icetray::name_of<ServiceType>();
			mesg += "' at key '";
			mesg += name;
			mesg += "'.";
			PyErr_SetString(PyExc_KeyError, mesg.c_str());
			throw boost::python::error_already_set();
		}
		
		return service;
	}
	
	bool
	PutInContext(ServicePtr self, I3Context &ctx, const std::string &name)
	{
		return ctx.Put(self, name);
	}
	
};

}} // namespace icetray::python

#endif



#include "clsim/backports/I3Matrix.h"
#include <boost/make_shared.hpp>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

namespace {

class pyerr_fmt : public std::ostringstream {
public:
	pyerr_fmt(PyObject *err = PyExc_TypeError) : err_(err) {};
	
	void raise(PyObject *err = PyExc_TypeError) throw(bp::error_already_set)
	{
		PyErr_SetString(err, this->str().c_str());
		throw bp::error_already_set();
	}
	
private:
	PyObject *err_;
};

std::string
type(bp::object &obj)
{
	bp::object typ(bp::handle<>(obj.ptr()->ob_type));
	return bp::extract<std::string>(bp::getattr(typ, "__name__"));
}

}

static bp::dict
array_interface(I3Matrix &m)
{
	bp::dict d;
	d["shape"]   = bp::make_tuple(m.size1(), m.size2());
	d["typestr"] = bp::str("<f8");
	d["data"]    = bp::make_tuple((long long)(&(m.data()[0])), false);
	
	return d;
}

static I3MatrixPtr
from_object(bp::object obj)
{
	pyerr_fmt err;
	if (!PyObject_HasAttrString(obj.ptr(), "__array_interface__")) {
		err << "'" << type(obj) << "' object does not support the array protocol!";
		err.raise(PyExc_TypeError);
	}
	
	bp::dict iface(bp::getattr(obj, "__array_interface__"));
	if (!(iface.has_key("shape") && iface.has_key("typestr") && iface.has_key("data"))) {
		err << "'" << type(obj) << "' object does not support the array protocol!";
		err.raise(PyExc_TypeError);
	}
	
	bp::tuple shape(iface["shape"]);
	if (bp::len(shape) != 2) {
		err << "Array must have 2 dimensions!";
		err.raise(PyExc_ValueError);
	}
	
	if (iface.has_key("strides") && iface["strides"]) {
		err << "I can't deal with strided arrays!";
		err.raise(PyExc_ValueError);
	}
	
	std::string typestr = bp::extract<std::string>(iface["typestr"]);
	if (typestr != "<f8") {
		err << "'" << type(obj) << "' object has unknown typecode: '" << typestr << "'";
		err.raise(PyExc_TypeError);
	}
	
	// NB: we ignore the read-only flag since we're making a copy anyway
	ptrdiff_t ptr = bp::extract<ptrdiff_t>(bp::tuple(iface["data"])[0]);
	size_t size1 = bp::extract<size_t>(shape[0]);
	size_t size2 = bp::extract<size_t>(shape[1]);
	
	I3MatrixPtr mat = boost::make_shared<I3Matrix>(size1, size2);
	memcpy(&(mat->data()[0]), (void*)ptr, size1*size2*sizeof(double));
	
	return mat;
}

void
register_I3Matrix()
{
	bp::class_<I3Matrix, boost::shared_ptr<I3Matrix>, bp::bases<I3FrameObject> >
	    ("I3Matrix", "I3Matrix has no dedicated python bindings. Use numpy.asarray(m) for data access.", bp::no_init)
	    .def(bp::init<size_t, size_t>("Create an unintialized matrix of the given size",
	        (bp::arg("size1"), "size2")))
	    .def(bp::init<size_t, size_t, const I3Matrix::value_type&>("Create and initialize a matrix of the given size",
	        (bp::arg("size1"), "size2", "value")))
	    .def("__init__", bp::make_constructor(from_object),
	        "Copy data from an object that supports the array protocol (e.g. a numpy array)")
	    .add_property("__array_interface__", array_interface)
	    .def(bp::dataclass_suite<I3Matrix>())
	;
	
	register_pointer_conversions<I3Matrix>();
}

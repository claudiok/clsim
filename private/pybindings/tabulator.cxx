
#include "clsim/tabulator/Axis.h"
#include "clsim/tabulator/Axes.h"

namespace bp = boost::python;

void register_Axis()
{
	using namespace clsim::tabulator;
	
	bp::class_<Axis, boost::shared_ptr<Axis>, boost::noncopyable>
	    ("Axis", bp::no_init)
	;
	
	bp::class_<LinearAxis, boost::shared_ptr<LinearAxis>, bp::bases<Axis> >
	    ("LinearAxis", bp::init<double,double,unsigned>((bp::arg("min"),"max","n_bins"),
	     "Create an axis with *n_bins* bins spaced linearly from *min* to *max*"))
	;
	
	bp::class_<PowerAxis, boost::shared_ptr<PowerAxis>, bp::bases<Axis> >
	    ("PowerAxis", bp::init<double,double,unsigned,unsigned>((bp::arg("min"),"max","n_bins",bp::arg("power")=1u),
	     "Create an axis with *n_bins* bins spaced from *min* to *max* in powers of *power*"))
	;
}

void register_Axes()
{
	using namespace clsim::tabulator;
	
	bp::class_<Axes, boost::shared_ptr<Axes>, boost::noncopyable>
	    ("Axes", bp::no_init)
	;
	
	// make it possible to construct an Axes instance with a list
	typedef std::vector<boost::shared_ptr<Axis> > AxisSeries;
	bp::class_<AxisSeries, boost::shared_ptr<AxisSeries> >("vector_AxisPtr");
	from_python_sequence<AxisSeries, variable_capacity_policy>(); 
	
	bp::class_<SphericalAxes, boost::shared_ptr<SphericalAxes>, bp::bases<Axes>, boost::noncopyable>
	    ("SphericalAxes", bp::init<const AxisSeries&>())
	;
	
	bp::class_<CylindricalAxes, boost::shared_ptr<CylindricalAxes>, bp::bases<Axes>, boost::noncopyable>
	    ("CylindricalAxes", bp::init<const AxisSeries&>())
	;
}

void register_tabulator()
{
	// Put all tabulator-related classes in a submodule
	bp::scope current;
	std::string submoduleName = bp::extract<std::string>(current.attr("__name__"));
	submoduleName.append(".tabulator");
	bp::object submodule(bp::borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("tabulator") = submodule;
	bp::scope submoduleScope = submodule;
	
	register_Axis();
	register_Axes();
}



#include "clsim/tabulator/Axis.h"
#include "clsim/I3CLSimHelperToFloatString.h"

#include <boost/foreach.hpp>

namespace clsim {

namespace tabulator {

Axis::Axis(double min, double max, unsigned n_bins)
    : min_(min), max_(max), n_bins_(n_bins)
{}

Axis::~Axis()
{}

std::string
Axis::GetIndexCode(const std::string &var) const {
	std::ostringstream ss;
	ss << "("<<var<<" >= "<<I3CLSimHelper::ToFloatString(max_)<<" ? "<<(n_bins_-1)<<"u"
	    <<" : "<<"("<<var<<" < "<<I3CLSimHelper::ToFloatString(min_)<<" ? 0u : (uint)("
	    << (n_bins_) << "*(" << GetInverseTransformCode(var) << " - " << GetInverseTransformCode(I3CLSimHelper::ToFloatString(min_)) << ")"
	    << "/(" << GetInverseTransformCode(I3CLSimHelper::ToFloatString(max_)) << " - " << GetInverseTransformCode(I3CLSimHelper::ToFloatString(min_)) << ")"
	    <<")))";
	;
	return ss.str();
}

std::vector<double>
Axis::GetBinEdges() const
{
	std::vector<double> edges(GetNBins()+1);

	double imin = InverseTransform(GetMin());
	double imax = InverseTransform(GetMax());
	double istep = (imax-imin)/GetNBins();
	for (unsigned i = 0; i <= GetNBins(); i++)
		edges[i] = Transform(imin + i*istep);

	return edges;
}

double
Axis::GetBinEdge(unsigned i) const
{
	double imin = InverseTransform(GetMin());
	double imax = InverseTransform(GetMax());
	double istep = (imax-imin)/GetNBins();
	return Transform(imin + i*istep);
}

LinearAxis::LinearAxis(double min, double max, unsigned n_bins)
    : Axis(min, max, n_bins)
{}

LinearAxis::~LinearAxis()
{}

double
LinearAxis::Transform(double value) const
{
	return value;
}

double
LinearAxis::InverseTransform(double value) const
{
	return value;
}

std::string
LinearAxis::GetTransformCode(const std::string &var) const
{
	return var;
}

std::string
LinearAxis::GetInverseTransformCode(const std::string &var) const
{
	return var;
}

PowerAxis::PowerAxis(double min, double max, unsigned n_bins, unsigned power)
    : Axis(min, max, n_bins), power_(power)
{}

PowerAxis::~PowerAxis()
{}

double
PowerAxis::Transform(double value) const
{
	return std::pow(value, power_);
}

double
PowerAxis::InverseTransform(double value) const
{
	return std::pow(value, 1./power_);
}

std::string
PowerAxis::GetTransformCode(const std::string &var) const
{
	std::ostringstream ss;
	if (power_ == 0) {
		ss << 1;
	} else if (power_ < 5) {
		ss << var;
		for (unsigned i = 0; i < power_-1; i++)
			ss << "*" << var;
	} else {
		ss << "pow("<<var<<", "<<power_<< ")";
	}
	return ss.str();
}

std::string
PowerAxis::GetInverseTransformCode(const std::string &var) const
{
	std::ostringstream ss;
	switch (power_) {
		case 0:
			ss << 1;
			break;
		case 1:
			ss << var;
			break;
		case 2:
			ss << "sqrt(" << var << ")";
			break;
		case 3:
			ss << "cbrt(" << var << ")";
			break;
		default:
			ss << "pow(" << var << "," << I3CLSimHelper::ToFloatString(1./power_) << ")";
	}
	return ss.str();
}

}

}

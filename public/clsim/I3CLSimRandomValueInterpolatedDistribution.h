#ifndef I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED
#define I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

#include <vector>

/**
 * @brief A random value chosen according to a given distribution.
 * The distribution is linearly interpolated between the given data
 * points.
 */
static const unsigned i3clsimrandomvalueinterpolateddistribution_version_ = 0;

struct I3CLSimRandomValueInterpolatedDistribution : public I3CLSimRandomValue
{
public:
    
    // arbitrary x values
    I3CLSimRandomValueInterpolatedDistribution(const std::vector<double> &x,
                                               const std::vector<double> &y);

    // x values with constant spacing (more efficient)
    I3CLSimRandomValueInterpolatedDistribution(double xFirst, double xSpacing,
                                               const std::vector<double> &y);

    virtual ~I3CLSimRandomValueInterpolatedDistribution();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    std::string WriteTableCode(const std::string &prefix) const;
    
    I3CLSimRandomValueInterpolatedDistribution();

    std::vector<double> x_;
    std::vector<double> y_;
    double constantXSpacing_;
    double firstX_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueInterpolatedDistribution, i3clsimrandomvalueinterpolateddistribution_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueInterpolatedDistribution);

#endif //I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED

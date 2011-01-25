#ifndef I3CLSIMRANDOMVALUETABULATEDDISTRIBUTIONCOSANGLE_H_INCLUDED
#define I3CLSIMRANDOMVALUETABULATEDDISTRIBUTIONCOSANGLE_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

#include <vector>

/**
 * @brief A value (cos(theta)) chosen according to a tabulated distribution
 * of p(theta). Note that the input is given as a function of the angle,
 * not as a function of cos(angle)!
 */
static const unsigned i3clsimrandomvaluetabulateddistributioncosangle_version_ = 0;

struct I3CLSimRandomValueTabulatedDistributionCosAngle : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueTabulatedDistributionCosAngle(const std::vector<double> &angles,
                                                    const std::vector<double> &values,
                                                    double powerLawIndexBeforeFirstBin);

    virtual ~I3CLSimRandomValueTabulatedDistributionCosAngle();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    std::string WriteScatteringTableCode(const std::string &prefix) const;
    
    I3CLSimRandomValueTabulatedDistributionCosAngle();

    std::vector<double> angles_;
    std::vector<double> values_;
    double powerLawIndexBeforeFirstBin_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueTabulatedDistributionCosAngle, i3clsimrandomvaluetabulateddistributioncosangle_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueTabulatedDistributionCosAngle);

#endif //I3CLSIMRANDOMVALUETABULATEDDISTRIBUTIONCOSANGLE_H_INCLUDED

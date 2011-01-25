#ifndef I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED
#define I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

/**
 * @brief The cosine of the scattering angle chosen according to
 * Rayleigh scattering.
 *
 */
static const unsigned i3clsimrandomvaluerayleighscatteringcosangle_version_ = 0;

struct I3CLSimRandomValueRayleighScatteringCosAngle : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueRayleighScatteringCosAngle();
    virtual ~I3CLSimRandomValueRayleighScatteringCosAngle();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueRayleighScatteringCosAngle, i3clsimrandomvaluerayleighscatteringcosangle_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueRayleighScatteringCosAngle);

#endif //I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED

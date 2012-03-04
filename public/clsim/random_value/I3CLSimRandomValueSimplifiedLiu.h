#ifndef I3CLSIMRANDOMVALUESIMPLIFIEDLIU_H_INCLUDED
#define I3CLSIMRANDOMVALUESIMPLIFIEDLIU_H_INCLUDED

#include "clsim/random_value/I3CLSimRandomValue.h"

/**
 * @brief A value chosen according to the simplified Liu distribution.
 * (see Pingyu Liu, "A new phase function approximating to mie
 * scattering for radiative transport equations", Phys. Med. Biol.,
 * 39:1025, 1994)
 *
 * The distribution looks like this:
 *
 * p(x) propto (1 + x)^alpha
 *
 * where x=cos(theta), alpha=2g/(1-g) and g=<cos(theta)> [the mean cosine].
 */
static const unsigned i3clsimrandomvaluesimplifiedliu_version_ = 0;

struct I3CLSimRandomValueSimplifiedLiu : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueSimplifiedLiu(double meanCosine);
    virtual ~I3CLSimRandomValueSimplifiedLiu();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueSimplifiedLiu();

    double meanCosine_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueSimplifiedLiu, i3clsimrandomvaluesimplifiedliu_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueSimplifiedLiu);

#endif //I3CLSIMRANDOMVALUESIMPLIFIEDLIU_H_INCLUDED

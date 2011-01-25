#ifndef I3CLSIMRANDOMVALUEHENYEYGREENSTEIN_H_INCLUDED
#define I3CLSIMRANDOMVALUEHENYEYGREENSTEIN_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

/**
 * @brief A value chosen according to a Henyey-Greenstein distribution.
 *
 * The distribution looks like this:
 *
 * p(x) = 0.5 * (1-g^2)/((1 + g^2 -2g*x)^(3/2))
 *
 * where x=cos(theta) and g=<cos(theta)> [the mean cosine].
 */
static const unsigned i3clsimrandomvaluehenyeygreenstein_version_ = 0;

struct I3CLSimRandomValueHenyeyGreenstein : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueHenyeyGreenstein(double meanCosine);
    virtual ~I3CLSimRandomValueHenyeyGreenstein();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueHenyeyGreenstein();

    double meanCosine_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueHenyeyGreenstein, i3clsimrandomvaluehenyeygreenstein_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueHenyeyGreenstein);

#endif //I3CLSIMRANDOMVALUEHENYEYGREENSTEIN_H_INCLUDED

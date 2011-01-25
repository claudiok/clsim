#ifndef I3CLSIMRANDOMVALUE_H_INCLUDED
#define I3CLSIMRANDOMVALUE_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include <string>

/**
 * @brief A value chosen from a random distribution
 */
static const unsigned i3clsimrandomvalue_version_ = 0;

struct I3CLSimRandomValue 
{
public:
    
    I3CLSimRandomValue();
    virtual ~I3CLSimRandomValue();

    /**
     * If the OpenCL function will only use a single random number,
     * the random number can be passed directly as a value instead of
     * passing an random number generator. This may improve performance,
     * so it can be made known to the caller.
     */
    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const = 0;
    
    /**
     * Shall return an OpenCL-compatible function.
     * The declaration of the form "float {functionName}({functionArgs})"
     * is provided by the caller in functionDecl.
     * The function call to generate uniformly distributed
     * random numbers between 0 and 1 is
     * provided in uniformRandomCall_{co|oc}.
     * (co: closed-open, 0 included, 1 not included
     *  oc: open-closed, 0 not included, 1 included)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const = 0;

    /**
     * Shall compare the internal state another I3CLSimRandomValue object
     */
    virtual bool CompareTo(const I3CLSimRandomValue &other) const = 0;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimRandomValue& a, const I3CLSimRandomValue& b)
{
    return a.CompareTo(b);
}


BOOST_CLASS_VERSION(I3CLSimRandomValue, i3clsimrandomvalue_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValue);

#endif //I3CLSIMRANDOMVALUE_H_INCLUDED

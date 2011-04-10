#ifndef I3CLSIMRANDOMVALUEAPPLYFUNCTION_H_INCLUDED
#define I3CLSIMRANDOMVALUEAPPLYFUNCTION_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

#include <vector>
#include <string>

/**
 * @brief Chooses a random value according to a specified distribution
 * and applies a configured function to the return value.
 */
static const unsigned i3clsimrandomvalueapplyfunction_version_ = 0;

struct I3CLSimRandomValueApplyFunction : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueApplyFunction(I3CLSimRandomValuePtr randomDistUsed,
                                    const std::string &functionName);

    virtual ~I3CLSimRandomValueApplyFunction();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return randomDistUsed_->OpenCLFunctionWillOnlyUseASingleRandomNumber();}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueApplyFunction();

    I3CLSimRandomValuePtr randomDistUsed_;
    std::string applyFunctionName_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimRandomValueApplyFunction, i3clsimrandomvalueapplyfunction_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueApplyFunction);

#endif //I3CLSIMRANDOMVALUEAPPLYFUNCTION_H_INCLUDED

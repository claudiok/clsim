#ifndef I3CLSIMRANDOMVALUEMIXED_H_INCLUDED
#define I3CLSIMRANDOMVALUEMIXED_H_INCLUDED

#include "clsim/I3CLSimRandomValue.h"

/**
 * @brief A mix of two random values.
 */
static const unsigned i3clsimrandomvaluemixed_version_ = 0;

struct I3CLSimRandomValueMixed : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueMixed(double fractionOfFirstDistribution,
                            I3CLSimRandomValueConstPtr firstDistribution,
                            I3CLSimRandomValueConstPtr secondDistribution);
    virtual ~I3CLSimRandomValueMixed();

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const;

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueMixed();

    double fractionOfFirstDistribution_;
    I3CLSimRandomValueConstPtr firstDistribution_;
    I3CLSimRandomValueConstPtr secondDistribution_;
    
    friend class boost::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


BOOST_CLASS_VERSION(I3CLSimRandomValueMixed, i3clsimrandomvaluemixed_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueMixed);

#endif //I3CLSIMRANDOMVALUEMIXED_H_INCLUDED

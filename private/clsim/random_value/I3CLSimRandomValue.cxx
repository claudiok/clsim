#include <icetray/serialization.h>
#include <clsim/random_value/I3CLSimRandomValue.h>

I3CLSimRandomValue::I3CLSimRandomValue()
{ 
    
}

I3CLSimRandomValue::~I3CLSimRandomValue() 
{ 

}

template <class Archive>
void I3CLSimRandomValue::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvalue_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValue class.",version,i3clsimrandomvalue_version_);
}     


I3_SERIALIZABLE(I3CLSimRandomValue);

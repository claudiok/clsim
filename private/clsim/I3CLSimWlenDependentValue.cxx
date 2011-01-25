#include <icetray/serialization.h>
#include <clsim/I3CLSimWlenDependentValue.h>

I3CLSimWlenDependentValue::I3CLSimWlenDependentValue()
{ 
    
}

I3CLSimWlenDependentValue::~I3CLSimWlenDependentValue() 
{ 

}

template <class Archive>
void I3CLSimWlenDependentValue::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvalue_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValue class.",version,i3clsimwlendependentvalue_version_);
}     


I3_SERIALIZABLE(I3CLSimWlenDependentValue);

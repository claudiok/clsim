#include <icetray/serialization.h>
#include <clsim/function/I3CLSimFunction.h>

I3CLSimFunction::I3CLSimFunction()
{ 
    
}

I3CLSimFunction::~I3CLSimFunction() 
{ 

}

template <class Archive>
void I3CLSimFunction::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunction_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunction class.",version,i3clsimfunction_version_);
}     


I3_SERIALIZABLE(I3CLSimFunction);

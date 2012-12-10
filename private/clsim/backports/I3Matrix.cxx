
#include <icetray/serialization.h>
#include "clsim/backports/I3Matrix.h"

template <typename Archive>
void
I3Matrix::serialize(Archive &ar, unsigned int version)
{	
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Matrix", base_object<I3Matrix::base>(*this));
}

I3_SERIALIZABLE(I3Matrix);

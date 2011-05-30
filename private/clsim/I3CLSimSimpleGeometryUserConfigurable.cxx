#include "clsim/I3CLSimSimpleGeometryUserConfigurable.h"


I3CLSimSimpleGeometryUserConfigurable::
I3CLSimSimpleGeometryUserConfigurable(double OMRadius,
                                      std::size_t numOMs)
:
OMRadius_(OMRadius),
numOMs_(numOMs),
stringIDs_(numOMs, 0),
domIDs_(numOMs, 0),
posX_(numOMs, NAN),
posY_(numOMs, NAN),
posZ_(numOMs, NAN),
subdetectors_(numOMs, "")
{
    
}

I3CLSimSimpleGeometryUserConfigurable::
~I3CLSimSimpleGeometryUserConfigurable()
{
    
}

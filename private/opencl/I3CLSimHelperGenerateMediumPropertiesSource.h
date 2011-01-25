#ifndef I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_H_INCLUDED
#define I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_H_INCLUDED

#include <string>

#include "clsim/I3CLSimMediumProperties.h"

namespace I3CLSimHelper
{
    /**
     * generates the OpenCL source code for a given mediumProperties object.
     */
    std::string GenerateMediumPropertiesSource(const I3CLSimMediumProperties &mediumProperties);

};

#endif //I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_H_INCLUDED

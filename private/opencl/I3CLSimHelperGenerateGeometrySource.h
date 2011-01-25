#ifndef I3CLSIMHELPERGENERATEGEOMETRYSOURCE_H_INCLUDED
#define I3CLSIMHELPERGENERATEGEOMETRYSOURCE_H_INCLUDED

#include <string>
#include <vector>

#include "clsim/I3CLSimSimpleGeometry.h"

namespace I3CLSimHelper
{
    /**
     * generates the OpenCL source code for a given I3CLSimSimpleGeometry object.
     */
    std::string GenerateGeometrySource(const I3CLSimSimpleGeometry &geometry,
                                       std::vector<unsigned char> &geoLayerToOMNumIndexPerStringSetBuffer);

};

#endif //I3CLSIMHELPERGENERATEGEOMETRYSOURCE_H_INCLUDED

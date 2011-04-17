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
                                       std::vector<unsigned char> &geoLayerToOMNumIndexPerStringSetBuffer,
                                       std::vector<int> &stringIndexToStringIDBuffer,
                                       std::vector<std::vector<unsigned int> > &domIndexToDomIDBuffer_perStringIndex);

};

#endif //I3CLSIMHELPERGENERATEGEOMETRYSOURCE_H_INCLUDED

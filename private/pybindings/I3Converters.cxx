//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <sstream>

#include "clsim/converter/I3PhotonConverter.h"
#include "clsim/converter/I3MCHitConverterWithIDs.h"
#include "tableio/converter/pybindings.h"

void register_I3Converters()
{
    I3CONVERTER_NAMESPACE(clsim);

    I3CONVERTER_EXPORT_DEFAULT(I3PhotonConverter, "Dumps a single I3Photon to a table column");
    I3_MAP_CONVERTER_EXPORT_DEFAULT(I3PhotonSeriesMapConverter,"Dumps all I3Photons in a I3PhotonSeriesMap");

    // this is not the default converter, the default one is in dataclasses
    I3_MAP_CONVERTER_EXPORT(I3MCHitSeriesMapConverterWithIDs, "Dumps all I3MCHits from a I3MCHitSeriesMap");

}

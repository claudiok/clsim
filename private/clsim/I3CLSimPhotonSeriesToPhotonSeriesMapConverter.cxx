#include "clsim/I3CLSimPhotonSeriesToPhotonSeriesMapConverter.h"

#include <boost/foreach.hpp>


I3CLSimPhotonSeriesToPhotonSeriesMapConverter::
I3CLSimPhotonSeriesToPhotonSeriesMapConverter()
{
    
}

I3CLSimPhotonSeriesToPhotonSeriesMapConverter::
~I3CLSimPhotonSeriesToPhotonSeriesMapConverter()
{
    
}

void
I3CLSimPhotonSeriesToPhotonSeriesMapConverter::
SetPhotonSeries(I3CLSimPhotonSeriesConstPtr photonSeries)
{
    currentPhotonSeries_ = photonSeries;

}

void
I3CLSimPhotonSeriesToPhotonSeriesMapConverter::
SetGeometry(I3CLSimSimpleGeometryFromI3GeometryConstPtr geometry)
{
    currentGeometry_ = geometry;
}

I3CLSimPhotonSeriesMapPtr
I3CLSimPhotonSeriesToPhotonSeriesMapConverter::
GetPhotonSeriesMapForIdentifier(uint32_t identifier)
{
    if (!currentPhotonSeries_)
        throw I3CLSimPhotonSeriesToPhotonSeriesMapConverter_exception
        ("no I3CLSimPhotonSeries set before call to GetPhotonSeriesMapForIdentifier().");

    if (!currentGeometry_)
        throw I3CLSimPhotonSeriesToPhotonSeriesMapConverter_exception
        ("no I3CLSimSimpleGeometryFromI3Geometry set before call to GetPhotonSeriesMapForIdentifier().");

    // return value
    I3CLSimPhotonSeriesMapPtr retMap = I3CLSimPhotonSeriesMapPtr(new I3CLSimPhotonSeriesMap());
    
    BOOST_FOREACH(const I3CLSimPhoton &photon, *currentPhotonSeries_)
    {
        if (photon.identifier != identifier) continue; // skip photons with different IDs
        int32_t stringAndDomID = photon.dummy;

        int32_t stringID;
        uint32_t domID;
        
        if (stringAndDomID >= 0)
        {
            stringID = stringAndDomID / 1000;
            domID = stringAndDomID % 1000;
        }
        else
        {
            // convention for negative string IDs
            stringID = -((-stringAndDomID) / 1000);
            domID = (-stringAndDomID) % 1000;
        }
    
        // now the OMKey can be generated
        const OMKey key(stringID, domID);
        
        // this either inserts a new vector or retrieves an existing one
        I3CLSimPhotonSeries &outputPhotonSeries = retMap->insert(std::make_pair(key, I3CLSimPhotonSeries())).first->second;
        
        // append a copy of this photon to the output vector
        outputPhotonSeries.push_back(photon);
    }
    
    return retMap;
}

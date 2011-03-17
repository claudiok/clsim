#ifndef I3CLSIMPHOTONSERIESTOPHOTONSERIESMAPCONVERTER_H_INCLUDED
#define I3CLSIMPHOTONSERIESTOPHOTONSERIESMAPCONVERTER_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include "dataclasses/geometry/I3Geometry.h"

#include "clsim/I3CLSimSimpleGeometry.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimPhoton.h"

#include <boost/noncopyable.hpp>

#include <map>
#include <string>
#include <stdexcept>

/**
 * @brief Simple helper class that takes a vector of I3CLSimPhotons
 * and uses their "dummy" field to get reconstruct OMKey for each photon.
 * The initial photon vector may contain data from different events
 * which are distinguished by the I3CLSimPhoton "identifier" field.
 */

class I3CLSimPhotonSeriesToPhotonSeriesMapConverter_exception : public std::runtime_error
{
public:
    I3CLSimPhotonSeriesToPhotonSeriesMapConverter_exception(const std::string &msg)
    :std::runtime_error(msg)
    {;}
};

struct I3CLSimPhotonSeriesToPhotonSeriesMapConverter : private boost::noncopyable
{
public:
    I3CLSimPhotonSeriesToPhotonSeriesMapConverter();
    ~I3CLSimPhotonSeriesToPhotonSeriesMapConverter();

    /**
     * Sets the current photon series. All output values will be re-set.
     */
    void SetPhotonSeries(I3CLSimPhotonSeriesConstPtr photonSeries);

    /**
     * Sets the geometry.
     * All output values will be re-set.
     */
    void SetGeometry(I3CLSimSimpleGeometryFromI3GeometryConstPtr geometry);
    
    /**
     * Retrieves a I3CLSimPhotonSeriesMap containing photons with a
     * given identifier. The map may be empty.
     * Will throw if geometry or input photon series are not set.
     */
    I3CLSimPhotonSeriesMapPtr GetPhotonSeriesMapForIdentifier(uint32_t identifier);
    
private:
    I3CLSimPhotonSeriesConstPtr currentPhotonSeries_;
    I3CLSimSimpleGeometryFromI3GeometryConstPtr currentGeometry_;
};

I3_POINTER_TYPEDEFS(I3CLSimPhotonSeriesToPhotonSeriesMapConverter);

#endif //I3CLSIMPHOTONSERIESTOPHOTONSERIESMAPCONVERTER_H_INCLUDED

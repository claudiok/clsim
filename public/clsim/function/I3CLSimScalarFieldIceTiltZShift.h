/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimScalarFieldIceTiltZShift.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSCALARFIELDICETILTZSHIFT_H_INCLUDED
#define I3CLSIMSCALARFIELDICETILTZSHIFT_H_INCLUDED

#include "clsim/function/I3CLSimScalarField.h"

#include <vector>

#include "dataclasses/I3Matrix.h"

/**
 * @brief Returns a local z-correction for the ice layer positions
 * as a function of (x,y,z) for ice tilt handling. Implements the
 * the PPC algorithm.
 */
static const unsigned i3clsimscalarfieldicetiltzshift_version_ = 0;

struct I3CLSimScalarFieldIceTiltZShift : public I3CLSimScalarField
{
public:
    static const double default_directionOfTiltAzimuth;

    I3CLSimScalarFieldIceTiltZShift(
        const std::vector<double> &distancesFromOriginAlongTilt,
        const std::vector<double> &zCoordinates,
        const I3Matrix &zCorrections,
        double directionOfTiltAzimuth=default_directionOfTiltAzimuth);

    virtual ~I3CLSimScalarFieldIceTiltZShift();
    
    /**
     * if this is true, it is assumed that GetValue() returns a
     * meaningful value. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const;
    
    /**
     * return the value at a requested 3-vector
     */
    virtual double GetValue(double x, double y, double z) const;

    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;

    /**
     * Shall compare to another I3CLSimScalarField object
     */
    virtual bool CompareTo(const I3CLSimScalarField &other) const;
    
private:
    I3CLSimScalarFieldIceTiltZShift();

    std::vector<double> distancesFromOriginAlongTilt_;
    std::vector<double> zCoordinates_;
    I3Matrix zCorrections_;
    double directionOfTiltAzimuth_;

    double firstZCoordinate_;
    double zCoordinateSpacing_;

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimScalarFieldIceTiltZShift, i3clsimscalarfieldicetiltzshift_version_);

I3_POINTER_TYPEDEFS(I3CLSimScalarFieldIceTiltZShift);

#endif //I3CLSIMSCALARFIELDICETILTZSHIFT_H_INCLUDED

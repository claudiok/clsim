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
 * @file I3CLSimFunctionScatLenPartic.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED
#define I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <limits>

/**
 * @brief The scattering length as defined by the "partic"
 * scattering model.
 */
static const unsigned i3clsimfunctionscatlenpartic_version_ = 0;

struct I3CLSimFunctionScatLenPartic : public I3CLSimFunction
{
public:
    static const double default_volumeConcentrationSmallParticles;
    static const double default_volumeConcentrationLargeParticles;
    
    I3CLSimFunctionScatLenPartic(double volumeConcentrationSmallParticles=default_volumeConcentrationSmallParticles,   // fraction (e.g. 0.0075*I3Units::perMillion)
                                 double volumeConcentrationLargeParticles=default_volumeConcentrationLargeParticles    // fraction (e.g. 0.0075*I3Units::perMillion)
                                );
    virtual ~I3CLSimFunctionScatLenPartic();
    
    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const {return true;};
    
    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const {return false;};
    
    /**
     * Shall return the value at a requested wavelength (n)
     */
    virtual double GetValue(double wlen) const;
    
    /**
     * Shall return the minimal supported wavelength (possibly -inf)
     */
    virtual double GetMinWlen() const {return -std::numeric_limits<double>::infinity();}
    
    /**
     * Shall return the maximal supported wavelength (possibly +inf)
     */
    virtual double GetMaxWlen() const {return std::numeric_limits<double>::infinity();}
    
    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;
    
    /**
     * Shall compare to another I3CLSimFunction object
     */
    virtual bool CompareTo(const I3CLSimFunction &other) const;
    
private:
    double volumeConcentrationSmallParticles_;
    double volumeConcentrationLargeParticles_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionScatLenPartic, i3clsimfunctionscatlenpartic_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionScatLenPartic);

#endif //I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED

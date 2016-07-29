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
 * @file I3CLSimFunctionRefIndexIceCube.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONREFINDEXICECUBE_H_INCLUDED
#define I3CLSIMFUNCTIONREFINDEXICECUBE_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <limits>
#include <string>

/**
 * @brief The phase refractive index for IceCube glacial ice, taken from
 * "Role of Group and Phase Velocity in High-Energy Neutrino Observatories",
 * P.B. Price and K. Woschnagg, Astropart. Phys. 15 (2001) 97
 */
static const unsigned i3clsimfunctionrefindexicecube_version_ = 0;

struct I3CLSimFunctionRefIndexIceCube : public I3CLSimFunction
{
public:
    static const std::string default_mode;
    static const double default_n0;
    static const double default_n1;
    static const double default_n2;
    static const double default_n3;
    static const double default_n4;
    static const double default_g0;
    static const double default_g1;
    static const double default_g2;
    static const double default_g3;
    static const double default_g4;
    
    
    I3CLSimFunctionRefIndexIceCube(std::string mode = default_mode,
                                             double n0 = default_n0,  // coefficients 0-4 (for phase ref index)
                                             double n1 = default_n1,
                                             double n2 = default_n2,
                                             double n3 = default_n3,
                                             double n4 = default_n4,
                                             double g0 = default_g0,  // coefficients 0-4 (for group ref index)
                                             double g1 = default_g1,
                                             double g2 = default_g2,
                                             double g3 = default_g3,
                                             double g4 = default_g4
                                             );
    virtual ~I3CLSimFunctionRefIndexIceCube();
    
    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const {return true;};
    
    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const {return true;};
    
    /**
     * Shall return the value at a requested wavelength (n)
     */
    virtual double GetValue(double wlen) const;
    
    /**
     * Shall return the derivative at a requested wavelength (dn/dlambda)
     */
    virtual double GetDerivative(double wlen) const;
    
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
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const;
    
    /**
     * Shall compare to another I3CLSimFunction object
     */
    virtual bool CompareTo(const I3CLSimFunction &other) const;
    
private:
    std::string mode_;
    double n0_;
    double n1_;
    double n2_;
    double n3_;
    double n4_;
    double g0_;
    double g1_;
    double g2_;
    double g3_;
    double g4_;
    
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionRefIndexIceCube, i3clsimfunctionrefindexicecube_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionRefIndexIceCube);

#endif //I3CLSIMFUNCTIONREFINDEXICECUBE_H_INCLUDED

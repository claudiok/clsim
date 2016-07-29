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
 * @file I3CLSimFunctionPolynomial.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONPOLYNOMIAL_H_INCLUDED
#define I3CLSIMFUNCTIONPOLYNOMIAL_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <vector>

/**
 * @brief A simple polynomial
 */
static const unsigned i3clsimfunctionpolynomial_version_ = 0;

struct I3CLSimFunctionPolynomial : public I3CLSimFunction
{
public:
    /**
     * Constructor for polynomials that are defined for the whole x-range
     **/
    I3CLSimFunctionPolynomial(const std::vector<double> &coeffs);
    
    /**
     * Constructor for polynomials that have a limited range
     * If it is called with an x outside the range, it will return the
     * values at rangemin/rangemax.
     */
    I3CLSimFunctionPolynomial(const std::vector<double> &coeffs,
                                        double rangemin,
                                        double rangemax);
    
    /**
     * Constructor for polynomials that have a limited range
     * If it is called with an x outside the range,
     * it will return the underflow/overflow value
     */
    I3CLSimFunctionPolynomial(const std::vector<double> &coeffs,
                                        double rangemin, double rangemax,
                                        double underflow, double overflow);
    
    /**
     * Destructor
     */
    virtual ~I3CLSimFunctionPolynomial();
    
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
    
    
    // access to the internal state
    inline const std::vector<double> &GetCoefficients() const {return coefficients_;};
    
private:
    I3CLSimFunctionPolynomial();
    
    std::vector<double> coefficients_;
    
    double rangemin_, rangemax_; //The range for which the polynomial can be calculated
    double underflow_, overflow_; //Value that is returned if the polynomial is called with an argument outside the range
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionPolynomial, i3clsimfunctionpolynomial_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionPolynomial);

#endif //I3CLSIMFUNCTIONPOLYNOMIAL_H_INCLUDED

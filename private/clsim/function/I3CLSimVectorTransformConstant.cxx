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
 * @file I3CLSimVectorTransformConstant.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimVectorTransformConstant.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <sstream>
#include <stdexcept>

I3CLSimVectorTransformConstant::
I3CLSimVectorTransformConstant()
{
}

I3CLSimVectorTransformConstant::~I3CLSimVectorTransformConstant() 
{;}

bool I3CLSimVectorTransformConstant::HasNativeImplementation() const 
{
    return true;
}

std::vector<double> I3CLSimVectorTransformConstant::ApplyTransform(const std::vector<double> &vec) const
{
    if (vec.size() != 3)
        throw std::range_error("vector must contain excatly 3 elements!");

    // return a copy
    return vec;
}

std::string I3CLSimVectorTransformConstant::GetOpenCLFunction(const std::string &functionName) const
{
    // the OpenCL interface takes a pointer to a float4, but ignores the fourth component

    std::string funcDef = 
    std::string("inline void ") + functionName + std::string("(float4 *vec)");
    
    std::string funcBody = std::string() + 
    "{\n"
    "    return; // does nothing\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + "\n" + funcBody;
}

bool I3CLSimVectorTransformConstant::CompareTo(const I3CLSimVectorTransform &other) const
{
    try
    {
        // this would be used if this class had any members to compare
        //const I3CLSimVectorTransformConstant &other_ = dynamic_cast<const I3CLSimVectorTransformConstant &>(other);
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}



template <class Archive>
void I3CLSimVectorTransformConstant::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimvectortransformconstant_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimVectorTransformConstant class.",version,i3clsimvectortransformconstant_version_);

    ar & make_nvp("I3CLSimVectorTransform", base_object<I3CLSimVectorTransform>(*this));
}     


I3_SERIALIZABLE(I3CLSimVectorTransformConstant);

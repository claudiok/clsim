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
 * @file I3CLSimVectorTransformMatrix.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimVectorTransformMatrix.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <sstream>
#include <stdexcept>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimVectorTransformMatrix::
I3CLSimVectorTransformMatrix(
    const I3Matrix &matrix,
    bool renormalize)
:
matrix_(matrix),
renormalize_(renormalize)
{
    if ((matrix_.size1() != 3) || (matrix_.size2() != 3))
        throw std::range_error("matrix must be 3x3!");
}

I3CLSimVectorTransformMatrix::I3CLSimVectorTransformMatrix() {;}

I3CLSimVectorTransformMatrix::~I3CLSimVectorTransformMatrix() 
{;}

bool I3CLSimVectorTransformMatrix::HasNativeImplementation() const 
{
    return true;
}

std::vector<double> I3CLSimVectorTransformMatrix::ApplyTransform(const std::vector<double> &vec) const
{
    if (vec.size() != 3)
        throw std::range_error("vector must contain excatly 3 elements!");

    boost::numeric::ublas::vector<double> uvec_in(3);
    std::copy(vec.begin(), vec.end(), uvec_in.begin());
    boost::numeric::ublas::vector<double> uvec_out(3);

    boost::numeric::ublas::noalias(uvec_out) = boost::numeric::ublas::prod(matrix_, uvec_in);

    std::vector<double> out_vec(3, NAN);

    double norm=1.;

    if (renormalize_) {
        norm=0.;
        for (std::size_t i=0;i<3;++i)
        {
            norm += (uvec_out(i)*uvec_out(i));
        }
        norm = std::sqrt(norm);

        for (std::size_t i=0;i<3;++i)
        {
            out_vec[i] = uvec_out(i)/norm;
        }
    } else {
        for (std::size_t i=0;i<3;++i)
        {
            out_vec[i] = uvec_out(i);
        }
    }

    return out_vec;
}

std::string I3CLSimVectorTransformMatrix::GetOpenCLFunction(const std::string &functionName) const
{
    // the OpenCL interface takes a pointer to a float4, but ignores the fourth component

    std::string funcDef = 
    std::string("inline void ") + functionName + std::string("(float4 *vec)");
    
    std::string funcBody = std::string() + 
    "{\n";

    funcBody = funcBody +
    "    *vec = (float4)\n"
    "    (\n"
    "        (" + ToFloatString(matrix_(0,0)) + "*(*vec).x)+(" + ToFloatString(matrix_(0,1)) + "*(*vec).y)+(" + ToFloatString(matrix_(0,2)) + "*(*vec).z),\n"
    "        (" + ToFloatString(matrix_(1,0)) + "*(*vec).x)+(" + ToFloatString(matrix_(1,1)) + "*(*vec).y)+(" + ToFloatString(matrix_(1,2)) + "*(*vec).z),\n"
    "        (" + ToFloatString(matrix_(2,0)) + "*(*vec).x)+(" + ToFloatString(matrix_(2,1)) + "*(*vec).y)+(" + ToFloatString(matrix_(2,2)) + "*(*vec).z),\n"
    "        (*vec).w\n"
    "    );\n";

    if (renormalize_) {
        funcBody = funcBody +
        "#ifdef USE_NATIVE_MATH\n"
        "    (*vec).xyz = fast_normalize((*vec).xyz);\n"
        "#else\n"
        "    const float norm = rsqrt((*vec).x*(*vec).x + (*vec).y*(*vec).y + (*vec).z*(*vec).z);\n"
        "    (*vec).xyz = (*vec).xyz*norm;\n"
        "#endif\n";
    }

    funcBody = funcBody +
    "}\n";
    
    return funcDef + ";\n\n" + funcDef + "\n" + funcBody;
}

bool I3CLSimVectorTransformMatrix::CompareTo(const I3CLSimVectorTransform &other) const
{
    try
    {
        const I3CLSimVectorTransformMatrix &other_ = dynamic_cast<const I3CLSimVectorTransformMatrix &>(other);

        if ((other_.matrix_.size1() != matrix_.size1()))
            return false;
        if ((other_.matrix_.size2() != matrix_.size2()))
            return false;
        for (std::size_t i=0;i<matrix_.size1();++i)
            for (std::size_t j=0;j<matrix_.size2();++j)
            {
                if ((other_.matrix_(i,j) != matrix_(i,j)))
                    return false;
            }

        if ((other_.renormalize_ != renormalize_))
            return false;

        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}



template <class Archive>
void I3CLSimVectorTransformMatrix::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimvectortransformmatrix_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimVectorTransformMatrix class.",version,i3clsimvectortransformmatrix_version_);

    ar & make_nvp("I3CLSimVectorTransform", base_object<I3CLSimVectorTransform>(*this));
    ar & make_nvp("matrix", matrix_);
    ar & make_nvp("renormalize", renormalize_);
}     


I3_SERIALIZABLE(I3CLSimVectorTransformMatrix);

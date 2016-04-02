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
 * @file I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.h"

#include <sstream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <boost/preprocessor.hpp>

#include "icetray/I3Units.h"

#include "clsim/function/I3CLSimFunctionAbsLenIceCube.h"
#include "clsim/function/I3CLSimFunctionScatLenIceCube.h"

#include "clsim/I3CLSimHelperToFloatString.h"


namespace I3CLSimHelper
{
    
    
    
    // takes a vector of boost::shared_ptrs-to-const and dynamic_casts all
    // contents to a new type. only returns a filled vector if the
    // cast succeeds for all elements. "worked" is set to true
    // if the cast worked and "false" if it did not.
    template<class T, class U>
    std::vector<boost::shared_ptr<const T> > dynamic_vector_contents_cast(const std::vector<boost::shared_ptr<const U> > &inVec, bool &worked)
    {
        std::vector<boost::shared_ptr<const T> > retVec;
        
        BOOST_FOREACH(const boost::shared_ptr<const U> &inputValue, inVec)
        {
            if (!inputValue) {
                retVec.push_back(boost::shared_ptr<const T>());
                continue;
            }
            
            boost::shared_ptr<const T> convertedPtr = boost::dynamic_pointer_cast<const T>(inputValue);
            if (!convertedPtr) {worked=false; return std::vector<boost::shared_ptr<const T> >();}
            
            retVec.push_back(convertedPtr);
        }
        
        worked=true;
        return retVec;
    }
    
    template<typename T, class V>
    bool is_constant(const std::vector<T> &vect, const V &visitor)
    {
        if (vect.size() <= 1) return true;
        
        for (std::size_t i = 1; i<vect.size(); ++i)
        {
            if (!(visitor(vect[0]) == visitor(vect[i]))) return false;
        }
        
        return true;
    }
    
    
    namespace {
        #define CHECK_THESE_FOR_CONSTNESS_ABSLEN                        \
        (Kappa)(A)(B)(D)(E)

        #define CHECK_THESE_FOR_CONSTNESS_SCATLEN                       \
        (Alpha)

        #define GEN_VISITOR(r, data, a)                                 \
        struct BOOST_PP_CAT(visit_,a)                                   \
        {                                                               \
            template<typename T>                                        \
            double operator()(const boost::shared_ptr<T> &ptr) const           \
            {                                                           \
                if (!ptr) return NAN;                                   \
                return BOOST_PP_CAT(ptr->Get,a)();                      \
            }                                                           \
                                                                        \
            template<typename T>                                        \
            double operator()(const boost::shared_ptr<const T> &ptr) const     \
            {                                                           \
                if (!ptr) return NAN;                                   \
                return BOOST_PP_CAT(ptr->Get,a)();                      \
            }                                                           \
        };

        BOOST_PP_SEQ_FOR_EACH(GEN_VISITOR, ~, CHECK_THESE_FOR_CONSTNESS_ABSLEN)
        BOOST_PP_SEQ_FOR_EACH(GEN_VISITOR, ~, CHECK_THESE_FOR_CONSTNESS_SCATLEN)
    }
    
    
    // optimizers (special converters for functions
    // where we can generate more optimized code
    // for layered values)
    std::string GenerateOptimizedCodeFor_I3CLSimFunctionAbsLenIceCube(const std::vector<I3CLSimFunctionConstPtr> &layeredFunction,
                                                                                const std::string &fullName,
                                                                                const std::string &functionName,
                                                                                bool &worked)
    {
        // are all of them of type I3CLSimFunctionAbsLenIceCube?
        const std::vector<I3CLSimFunctionAbsLenIceCubeConstPtr> layeredFunctionIceCube =
        dynamic_vector_contents_cast<I3CLSimFunctionAbsLenIceCube>(layeredFunction, worked);
        if (layeredFunctionIceCube.size() <= 1) worked=false;
        if (!worked) return std::string("");

        {
            // this currently only works if kappa, A, B, D and E are constant
            bool isConst=true;
            #define CONST_CHECK(r, data, a) isConst = isConst && is_constant(layeredFunctionIceCube, BOOST_PP_CAT(visit_,a)());
            BOOST_PP_SEQ_FOR_EACH(CONST_CHECK, ~, CHECK_THESE_FOR_CONSTNESS_ABSLEN)
            #undef CONST_CHECK
            if (!isConst) {worked=false; return "";}
        }        
        
        std::ostringstream code;

        code << "///////////////// START " << fullName << " (optimized) ////////////////\n";
        code << "\n";

        code << "__constant float " << functionName << "_aDust400[" << layeredFunctionIceCube.size() << "] = {\n";
        BOOST_FOREACH(const I3CLSimFunctionAbsLenIceCubeConstPtr &function, layeredFunctionIceCube)
        {
            code << "    " << ToFloatString(function->GetADust400()) << ",\n";
        }
        code << "};\n";
        code << "\n";

        code << "__constant float " << functionName << "_deltaTau[" << layeredFunctionIceCube.size() << "] = {\n";
        BOOST_FOREACH(const I3CLSimFunctionAbsLenIceCubeConstPtr &function, layeredFunctionIceCube)
        {
            code << "    " << ToFloatString(function->GetDeltaTau()) << ",\n";
        }
        code << "};\n";
        code << "\n";

        code << "inline float " << functionName << "(unsigned int layer, float wlen);\n\n";
        code << "inline float " << functionName << "(unsigned int layer, float wlen)\n";
        code << "{\n";
        
        code << "    const float kappa = " << ToFloatString(layeredFunctionIceCube[0]->GetKappa()) << ";\n";
        code << "    const float A = " << ToFloatString(layeredFunctionIceCube[0]->GetA()) << ";\n";
        code << "    const float B = " << ToFloatString(layeredFunctionIceCube[0]->GetB()) << ";\n";
        code << "    const float D = " << ToFloatString(layeredFunctionIceCube[0]->GetD()) << ";\n";
        code << "    const float E = " << ToFloatString(layeredFunctionIceCube[0]->GetE()) << ";\n";
        code << "    \n";
        code << "    const float x = wlen/" << ToFloatString(I3Units::nanometer) << ";\n";
        code << "    \n";
        code << "#ifdef USE_NATIVE_MATH\n";
        code << "    return " << ToFloatString(I3Units::meter) << "*native_recip( (D*" << functionName << "_aDust400[layer]+E) * native_powr(x, -kappa)  +  A*native_exp(-B/x) * (1.f + 0.01f*" << functionName << "_deltaTau[layer]) );\n";
        code << "#else\n";
        code << "    return " << ToFloatString(I3Units::meter) << "/( (D*" << functionName << "_aDust400[layer]+E) * powr(x, -kappa)  +  A*exp(-B/x) * (1.f + 0.01f*" << functionName << "_deltaTau[layer]) );\n";
        code << "#endif\n";
        
        code << "}\n";
        code << "\n";
        
        code << "///////////////// END " << fullName << " (optimized) ////////////////\n";
        code << "\n";        

        worked=true;
        return code.str();
    }
    
    
    
    
    std::string GenerateOptimizedCodeFor_I3CLSimFunctionScatLenIceCube(const std::vector<I3CLSimFunctionConstPtr> &layeredFunction,
                                                                                 const std::string &fullName,
                                                                                 const std::string &functionName,
                                                                                 bool &worked)
    {
        // are all of them of type I3CLSimFunctionScatLenIceCube?
        const std::vector<I3CLSimFunctionScatLenIceCubeConstPtr> layeredFunctionIceCube =
        dynamic_vector_contents_cast<I3CLSimFunctionScatLenIceCube>(layeredFunction, worked);
        if (layeredFunctionIceCube.size() <= 1) worked=false;
        if (!worked) return std::string("");
        
        {
            // this currently only works if alpha is constant
            bool isConst=true;
#define CONST_CHECK(r, data, a) isConst = isConst && is_constant(layeredFunctionIceCube, BOOST_PP_CAT(visit_,a)());
            BOOST_PP_SEQ_FOR_EACH(CONST_CHECK, ~, CHECK_THESE_FOR_CONSTNESS_SCATLEN)
#undef CONST_CHECK
            if (!isConst) {worked=false; return "";}
        }        
        
        std::ostringstream code;
        
        code << "///////////////// START " << fullName << " (optimized) ////////////////\n";
        code << "\n";
        
        code << "__constant float " << functionName << "_b400[" << layeredFunctionIceCube.size() << "] = {\n";
        BOOST_FOREACH(const I3CLSimFunctionScatLenIceCubeConstPtr &function, layeredFunctionIceCube)
        {
            code << "    " << ToFloatString(function->GetB400()) << ",\n";
        }
        code << "};\n";
        code << "\n";
        
        code << "inline float " << functionName << "(unsigned int layer, float wlen);\n\n";
        code << "inline float " << functionName << "(unsigned int layer, float wlen)\n";
        code << "{\n";
        
        const std::string refWlenAsString = ToFloatString(1./(400.*I3Units::nanometer));

        code << "    const float alpha = " << ToFloatString(layeredFunctionIceCube[0]->GetAlpha()) << ";\n";
        code << "    \n";
        code << "#ifdef USE_NATIVE_MATH\n";
        code << "    return " << ToFloatString(I3Units::meter) << "*native_recip( " << functionName << "_b400[layer] * native_powr(wlen*" + refWlenAsString + ", -alpha) );\n";
        code << "#else\n";
        code << "    return " << ToFloatString(I3Units::meter) << "/( " << functionName << "_b400[layer] * powr(wlen*" + refWlenAsString + ", -alpha) );\n";
        code << "#endif\n";
        
        code << "}\n";
        code << "\n";
        
        code << "///////////////// END " << fullName << " (optimized) ////////////////\n";
        code << "\n";        

        worked=true;
        return code.str();
    }
    
};


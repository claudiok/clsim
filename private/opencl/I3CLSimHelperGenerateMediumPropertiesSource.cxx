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
 * @file I3CLSimHelperGenerateMediumPropertiesSource.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"

#include <string>
#include <sstream>
#include <stdexcept>

#include "dataclasses/I3Constants.h"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>


#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.h"

#include "clsim/I3CLSimHelperToFloatString.h"


namespace I3CLSimHelper
{
    
    template<class T>
    bool is_in_vector(const T &compVal, const std::vector<T> &vect, std::size_t &pos)
    {
        for (std::size_t i=0;i<vect.size();++i)
        {
            if (vect[i]==compVal)
            {
                pos=i;
                return true;
            }
            
        }
        return false;
    }
    
    template<class T>
    void OptimizeLayeredValue(const std::vector<boost::shared_ptr<const T> > &inputValues, // input
                              std::vector<boost::shared_ptr<const T> > &outputFunctions,   // output
                              std::vector<std::size_t> &outputFunctionForLayer)     // output
    {
        outputFunctions.clear();
        outputFunctionForLayer.clear();
        
        BOOST_FOREACH(const boost::shared_ptr<const T> &inputValue, inputValues)
        {
            std::size_t pos;
            if (is_in_vector(inputValue, outputFunctions, pos))
            {
                outputFunctionForLayer.push_back(pos);
            }
            else
            {
                outputFunctions.push_back(inputValue);
                outputFunctionForLayer.push_back(outputFunctions.size()-1);
            }
        }
        
    }
    
    std::string GenerateLayeredWlenDependentFunctions_WriteMainFuncCode(const std::vector<std::size_t> &outputFunctionForLayer,
                                                                        const std::vector<I3CLSimFunctionConstPtr> &functionsToGenerate,
                                                                        const std::string &fullName,
                                                                        const std::string &functionName)
    {
        std::ostringstream code;

        if (functionsToGenerate.size()==1)
        {
            code << "#define FUNCTION_" << functionName << "_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
        }
        
        code << "inline float " << functionName << "(unsigned int layer, float wavelength);\n\n";
        code << "inline float " << functionName << "(unsigned int layer, float wavelength)\n";
        code << "{\n";
        if (functionsToGenerate.size()==1) {
            // only a single function
            code << "    // " << fullName << " does not have a layer structure\n";
            code << "    return " << functionName << "_func0(wavelength);\n";
        } else {
            code << "    switch(layer)\n";
            code << "    {\n";
            for (std::size_t i=0;i<outputFunctionForLayer.size();++i)
            {
                const std::string iAsString = boost::lexical_cast<std::string>(i);
                code << "        case " << iAsString << ": return " << functionName << "_func" << outputFunctionForLayer.at(i) << "(wavelength);\n";
            }
            code << "        default: return 0.;\n";
            code << "    }\n";
            code << "\n";
        }
        code << "}\n";
        
        return code.str();
    }
                                                                        
    
    std::string GenerateLayeredWlenDependentFunctions(const std::vector<I3CLSimFunctionConstPtr> &layeredFunction,
                                                      const std::string &fullName,
                                                      const std::string &functionName,
                                                      std::string derivativeFunctionName="")
    {
        // first, check if one of the optimizers work
        if (derivativeFunctionName=="")
        {
            bool worked;
            std::string ret;
            
            ret = GenerateOptimizedCodeFor_I3CLSimFunctionAbsLenIceCube(layeredFunction,
                                                                                  fullName,
                                                                                  functionName,
                                                                                  worked);
            if (worked) return ret;

            ret = GenerateOptimizedCodeFor_I3CLSimFunctionScatLenIceCube(layeredFunction,
                                                                                   fullName,
                                                                                   functionName,
                                                                                   worked);
            if (worked) return ret;
            
        }        

        
        std::ostringstream code;

        code << "///////////////// START " << fullName << " ////////////////\n";
        code << "\n";
        
        std::vector<I3CLSimFunctionConstPtr> functionsToGenerate;
        std::vector<std::size_t> outputFunctionForLayer;
        
        OptimizeLayeredValue(layeredFunction,
                             functionsToGenerate,
                             outputFunctionForLayer);
        
        for (uint32_t i=0;i<functionsToGenerate.size();++i)
        {
            I3CLSimFunctionConstPtr currentValues = functionsToGenerate[i];
            if (!currentValues) log_fatal("%s function %u is (null)", fullName.c_str(), static_cast<unsigned int>(i));
            
            code << currentValues->GetOpenCLFunction(functionName+"_func"+boost::lexical_cast<std::string>(i));
            code << "\n";
            
            if (derivativeFunctionName != "")
            {
                if (!currentValues->HasDerivative()) log_fatal("%s function %u does not have a derivative", fullName.c_str(), static_cast<unsigned int>(i));

                code << currentValues->GetOpenCLFunctionDerivative(derivativeFunctionName+"_func"+boost::lexical_cast<std::string>(i));
                code << "\n";
            }
        }
        code << "\n";

        code << GenerateLayeredWlenDependentFunctions_WriteMainFuncCode(outputFunctionForLayer,
                                                                        functionsToGenerate,
                                                                        fullName,
                                                                        functionName);
        code << "\n";

        if (derivativeFunctionName != "")
        {
            code << GenerateLayeredWlenDependentFunctions_WriteMainFuncCode(outputFunctionForLayer,
                                                                            functionsToGenerate,
                                                                            fullName,
                                                                            derivativeFunctionName);
            code << "\n";
        }
        
        code << "///////////////// END " << fullName << " ////////////////\n";
        code << "\n";        

        
        return code.str();
    }
    
    
    std::string GenerateMediumPropertiesSource(const I3CLSimMediumProperties &mediumProperties)
    {
        std::ostringstream code;
        
        code << "// ice/water properties, auto-generated by\n";
        code << "// I3CLSimHelper::GenerateMediumPropertiesSource()\n";
        code << "\n";
        
        
        code << "#define MEDIUM_LAYERS " << mediumProperties.GetLayersNum() << "\n";
        code << "\n";
        code << "#define MEDIUM_MIN_WLEN " << ToFloatString(mediumProperties.GetMinWavelength()) << "\n";
        code << "#define MEDIUM_MAX_WLEN " << ToFloatString(mediumProperties.GetMaxWavelength()) << "\n";
        code << "\n";
        code << "#define MEDIUM_MIN_RECIP_WLEN " << ToFloatString(1./mediumProperties.GetMaxWavelength()) << "\n";
        code << "#define MEDIUM_MAX_RECIP_WLEN " << ToFloatString(1./mediumProperties.GetMinWavelength()) << "\n";
        code << "\n";
        code << "// medium layer structure:\n";
        code << "#define MEDIUM_LAYER_BOTTOM_POS " << ToFloatString(mediumProperties.GetLayersZStart()) << "\n";
        code << "#define MEDIUM_LAYER_THICKNESS  " << ToFloatString(mediumProperties.GetLayersHeight()) << "\n";
        code << "\n";

        // phase refractive index
        if (mediumProperties.GetPhaseRefractiveIndices()[0]->HasDerivative())
        {
            code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetPhaseRefractiveIndices(),
                                                          "phase refractive index",
                                                          "getPhaseRefIndex",
                                                          "getDispersion");
        } else {
            code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetPhaseRefractiveIndices(),
                                                          "phase refractive index",
                                                          "getPhaseRefIndex",
                                                          "");
        }
        
        if (mediumProperties.GetGroupRefractiveIndexOverride(0))
        {
            // seems there is an override parameterization for the group refractive index.
            // Use that instead of the one calculated from the dispersion.

            // sanity check
            for (uint32_t i=0;i<mediumProperties.GetLayersNum();++i)
            {
                if (!mediumProperties.GetGroupRefractiveIndexOverride(i))
                    log_fatal("Medium property error: group refractive index overrides are not set for all layers! (unset for layer %" PRIu32 ")", i);
            }
            
            // phase refractive index
            code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetGroupRefractiveIndicesOverride(),
                                                          "group refractive index",
                                                          "getGroupRefIndex");

            code << "#ifdef FUNCTION_getGroupRefIndex_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
            code << "#define FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
            code << "#endif" << std::endl;
            code << "// group velocity from group refractive index\n";
            code << "inline float getGroupVelocity(unsigned int layer, float wavelength);\n\n";
            code << "inline float getGroupVelocity(unsigned int layer, float wavelength)\n";
            code << "{\n";
            code << "    const float c_light = " << ToFloatString(I3Constants::c) << ";\n";
            code << "    const float n_group = getGroupRefIndex(layer, wavelength);\n";
            code << "    \n";
            code << "    return c_light / n_group;\n";
            code << "}\n";
        }
        else
        {    
            // group velocity from dispersion
            code << "#ifdef FUNCTION_getPhaseRefIndex_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
            code << "//implies: FUNCTION_getDispersion_DOES_NOT_DEPEND_ON_LAYER"<< std::endl;
            code << "#define FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
            code << "#define FUNCTION_getGroupRefIndex_DOES_NOT_DEPEND_ON_LAYER" << std::endl;
            code << "#endif" << std::endl;

            code << "// group velocity from dispersion\n";
            code << "inline float getGroupVelocity(unsigned int layer, float wavelength);\n\n";
            code << "inline float getGroupVelocity(unsigned int layer, float wavelength)\n";
            code << "{\n";
            code << "    const float c_light = " << ToFloatString(I3Constants::c) << ";\n";
            code << "    \n";
            code << "#ifdef USE_NATIVE_MATH\n";
            code << "    const float n_inv = native_recip(getPhaseRefIndex(layer, wavelength));\n";
            code << "#else\n";
            code << "    const float n_inv = 1.f/getPhaseRefIndex(layer, wavelength);\n";
            code << "#endif\n";
            code << "    \n";
            code << "    const float y = getDispersion(layer, wavelength);\n";
            code << "    \n";
            code << "    return c_light * (1.0f + y*wavelength*n_inv) * n_inv;\n";
            code << "}\n";
            code << "\n";        
            code << "inline float getGroupRefIndex(unsigned int layer, float wavelength);\n\n";
            code << "inline float getGroupRefIndex(unsigned int layer, float wavelength)\n";
            code << "{\n";
            code << "    const float c_light = " << ToFloatString(I3Constants::c) << ";\n";
            code << "    const float groupvel = getGroupVelocity(layer, wavelength);\n";
            code << "    \n";
            code << "    return c_light / groupvel;\n";
            code << "}\n";
            code << "\n";        
        }
    
        // scattering length
        code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetScatteringLengths(),
                                                      "scattering length",
                                                      "getScatteringLength");
        
        // absorption length
        code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetAbsorptionLengths(),
                                                      "absorption length",
                                                      "getAbsorptionLength");
        
        
        // scattering angle distribution
        {
            code << "///////////////// START scattering angle distribution ////////////////\n";
            code << "\n";

            I3CLSimRandomValueConstPtr scatAngles = mediumProperties.GetScatteringCosAngleDistribution();
            if (!scatAngles) log_fatal("scattering angle function is (null).");
            
            
            code << scatAngles->GetOpenCLFunction("makeScatteringCosAngle", // name
                                                  
                                                  // these are all defined as macros by the rng code:
                                                  "RNG_ARGS",               // function arguments for rng
                                                  "RNG_ARGS_TO_CALL",       // if we call anothor function, this is how we pass on the rng state
                                                  "RNG_CALL_UNIFORM_CO",    // the call to the rng for creating a uniform number [0;1[
                                                  "RNG_CALL_UNIFORM_OC"     // the call to the rng for creating a uniform number ]0;1]
                                                  );
            
            
            code << "///////////////// END scattering angle distribution ////////////////\n";
            code << "\n";        
        }
        
        // directional absorption length correction
        {
            code << "///////////////// START directional absorption length correction function ////////////////\n";
            code << "\n";
            I3CLSimScalarFieldConstPtr dirAbsLenCorr = mediumProperties.GetDirectionalAbsorptionLengthCorrection();
            if (!dirAbsLenCorr) log_fatal("directional absorption length correction function is (null).");
            code << dirAbsLenCorr->GetOpenCLFunction("getDirectionalAbsLenCorrFactor"); // name
            code << "///////////////// END directional absorption length correction function ////////////////\n";
            code << "\n";        
        }

        // pre-scattering direction transformation
        {
            code << "///////////////// START pre-scattering direction transformation ////////////////\n";
            code << "\n";
            I3CLSimVectorTransformConstPtr preScatterDirTransform = mediumProperties.GetPreScatterDirectionTransform();
            if (!preScatterDirTransform) log_fatal("pre-scattering direction transformation is (null).");
            code << preScatterDirTransform->GetOpenCLFunction("transformDirectionPreScatter"); // name
            code << "///////////////// END pre-scattering direction transformation ////////////////\n";
            code << "\n";        
        }

        // post-scattering direction transformation
        {
            code << "///////////////// START post-scattering direction transformation ////////////////\n";
            code << "\n";
            I3CLSimVectorTransformConstPtr postScatterDirTransform = mediumProperties.GetPostScatterDirectionTransform();
            if (!postScatterDirTransform) log_fatal("post-scattering direction transformation is (null).");
            code << postScatterDirTransform->GetOpenCLFunction("transformDirectionPostScatter"); // name
            code << "///////////////// END post-scattering direction transformation ////////////////\n";
            code << "\n";        
        }

        // ice tilt z-shift
        {
            code << "///////////////// START ice tilt z-shift ////////////////\n";
            code << "\n";
            I3CLSimScalarFieldConstPtr iceTiltZShift = mediumProperties.GetIceTiltZShift();
            if (!iceTiltZShift) log_fatal("ice tilt z-shift (null).");
            code << iceTiltZShift->GetOpenCLFunction("getTiltZShift"); // name
            code << "///////////////// END ice tilt z-shift ////////////////\n";
            code << "\n";        
        }
        
        return code.str();
    }

namespace {
    std::string makeGenerateWavelengthMasterFunction(std::size_t num,
                                                     const std::string &functionName,
                                                     const std::string &functionArgs,
                                                     const std::string &functionArgsToCall)
    {
        std::string ret = std::string("inline float ") + functionName + "(uint number, " + functionArgs + ");\n\n";
        ret = ret + std::string("inline float ") + functionName + "(uint number, " + functionArgs + ")\n";
        ret = ret + "{\n";

        if (num==0) {
            ret = ret + "    return 0.f;\n";
            ret = ret + "}\n";
            return ret;
        }

        if (num==1) {
            ret = ret + "    return " + functionName + "_0(" + functionArgsToCall + ");\n";
            ret = ret + "}\n";
            return ret;
        }

        // num>=2:
        
        ret = ret + "    if (number==0) {\n";
        ret = ret + "        return " + functionName + "_0(" + functionArgsToCall + ");\n";
        
        for (std::size_t i=1;i<num;++i)
        {
            ret = ret + "    } else if (number==" + boost::lexical_cast<std::string>(i) + ") {\n";
            ret = ret + "        return " + functionName + "_" + boost::lexical_cast<std::string>(i) + "(" + functionArgsToCall + ");\n";
        }

        ret = ret + "    } else {\n";
        ret = ret + "        return 0.f;\n";
        ret = ret + "    }\n";
        
        
        ret = ret + "}\n\n";

        return ret;
    }
}

std::string
GenerateWavelengthGeneratorSource(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators)
{
    std::string wlenGeneratorSource;
    for (std::size_t i=0; i<wlenGenerators.size(); ++i)
    {
        const std::string generatorName = "generateWavelength_" + boost::lexical_cast<std::string>(i);
        const std::string thisGeneratorSource = 
        wlenGenerators[i]->GetOpenCLFunction(generatorName, // name
                                              // these are all defined as macros by the rng code:
                                              "RNG_ARGS",               // function arguments for rng
                                              "RNG_ARGS_TO_CALL",       // if we call anothor function, this is how we pass on the rng state
                                              "RNG_CALL_UNIFORM_CO",    // the call to the rng for creating a uniform number [0;1[
                                              "RNG_CALL_UNIFORM_OC"     // the call to the rng for creating a uniform number ]0;1]
                                              );
        
        wlenGeneratorSource = wlenGeneratorSource + thisGeneratorSource + "\n";
    }
    wlenGeneratorSource = wlenGeneratorSource+ makeGenerateWavelengthMasterFunction(wlenGenerators.size(),
                                                                                      "generateWavelength",
                                                                                      "RNG_ARGS",               // function arguments for rng
                                                                                      "RNG_ARGS_TO_CALL"        // if we call anothor function, this is how we pass on the rng state
                                                                                      );
    wlenGeneratorSource = wlenGeneratorSource + "\n";
    
    return wlenGeneratorSource;
}


};

#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"

#include <string>
#include <sstream>
#include <stdexcept>

#include "dataclasses/I3Constants.h"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>


#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.h"

#include "clsim/to_float_string.h"


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
    void OptimizeLayeredValue(const std::vector<shared_ptr<const T> > &inputValues, // input
                              std::vector<shared_ptr<const T> > &outputFunctions,   // output
                              std::vector<std::size_t> &outputFunctionForLayer)     // output
    {
        outputFunctions.clear();
        outputFunctionForLayer.clear();
        
        BOOST_FOREACH(const shared_ptr<const T> &inputValue, inputValues)
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
                                                                        const std::vector<I3CLSimWlenDependentValueConstPtr> &functionsToGenerate,
                                                                        const std::string &fullName,
                                                                        const std::string &functionName)
    {
        std::ostringstream code;

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
                                                                        
    
    std::string GenerateLayeredWlenDependentFunctions(const std::vector<I3CLSimWlenDependentValueConstPtr> &layeredFunction,
                                                      const std::string &fullName,
                                                      const std::string &functionName,
                                                      std::string derivativeFunctionName="")
    {
        // first, check if one of the optimizers work (TODO: take them from a registry of some kind)
        if (derivativeFunctionName=="")
        {
            bool worked;
            std::string ret;
            
            ret = GenerateOptimizedCodeFor_I3CLSimWlenDependentValueAbsLenIceCube(layeredFunction,
                                                                                  fullName,
                                                                                  functionName,
                                                                                  worked);
            if (worked) return ret;

            ret = GenerateOptimizedCodeFor_I3CLSimWlenDependentValueScatLenIceCube(layeredFunction,
                                                                                   fullName,
                                                                                   functionName,
                                                                                   worked);
            if (worked) return ret;
            
        }        

        
        std::ostringstream code;

        code << "///////////////// START " << fullName << " ////////////////\n";
        code << "\n";
        
        std::vector<I3CLSimWlenDependentValueConstPtr> functionsToGenerate;
        std::vector<std::size_t> outputFunctionForLayer;
        
        OptimizeLayeredValue(layeredFunction,
                             functionsToGenerate,
                             outputFunctionForLayer);
        
        for (uint32_t i=0;i<functionsToGenerate.size();++i)
        {
            I3CLSimWlenDependentValueConstPtr currentValues = functionsToGenerate[i];
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
        code << "#define MEDIUM_MIN_WLEN " << to_float_string(mediumProperties.GetMinWavelength()) << "\n";
        code << "#define MEDIUM_MAX_WLEN " << to_float_string(mediumProperties.GetMaxWavelength()) << "\n";
        code << "\n";
        code << "#define MEDIUM_MIN_RECIP_WLEN " << to_float_string(1./mediumProperties.GetMaxWavelength()) << "\n";
        code << "#define MEDIUM_MAX_RECIP_WLEN " << to_float_string(1./mediumProperties.GetMinWavelength()) << "\n";
        code << "\n";
        code << "// medium layer structure:\n";
        code << "#define MEDIUM_LAYER_BOTTOM_POS " << to_float_string(mediumProperties.GetLayersZStart()) << "\n";
        code << "#define MEDIUM_LAYER_THICKNESS  " << to_float_string(mediumProperties.GetLayersHeight()) << "\n";
        code << "\n";

        // phase refractive index
        code << GenerateLayeredWlenDependentFunctions(mediumProperties.GetPhaseRefractiveIndices(),
                                                      "phase refractive index",
                                                      "getPhaseRefIndex",
                                                      "getDispersion");

        // group velocity
        code << "// group velocity\n";
        code << "inline float getGroupVelocity(unsigned int layer, float wavelength)\n";
        code << "{\n";
        code << "    const float c_light = " << to_float_string(I3Constants::c) << ";\n";
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
        
        
        return code.str();
    }
    


};

#ifndef I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED
#define I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED

#include <string>
#include <vector>

#include "clsim/function/I3CLSimFunction.h"

namespace I3CLSimHelper
{

    // optimizers (special converters for functions
    // where we can generate more optimized code
    // for layered values)
    std::string GenerateOptimizedCodeFor_I3CLSimFunctionAbsLenIceCube(const std::vector<I3CLSimFunctionConstPtr> &layeredFunction,
                                                                                const std::string &fullName,
                                                                                const std::string &functionName,
                                                                                bool &worked);

    std::string GenerateOptimizedCodeFor_I3CLSimFunctionScatLenIceCube(const std::vector<I3CLSimFunctionConstPtr> &layeredFunction,
                                                                                 const std::string &fullName,
                                                                                 const std::string &functionName,
                                                                                 bool &worked);
    

};

#endif //I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED

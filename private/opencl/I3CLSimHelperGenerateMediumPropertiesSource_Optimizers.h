#ifndef I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED
#define I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED

#include <string>
#include <vector>

#include "clsim/I3CLSimWlenDependentValue.h"

namespace I3CLSimHelper
{

    // optimizers (special converters for functions
    // where we can generate more optimized code
    // for layered values)
    std::string GenerateOptimizedCodeFor_I3CLSimWlenDependentValueAbsLenIceCube(const std::vector<I3CLSimWlenDependentValueConstPtr> &layeredFunction,
                                                                                const std::string &fullName,
                                                                                const std::string &functionName,
                                                                                bool &worked);

    std::string GenerateOptimizedCodeFor_I3CLSimWlenDependentValueScatLenIceCube(const std::vector<I3CLSimWlenDependentValueConstPtr> &layeredFunction,
                                                                                 const std::string &fullName,
                                                                                 const std::string &functionName,
                                                                                 bool &worked);
    

};

#endif //I3CLSIMHELPERGENERATEMEDIUMPROPERTIESSOURCE_OPTIMIZERS_H_INCLUDED


#ifndef CLSIM_HELPER_MATH_H_INCLUDED
#define CLSIM_HELPER_MATH_H_INCLUDED

#include <string>

class I3CLSimOpenCLDevice;

namespace I3CLSimHelper {

std::string GetMathPreamble(const I3CLSimOpenCLDevice &device, bool double_precision=false);

}

#endif // CLSIM_HELPER_MATH_H_INCLUDED

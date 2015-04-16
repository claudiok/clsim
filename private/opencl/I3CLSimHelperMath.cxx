
#include "opencl/I3CLSimHelperMath.h"
#include "clsim/I3CLSimOpenCLDevice.h"

std::string
I3CLSimHelper::GetMathPreamble(const I3CLSimOpenCLDevice &device,
    bool doublePrecision)
{
	std::string preamble;

	// necessary OpenCL extensions
	preamble += \
	    "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"\
	    "#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable\n";
	
	// are we using double precision?
	if (doublePrecision) {
		preamble += "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"\
		            "typedef double floating_t;\n"                   \
		            "typedef double2 floating2_t;\n"                 \
		            "typedef double4 floating4_t;\n"                 \
		            "#define convert_floating_t convert_double\n"    \
		            "#define DOUBLE_PRECISION\n"                     \
		            "#define ZERO 0.\n"                              \
		            "#define ONE 1.\n";
		
		if (device.GetPlatformName()=="Apple") {
			preamble += "#define USE_FABS_WORKAROUND\n";
			log_info("enabled fabs() workaround for OpenCL double-precision on Apple");
		}
	} else {
		preamble += "typedef float floating_t;\n"                    \
		            "typedef float2 floating2_t;\n"                  \
		            "typedef float4 floating4_t;\n"                  \
		            "#define convert_floating_t convert_float\n"     \
		            "#define ZERO 0.f\n"                             \
		            "#define ONE 1.f\n";
	}
	preamble += "\n";
	
	return preamble;
}

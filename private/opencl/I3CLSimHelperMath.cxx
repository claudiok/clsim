
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
	preamble += "#ifdef DOUBLE_PRECISION\n"
	            "// can't have native_math with double precision\n"
	            "#ifdef USE_NATIVE_MATH\n"
	            "#undef USE_NATIVE_MATH\n"
	            "#endif\n"
	            "#endif\n"
	            "\n"
	            "#ifdef USE_NATIVE_MATH\n"
	            "inline floating_t my_divide(floating_t a, floating_t b) {return native_divide(a,b);}\n"
	            "inline floating_t my_recip(floating_t a) {return native_recip(a);}\n"
	            "inline floating_t my_powr(floating_t a, floating_t b) {return native_powr(a,b);}\n"
	            "inline floating_t my_sqrt(floating_t a) {return native_sqrt(a);}\n"
	            "inline floating_t my_rsqrt(floating_t a) {return native_rsqrt(a);}\n"
	            "inline floating_t my_cos(floating_t a) {return native_cos(a);}\n"
	            "inline floating_t my_sin(floating_t a) {return native_sin(a);}\n"
	            "inline floating_t my_log(floating_t a) {return native_log(a);}\n"
	            "inline floating_t my_exp(floating_t a) {return native_exp(a);}\n"
	            "#else\n"
	            "inline floating_t my_divide(floating_t a, floating_t b) {return a/b;}\n"
	            "inline floating_t my_recip(floating_t a) {return 1.f/a;}\n"
	            "inline floating_t my_powr(floating_t a, floating_t b) {return powr(a,b);}\n"
	            "inline floating_t my_sqrt(floating_t a) {return sqrt(a);}\n"
	            "inline floating_t my_rsqrt(floating_t a) {return rsqrt(a);}\n"
	            "inline floating_t my_cos(floating_t a) {return cos(a);}\n"
	            "inline floating_t my_sin(floating_t a) {return sin(a);}\n"
	            "inline floating_t my_log(floating_t a) {return log(a);}\n"
	            "inline floating_t my_exp(floating_t a) {return exp(a);}\n"
	            "#endif\n"
	            "\n"
	            "#ifdef USE_FABS_WORKAROUND\n"
	            "inline floating_t my_fabs(floating_t a) {return (a<ZERO)?(-a):(a);}\n"
	            "#else\n"
	            "inline floating_t my_fabs(floating_t a) {return fabs(a);}\n"
	            "#endif\n"
	            "inline floating_t sqr(floating_t a) {return a*a;}\n"
	;
	
	return preamble;
}

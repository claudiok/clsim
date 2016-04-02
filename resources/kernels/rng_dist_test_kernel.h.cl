__kernel void testKernel(__global float* randomNumbers,
                         __global ulong* MWC_RNG_x,
                         __global uint* MWC_RNG_a);

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable

// // disable dbg_printf for GPU
// #define dbg_printf(format, ...)
// 
// // enable printf for CPU
// //#pragma OPENCL EXTENSION cl_amd_printf : enable
// //#define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)


inline float my_divide(float a, float b);
inline float my_recip(float a);
inline float my_powr(float a, float b);
inline float my_sqrt(float a);
inline float my_cos(float a);
inline float my_sin(float a);
inline float my_log(float a);
inline float my_exp(float a);
inline float sqr(float a);

#ifdef USE_NATIVE_MATH
inline float my_divide(float a, float b) {return native_divide(a,b);}
inline float my_recip(float a) {return native_recip(a);}
inline float my_powr(float a, float b) {return native_powr(a,b);}
inline float my_sqrt(float a) {return native_sqrt(a);}
inline float my_cos(float a) {return native_cos(a);}
inline float my_sin(float a) {return native_sin(a);}
inline float my_log(float a) {return native_log(a);}
inline float my_exp(float a) {return native_exp(a);}
#else
inline float my_divide(float a, float b) {return a/b;}
inline float my_recip(float a) {return 1.f/a;}
inline float my_powr(float a, float b) {return powr(a,b);}
inline float my_sqrt(float a) {return sqrt(a);}
inline float my_cos(float a) {return cos(a);}
inline float my_sin(float a) {return sin(a);}
inline float my_log(float a) {return log(a);}
inline float my_exp(float a) {return exp(a);}
#endif

inline float sqr(float a) {return a*a;}

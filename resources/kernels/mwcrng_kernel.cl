// Multiply-with-carry random number generator for OpenCL using an
// implementation along the lines of the CUDAMCML code described here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf

// prototypes to make some compilers happy
inline float rand_MWC_co(ulong *x,uint *a);
inline float rand_MWC_oc(ulong *x,uint *a);

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1) 
//////////////////////////////////////////////////////////////////////////////
inline float rand_MWC_co(ulong *x,uint *a)
{
  *x=(*x&0xfffffffful)*(*a)+(*x>>32);
#ifdef USE_NATIVE_MATH
  return native_divide(convert_float_rtz((uint)(*x&0xfffffffful)),(float)0x100000000); // OpenCL - native divide
#else
  return (convert_float_rtz((uint)(*x&0xfffffffful))/(float)0x100000000); // OpenCL
#endif
} 

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 (0,1]
//////////////////////////////////////////////////////////////////////////////
inline float rand_MWC_oc(ulong *x,uint *a)
{
  return 1.0f-rand_MWC_co(x,a);
} 

// typedefs for later use
#define RNG_ARGS ulong *rnd_x,uint *rnd_a
#define RNG_ARGS_TO_CALL rnd_x,rnd_a
#define RNG_CALL_UNIFORM_CO rand_MWC_co(rnd_x,rnd_a)
#define RNG_CALL_UNIFORM_OC rand_MWC_oc(rnd_x,rnd_a)

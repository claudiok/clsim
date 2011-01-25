// This code is part of IceTray. It has been taken from
// GPUMCML. Here is its original copyright notice:

/*	 
 *   This file is part of GPUMCML.
 * 
 *   GPUMCML is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   GPUMCML is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with GPUMCML.  If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1) 
//////////////////////////////////////////////////////////////////////////////
inline float rand_MWC_co(ulong *x,uint *a)
{
  //*x=(*x&0xffffffffull)*(*a)+(*x>>32); // CUDA
  *x=(*x&0xfffffffful)*(*a)+(*x>>32);  // OpenCL
  //return __fdividef(__uint2float_rz((UINT32)(*x)),(FLOAT)0x100000000);  // CUDA
#ifdef USE_NATIVE_MATH
  return native_divide(convert_float_rtz((uint)(*x&0xfffffffful)),(float)0x100000000); // OpenCL - native divide
#else
  return (convert_float_rtz((uint)(*x&0xfffffffful))/(float)0x100000000); // OpenCL
#endif

  // The typecast will truncate the x so that it is 0<=x<(2^32-1),
  // __uint2FLOAT_rz ensures a round towards zero since 32-bit FLOATing point 
  // cannot represent all integers that large. 
  // Dividing by 2^32 will hence yield [0,1)
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



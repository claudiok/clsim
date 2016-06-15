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
 * @file propagation_kernel.h.cl
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

/////////////////// preprocessor defines

// debug mode: store all photons right after they are generated on string0/OM0
//#define DEBUG_STORE_GENERATED_PHOTONS

// this will be (optinally) defined by the main code:
//#define STOP_PHOTONS_ON_DETECTION

//#define PRINTF_ENABLED

// // disable dbg_printf for GPU
// //#define dbg_printf(format, ...)

// #ifdef PRINTF_ENABLED
// enable printf for CPU
// #pragma OPENCL EXTENSION cl_amd_printf : enable
// #define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)
// #else
// #define dbg_printf(format, ...)
// #endif

// ZERO and ONE will be defined as either 0.f/1.f or 0./1. depending on DOUBLE_PRECISION

/////////////////// struct definitions

struct __attribute__ ((packed)) I3CLSimStep 
{
    float4 posAndTime;   // x,y,z,time                      // 4x 32bit float
    float4 dirAndLengthAndBeta; // theta,phi,length,beta    // 4x 32bit float
    uint numPhotons;                                        //    32bit unsigned
    float weight;                                           //    32bit float
    uint identifier;                                        //    32bit unsigned
    uchar sourceType;                                       //     8bit unsigned
    uchar dummy1;                                           //     8bit unsigned
    ushort dummy2;                                          //    16bit unsigned
                                                            // total: 12x 32bit float = 48 bytes
};

struct __attribute__ ((packed)) I3CLSimPhoton 
{
    float4 posAndTime;   // x,y,z,time                      // 4x 32bit float
    float2 dir; // theta,phi                                // 2x 32bit float
    float wavelength; // photon wavelength                  //    32bit float
    float cherenkovDist; // Cherenkov distance travelled    //    32bit float
    uint numScatters; // number of scatters                 //    32bit unsigned
    float weight;                                           //    32bit float
    uint identifier;                                        //    32bit unsigned
    short stringID;                                         //    16bit signed
    ushort omID;                                            //    16bit unsigned
    float4 startPosAndTime;                                 // 4x 32bit float
    float2 startDir;                                        // 2x 32bit float
    float groupVelocity;                                    //    32bit float
    float distInAbsLens;                                    //    32bit float
                                                            // total: 20x 32bit float = 80 bytes
};

struct __attribute__ ((packed)) I3CLSimTableEntry
{
    uint index;
    float weight;
};

struct __attribute__ ((packed)) I3CLSimReferenceParticle
{
    float4 posAndTime;   // x,y,z,time
    float4 dir;          // dx,dy,dz,0
    float4 perpDir;
};

///////////////// forward declarations

inline int findLayerForGivenZPos(floating_t posZ);

inline floating_t mediumLayerBoundary(int layer);

void scatterDirectionByAngle(floating_t cosa,
    floating_t sina,
    floating4_t *direction,
    floating_t randomNumber);

inline void createPhotonFromTrack(struct I3CLSimStep *step,
    const floating4_t stepDir,
    RNG_ARGS,
    floating4_t *photonPosAndTime,
    floating4_t *photonDirAndWlen);

#ifdef DOUBLE_PRECISION
inline float2 sphDirFromCar(double4 carDir);
#else
inline float2 sphDirFromCar(float4 carDir);
#endif

inline void saveHit(
    const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
    const floating_t thisStepLength,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    floating_t distanceTraveledInAbsorptionLengths,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    unsigned short hitOnString,
    unsigned short hitOnDom,
    __global uint* hitIndex,
    uint maxHitIndex,
    __global struct I3CLSimPhoton *outputPhotons
#ifdef SAVE_PHOTON_HISTORY
  , __global float4 *photonHistory,
    float4 *currentPhotonHistory
#endif
    );

///////////////////////// some constants

#ifdef DOUBLE_PRECISION
__constant double speedOfLight = 0.299792458; // [m/ns]
__constant double recip_speedOfLight = 3.33564095; // [ns/m]
__constant double PI = 3.14159265359;
#else
__constant float speedOfLight = 0.299792458f; // [m/ns]
__constant float recip_speedOfLight = 3.33564095f; // [ns/m]
__constant float PI = 3.14159265359f;
#endif

///////////////////////////


inline floating_t my_divide(floating_t a, floating_t b);
inline floating_t my_recip(floating_t a);
inline floating_t my_powr(floating_t a, floating_t b);
inline floating_t my_sqrt(floating_t a);
inline floating_t my_rsqrt(floating_t a);
inline floating_t my_cos(floating_t a);
inline floating_t my_sin(floating_t a);
inline floating_t my_log(floating_t a);
inline floating_t my_exp(floating_t a);
inline floating_t my_fabs(floating_t a);
inline floating_t sqr(floating_t a);

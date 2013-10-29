// This code is based on the CUDA code described in
// the documentation of the "CUDAMCML" package which can be found here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf

#ifndef MWCRNG_INIT_H
#define MWCRNG_INIT_H

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include <cstdio>
#include <cstring>
#include <stdint.h>

#include <string>

#include "phys-services/I3RandomService.h"

// Initialize random number generator 
inline int init_MWC_RNG(uint64_t *x, uint32_t *a, 
                 const uint32_t n_rng,
                 std::string safeprimes_file,
                 I3RandomServicePtr randomService)
{
    FILE *fp;
    uint32_t fora,tmp1,tmp2;

    if (safeprimes_file == "")
    {
        // Try to find it in the local directory
        safeprimes_file = "safeprimes_base32.txt";
    }

    fp = fopen(safeprimes_file.c_str(), "r");

    if(fp == NULL)
    {
        log_error("Could not find the safeprimes file (%s)! Terminating!", safeprimes_file.c_str());
        return 1;
    }

    // Here we set up a loop, using the first multiplier in the file to generate x's and c's
    // There are some restictions to these two numbers:
    // 0<=c<a and 0<=x<b, where a is the multiplier and b is the base (2^32)
    // also [x,c]=[0,0] and [b-1,a-1] are not allowed.

    for (uint32_t i=0;i < n_rng;i++)
    {
        int ret = fscanf(fp,"%u %u %u",&fora,&tmp1,&tmp2);
        if (ret==EOF) 
        {
            log_error("Could not parse data! Terminating!");
            return 1;
        }

        a[i]=fora; // primes from file go to a[]

        // generate x[] from the supplied rng
        x[i]=0;
        while( (x[i]==0) | (((uint32_t)(x[i]>>32))>=(fora-1)) | (((uint32_t)x[i])>=0xfffffffful))
        {
            // generate a random numbers for x and x (both are stored in the "x" array
            x[i] = static_cast<uint32_t>(randomService->Integer(0xffffffff));
            x[i]=x[i]<<32;
            x[i] += static_cast<uint32_t>(randomService->Integer(0xffffffff));
        }
    }
    fclose(fp);

    return 0;
}

#endif

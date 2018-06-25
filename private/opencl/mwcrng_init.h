// Initialization of the multiply-with-carry random number generator for OpenCL using an
// implementation along the lines of the CUDAMCML code described here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf
//
// This code can generate a "safeprimes_base32.txt" file compatible with their implementation,
// but a binary file format with much smaller file sizes is also supported.

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
#include <boost/filesystem.hpp>
#include <icetray/open.h>

#include "phys-services/I3RandomService.h"

// Initialize random number generator 
inline int init_MWC_RNG(uint64_t *x, uint32_t *a, 
                 const uint32_t n_rng,
                 I3RandomServicePtr randomService,
                 std::string safeprimes_file="")
{
    int64_t multiplier;

    if (safeprimes_file == "")
    {
        bool success = 0;
        namespace fs = boost::filesystem;
        safeprimes_file = "safeprimes_base32.txt"; 
        if (getenv("I3_BUILD")) {
            const fs::path I3_BUILD(getenv("I3_BUILD"));
            if (fs::exists(I3_BUILD/"clsim/resources/safeprimes_base32.gz")) {
                safeprimes_file = (I3_BUILD/"/clsim/resources/safeprimes_base32.gz").string();
                success = 1;
            }
            else if (fs::exists(I3_BUILD/"clsim/resources/safeprimes_base32.txt")) {
              safeprimes_file = (I3_BUILD/"/clsim/resources/safeprimes_base32.txt").string();
              success = 1;
            }
          }
        if (!success && getenv("I3_DATA")) {
            const fs::path I3_DATA(getenv("I3_DATA"));
            if (fs::exists(I3_DATA/"safeprimes_base32.gz")) {
                safeprimes_file = (I3_DATA/"safeprimes_base32.gz").string();
            }
            else if (fs::exists(I3_DATA/"safeprimes_base32.txt")) {
              safeprimes_file = (I3_DATA/"safeprimes_base32.txt").string();
            }
          }
    }
    
    boost::iostreams::filtering_istream ifs;
    I3::dataio::open(ifs, safeprimes_file);
    if (!ifs.good()) {
        log_error("Could not find the safeprimes file (%s)! Terminating!", safeprimes_file.c_str());
        return 1;
    }
    
    bool plaintext = false;
    {
        // Detect newer binary file format
        char tag[18];
        ifs.read(tag, 17);
        tag[17] = '\0';
        if (strcmp(tag, "safeprimes_base32") != 0) {
            plaintext = true;
            I3::dataio::open(ifs, safeprimes_file);
        }
    }

    for (uint32_t i=0;i < n_rng;i++) {
        if (ifs.eof())
            log_error("File ended before %u primes could be read!", i+1);
        if (plaintext) {
            ifs >> multiplier;
            if (ifs.fail()) {
                log_error("Couldn't parse prime at line %u", i+1);
                return 1;
            }
            char c;
            while (!ifs.eof() && ifs.peek() != '\n')
                ifs.read(&c, 1);
        } else {
            ifs.read(reinterpret_cast<char*>(&multiplier), sizeof(multiplier));
            if (ifs.fail()) {
                log_error("Couldn't read prime #%u", i+1);
                return 1;
            }
        }

        if ((multiplier < std::numeric_limits<uint32_t>::min()) || (multiplier > std::numeric_limits<uint32_t>::max())) {
            log_error("Prime #%u (%" PRIi64 ") is out of range!", i+1, multiplier);
            return 1;
        }
        
        a[i]=multiplier; // primes from file go to a[]

        // generate x[] from the supplied rng
        x[i]=0;
        while( (x[i]==0) | (((uint32_t)(x[i]>>32))>=(a[i]-1)) | (((uint32_t)x[i])>=0xfffffffful))
        {
            x[i] = static_cast<uint32_t>(randomService->Integer(0xffffffff));
            x[i]=x[i]<<32;
            x[i] += static_cast<uint32_t>(randomService->Integer(0xffffffff));
        }
    }
    
    return 0;
}

#endif

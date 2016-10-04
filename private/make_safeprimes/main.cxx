// Initialization of the multiply-with-carry random number generator for OpenCL using an
// implementation along the lines of the CUDAMCML code described here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf
//
// the output file should be compatible to CUDAMCML's output

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <iostream>
#include <fstream>
#include <cstddef>

// this uses GMP for arbitrary precision math and Miller-Rabin probabilistic primality tests
#include "gmp.h"

// Probabilistic test if "n" is prime.
bool is_prime(mpz_t n)
{
    int result = mpz_probab_prime_p(n, 100);

    // 2: "definitely prime"
    // 1: "probably prime"
    // 0: "definitely composite"
    return (result==2) || (result==1);
}

// Make a CUDAMCML-compatible safeprimes file (only the multipliers are actually
// needed).
int main()
{
    std::ofstream output_file("safeprimes_base32.txt", std::ios_base::out);

    const std::size_t num_entries = 131072;
    const std::size_t stringbuf_size = 1000;

    mpz_t n1,n2,a;
    char stringbuf[stringbuf_size];

    mpz_init(n1);
    mpz_init(n2);
    mpz_init(a);

    mpz_set_str(a,"4294967118",0);

    std::size_t i=0;

    std::cout << "making " << num_entries << " multipliers." << std::endl;

    for (;;)
    {
        mpz_mul_2exp(n2,a,32);
        mpz_sub_ui(n2,n2,1); // n2 = n2-1

        if (!is_prime(n2))
        {
            // n2 is not a prime number. skip it.
            mpz_sub_ui(a,a,1); // a = a-1
            continue;
        }

        mpz_sub_ui(n1,n2,1); // n1 = n2-1
        mpz_fdiv_q_2exp(n1,n1,1);

        if(!is_prime(n1))
        {
            // n1 is not a prime number. skip it.
            mpz_sub_ui(a,a,1); // a = a-1
            continue;
        }

        // we have a valid value for the multiplier a, save it (and the two primes)
        gmp_snprintf(stringbuf, stringbuf_size, "%Zd", a);
        output_file << stringbuf << " ";

        gmp_snprintf(stringbuf, stringbuf_size, "%Zd", n2);
        output_file << stringbuf << " ";

        gmp_snprintf(stringbuf, stringbuf_size, "%Zd", n1);
        output_file << stringbuf << std::endl;

        i++;
        if (i >= num_entries)
        {
            // done.
            break;
        }

        if (i % 1000 == 0)
        {
            const double percentage = 100.*static_cast<double>(i)/static_cast<double>(num_entries);
            std::cout << percentage << "% done. (" << i << "/" << num_entries << ")" << std::endl;
        }

        mpz_sub_ui(a,a,1lu); // a = a-1
    }

    output_file.close();

    std::cout << "finished." << std::endl;

    return 0;
}

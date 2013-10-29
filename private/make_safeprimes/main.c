// This code is based on the CUDA code described in
// the documentation of the "CUDAMCML" package which can be found here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gmp.h"

int isprime(mpz_t n)
{
    int i;
    int class;
    int alist[]= {2,3,5,7,11,13,17,19,23,29,31,37,41};
    for(i=0; i<13; i++)
    {
        class = mpz_probab_prime_p (n, alist[i]);
        if(class==0) return 0;
    }
    return 1;
}

int main()
{
    FILE* file;
    mpz_t n1,n2,a;
    unsigned long long i=0;
    mpz_init(n1);
    mpz_init(n2);
    mpz_init(a);

    mpz_set_str(a,"4294967118",0);

    file = fopen("safeprimes_base32.txt","w");

    while(i<15000000)
    {
        mpz_mul_2exp(n2,a,32);
        mpz_sub_ui(n2,n2,1lu);
        if(isprime(n2))

        {

            mpz_sub_ui(n1,n2,1lu);
            mpz_fdiv_q_2exp(n1,n1,1lu);
            if(isprime(n1))
            {
                //We have found our safeprime, calculate a and print to file
                mpz_out_str (file, 10, a);
                fprintf(file," ");
                mpz_out_str (file, 10, n2);
                fprintf(file," ");
                mpz_out_str (file, 10, n1);
                fprintf(file,"\n");
                printf("%llu\n",i++);
            }
        }
        mpz_sub_ui(a,a,1lu);
    }
    fclose(file);
    exit (0);
}

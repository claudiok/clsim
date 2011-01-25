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


#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gmp/gmp.h"

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
    
    //while(i<150000)
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



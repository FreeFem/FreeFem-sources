// SUMMARY  :   matrix manipulation
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : P. Jolivet 
// E-MAIL   : Pierre Jolivet <pierre.jolivet@ljll.math.upmc.fr>
//

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 */

#include <vector>
//#include <mpi.h>
#include <iostream>
#include <numeric>

struct step {
    public:
        step(int x, int y) : x(x), y(y) { }
        int operator()() { return x += y; }

    private:
        int x, y;
};

template<unsigned char M, unsigned char S>
static void CSR2COO(unsigned int n, int* compressedI, int* uncompressedI) {
    if(S == 'U') {
        for(int i = n - 1; i > -1; --i) {
            if(M == 'F')
                std::fill(uncompressedI + compressedI[i] - 1, uncompressedI + compressedI[i + 1] - 1, i + 1);
            else
                std::fill(uncompressedI + compressedI[i], uncompressedI + compressedI[i + 1], i + 1);
        }
    }
    else if(S == 'L') {
        for(int i = 1; i < n; ++i) {
            if(M == 'F')
                std::fill(uncompressedI + compressedI[i] - i - 1, uncompressedI + compressedI[i + 1] - i - 1, i + 1);
            else
                std::fill(uncompressedI + compressedI[i] - i, uncompressedI + compressedI[i + 1] - i - 1, i + 1);
        }
    }
};

template<bool WithDiagonal, unsigned char N,typename Scalar>
static unsigned int trimCSR(unsigned int n, int* trimmedI, int* untrimmedI, int* trimmedJ, int* untrimmedJ, Scalar* trimmedC, Scalar* untrimmedC) {
    unsigned int upper = 0;
    for(unsigned int i = 0; i < n - WithDiagonal; ++i) {
        trimmedI[i] = upper + (N == 'F');
        int* jIndex = lower_bound(untrimmedJ + untrimmedI[i], untrimmedJ + untrimmedI[i + 1], i + !WithDiagonal);
        unsigned int j = untrimmedI[i] + jIndex - (untrimmedJ + untrimmedI[i]);
        if(N == 'F') {
            for(unsigned int k = j; k < untrimmedI[i + 1]; ++k)
                trimmedJ[upper + k - j] = untrimmedJ[k] + 1;

        } else
            std::copy(untrimmedJ + j, untrimmedJ + untrimmedI[i + 1], trimmedJ + upper);
        std::copy(untrimmedC + j, untrimmedC + untrimmedI[i + 1], trimmedC + upper);
        upper += untrimmedI[i + 1] - j;
    }
    if(WithDiagonal) {
        trimmedI[n - 1] = upper + (N == 'F');
        trimmedI[n] = trimmedI[n - 1] + 1;
        trimmedJ[upper] = n - (N == 'C');
        trimmedC[upper] = untrimmedC[untrimmedI[n] - 1];
        return trimmedI[n];
    }
    else {
        trimmedI[n] = trimmedI[n - 1];
        return trimmedI[n + 1];
    }
};

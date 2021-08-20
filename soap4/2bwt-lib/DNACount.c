/*

   DNACount.c        DNA Count

   This module contains DNA occurrence counting functions. The DNA must be
   in word-packed format.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include "DNACount.h"
#include "HSP.h"
#include "MiscUtilities.h"


void GenerateDNAOccCountTable(unsigned int *dnaDecodeTable) {

    unsigned int i, j, c, t;

    for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
        dnaDecodeTable[i] = 0;
        c = i;
        for (j=0; j<DNA_OCC_CNT_TABLE_SIZE; j++) { // log(alphabet_size, DNA_OCC_CNT_TABLE_SIZE_IN_WORD)
            t = c & 0x00000003;
            dnaDecodeTable[i] += 1 << (t * DNA_OCC_CNT_TABLE_SIZE);
            c >>= DNA_OCC_CNT_TABLE_SIZE / ALPHABET_SIZE; // log(alphabet_size, DNA_OCC_CNT_TABLE_SIZE_IN_WORD) / alphabet_size
        } 
    }

}

void ForwardDNAAllOccCountNoLimit(const unsigned int*  dna, const unsigned long long index, unsigned long long* __restrict occCount, const unsigned int*  dnaDecodeTable) {

    static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
                                               0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
                                               0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
                                               0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

    unsigned long long iteration, wordToCount, charToCount;
    unsigned long long i, j, c;
    unsigned long long sum;

    occCount[0] = 0;
    occCount[1] = 0;
    occCount[2] = 0;
    occCount[3] = 0;

    iteration = index / 256;
    wordToCount = (index - iteration * 256) / 16;
    charToCount = index - iteration * 256 - wordToCount * 16;

    for (i=0; i<iteration; i++) {

        sum = 0;
        for (j=0; j<16; j++) {
            sum += dnaDecodeTable[*dna >> 16];
            sum += dnaDecodeTable[*dna & 0x0000FFFF];
            dna++;
        }
        if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            occCount[0] += sum & 0x000000FF;    sum >>= 8;
            occCount[1] += sum & 0x000000FF;    sum >>= 8;
            occCount[2] += sum & 0x000000FF;    sum >>= 8;
            occCount[3] += sum;
        } else {
            // only some or all of the 3 bits are on
            // in reality, only one of the four cases are possible
            if (sum == 0x00000100) {
                occCount[0] += 256;
            } else if (sum == 0x00010000) {
                occCount[1] += 256;
            } else if (sum == 0x01000000) {
                occCount[2] += 256;
            } else if (sum == 0x00000000) {
                occCount[3] += 256;
            } else {
                fprintf(stderr, "ForwardDNAAllOccCountNoLimit(): DNA occ sum exception!\n");
                exit(1);
            }
        }

    }

    sum = 0;
    for (j=0; j<wordToCount; j++) {
        sum += dnaDecodeTable[*dna >> 16];
        sum += dnaDecodeTable[*dna & 0x0000FFFF];
        dna++;
    }

    if (charToCount > 0) {
        c = *dna & truncateRightMask[charToCount];    // increase count of 'a' by 16 - c;
        sum += dnaDecodeTable[c >> 16];
        sum += dnaDecodeTable[c & 0xFFFF];
        sum += charToCount - 16;    // decrease count of 'a' by 16 - positionToProcess
    }

    occCount[0] += sum & 0x000000FF;    sum >>= 8;
    occCount[1] += sum & 0x000000FF;    sum >>= 8;
    occCount[2] += sum & 0x000000FF;    sum >>= 8;
    occCount[3] += sum;

}
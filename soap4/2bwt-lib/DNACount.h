/*

   DNACount.h        DNA Count

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

#ifndef __DNA_COUNT_H__
#define __DNA_COUNT_H__

#include <stdint.h>
#include "TypeNLimit.h"

// DNA
#define DNA_ALPHABET_SIZE            4
#define DNA_CHAR_PER_WORD            16
#define DNA_BIT_PER_CHAR            2

// DNA occurrence count table
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD    65536 // this is pow(ALPHABET_SIZE, DNA_OCC_CNT_TABLE_SIZE)
#define DNA_OCC_CNT_TABLE_SIZE   8
#define DNA_OCC_SUM_EXCEPTION(sum)            ((sum & 0xfefefeff) == 0)


void GenerateDNAOccCountTable(unsigned int *dnaDecodeTable);
void ForwardDNAAllOccCountNoLimit(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned int* dnaDecodeTable);

#endif


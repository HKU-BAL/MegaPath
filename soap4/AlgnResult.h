/*
 *
 *    AlgnResult.h
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

///////////////////////////////////////////////////////////
// Structures for storing single-end alignment result    //
///////////////////////////////////////////////////////////

#ifndef __ALIGN_RESULT_H__
#define __ALIGN_RESULT_H__


#define INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN 10485760 // 10M

#define INITIAL_SIZE_SA_LIST_FOR_SINGLE 10485760 // 10M
#define INITIAL_SIZE_OCC_LIST_FOR_SINGLE 10485760 // 10M
#define INITIAL_SIZE_READ_FOR_SINGLE 1048576 // 1M

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"

// To hold one occurrence.
typedef struct OccRecord
{
    unsigned long long pos;
    // char paired;          // 0 = not paired, 1 = paired
    unsigned int posDiff; // same as posDiff in SARecord
    unsigned int seedLength;       // for statistic, seed alignment length
    // uint multiplicity;
    char strand;          // value 1: '+'
} OccRecord;

// To hold one SA range.
typedef struct SARecord
{
    unsigned long long saLeft;
    unsigned long long saRight;
    unsigned char strand; // value 1: '+'
    bool isNBMOnly;
    bool isExact;
    unsigned int posDiff; // for variable size seed, 
                          // most significant 16 bits represent the position difference of front seed position offset
                          // less significant 16 bits represent the position difference of latter seed position offset\
                          // 
                          //        Seed Alignment Example
                          //      Seed pos           Seed end
                          //         v                  v
                          // --------|XXXOOOOOOOOOOOXXXX|--------------
                          //             ^         ^
                          //         variable size position
                          // in this case, posDiff >> 16 = 3; posDiff & 0xffff = 4
    unsigned int seedLength;       // for variable size seed, maximun seed length is used
    unsigned int seedOffset;      
} SARecord;

// DP Info for each read
// used for deduplication
typedef struct DPInfoForReads
{
    unsigned long long * startPositions;
    unsigned short * strand_dpLengths;
    uint * peLeftAnchors;
    uint * peRightAnchors;
    uint size;
} DPInfoForReads;

DPInfoForReads * constructDPInfoForReads ( unsigned long long size );
void freeDPInfoForReads ( DPInfoForReads * dpInfoForReads );

// To contain the first (even) read id of the both-unaligned pairs
// (i.e. both ends cannot be aligned)
typedef struct BothUnalignedPairs
{
    unsigned int * readIDs; // contain the first read id of the unaligned pairs
    unsigned int totalNum; // the number of read id inside the array "readIDs".
    unsigned int size; // available size of the array
} BothUnalignedPairs;

typedef struct BothUnalignedPairsArrays
{
    BothUnalignedPairs ** array;
    unsigned int arrayNum; // number of arrays
} BothUnalignedPairsArrays;

typedef BothUnalignedPairs UnalignedSingles;
typedef BothUnalignedPairsArrays UnalignedSinglesArrays;


BothUnalignedPairs * constructBothUnalignedPairs ();
void freeBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs );
void resetBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArrays );
void addReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int readID );
// Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
void addAllFirstReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int totalReadNum );
BothUnalignedPairsArrays * constructBothUnalignedPairsArrays ( unsigned int numCPUThreads );
void freeBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArray );

#endif


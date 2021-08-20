/*
 *
 *    CPUfunctions.h
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

#ifndef _CPUFUNCTIONS_H_
#define _CPUFUNCTIONS_H_

#include <pthread.h>
#include <math.h>
#include <ctype.h>
#include <sys/mman.h>
#include "2bwt-lib/dictionary.h"
#include "2bwt-lib/iniparser.h"
#include "2bwt-lib/Timing.h"
#include "definitions.h"
#include "PE.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"
#include "SAM.h"
#include "IniParam.h"
#include "UsageInterface.h"
#include "IndexHandler.h"
#include "SeedPool.h"

// The offset mask to retrieve the least significant 24-bit from the 32-bit word.
#define  BGS_GPU_ANSWER_OFFSET_LENGTH   24
#define  BGS_GPU_ANSWER_OFFSET_MASK     ((1<<BGS_GPU_ANSWER_OFFSET_LENGTH)-1)

// for user does not specify # of mismatches and dp is disabled
// if the read length < 50, then return DEFAULT_NUM_MISMATCH_FOR_SHORT_READ
// else return DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ
int getDefaultMismatchNum ( uint read_length );

// get the max hit # for default DP
int getMaxHitNumForDefaultDP ( uint read_length );

// get the parameters for All DP
void getParameterForAllDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params );

// get the parameters for Default DP
void getParameterForNewDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 );

// get the parameters for Default DP
void getParameterForDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 );

// get the parameters for Deep DP
void getParameterForDeepDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2, int maxReadLength, int roundNum );

// get the parameters for Single-end DP
void getParameterForSingleDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length );
                                     
// pack the reads which are unpaired
void packUnPairedReads ( uint * queries, uint * readIDs, uint * readLengths, uint * unAlignedPair,
                         uint wordPerQuery, uint numOfUnPaired, ullint maxBatchSize );


// pack the reads with no alignment together
// return number of reads with no alignment
uint packReads ( uint * queries, uint * readIDs, uint * readLengths, uint * seedLengths, unsigned char * noAlignment,
                 uint wordPerQuery, uint numQueries );

// pack the reads with no alignment together
// return readIDS of the unaligned reads
uint * packReads2 ( uint * queries, uint * readLengths, unsigned char * noAlignment,
                    uint wordPerQuery, uint numQueries, uint & numUnAligned );


// repack the reads
// no read will be removed, but
// the reads which need to be processed in next-round by soap3 will be duplicated
// to the front of the list. The readIDs are stored inside the array called "needProcessPair"
// the corresponding readIDs inside "readInputForDP", "readInputForNewDP" and
// "bothUnalignedPairs" need to be updated correspondingly.
void repackUnPairedReads ( uint ** queries, uint ** readIDs, uint ** readLengths, uint * needProcessPair,
                           uint wordPerQuery, uint numOfReadsToProcess, ullint numOfTotalReads,
                           ReadInputForDPArrays * readInputForDPall,
                           ReadInputForDPArrays * readInputForNewDPall,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays );

// convert 2-bit-per-character array to 1-byte-per-character array
void retrieve2BWTQueryFormat ( unsigned int * packedPatterns,
                               unsigned int queryIdx, unsigned int wordPerQuery, unsigned char * unpackedPattern );

// For initialization of g_log_n array. This array will be used for calculation of MAPQ of
static inline void bwase_initialize ( int * g_log_n ) {
    for (int i = 1; i < 256; ++i) {
        g_log_n[i] = int(4.343 * log(i) + 0.5);
    }
}

#endif

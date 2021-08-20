/*
 *
 *    alignment.h
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

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <pthread.h>
#include "definitions.h"
#include "CPUfunctions.h"
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "AlgnResult.h"
#include "DV-SemiDP.h"
#include "DV-DPForBothUnalign.h"
#include "DV-DPForSingleReads.h"
#include "OutputDPResult.h"
#include "IndexHandler.h"
#include "SeedPool.h"

// Perform SOAP3-DP Paired-End Alignment
void soap3_dp_pair_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                           ullint maxBatchSize, uint numQueries, uint accumReadNum,
                           Soap3Index * index,
                           IniParams ini_params, InputOptions input_options,
                           uint maxReadLength, uint detected_read_length, uint detected_read_length2,
                           char * upkdQualities,
                           uint * readIDs, char ** queryNames, char ** queryComments,
                           DPInfoForReads * dpInfoForReads,
                           samfile_t ** currSamOutputFilePtr, samfile_t * samOutputDPFilePtr,
                           samfile_t * samOutputUnpairFilePtr,
                           unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays, SeedPool * seedPool,
                           double startTime, double & lastEventTime, double & totalAlignmentTime);

#endif

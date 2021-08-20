/*
 *
 *    OutputDPResult.h
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

#ifndef _OUTPUTDPRESULT_H_
#define _OUTPUTDPRESULT_H_

#include "PEAlgnmt.h"
#include "CPUfunctions.h"
#include "BGS-IO.h"
#include "definitions.h"
#include "SeedPool.h"

#define MC_CeilDivide16(x) ((x+15)>>4)

// output a set of deep DP alignment results
// the alignments of the same pair of reads have to be together
void outputDeepDPResult ( DeepDPAlignResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char ** queryNames, char ** queryComments,
                          unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                          samfile_t * samOutputDPFilePtr, HSP * hsp, int peStrandLeftLeg, int peStrandRightLeg );

// output best alignment result pointer
// pointer is pointing to data of algnResult, no free to the output
DeepDPAlignResult * outputDeepDPResult2 ( DeepDPAlignResult * algnResult, unsigned int algnNum,
                           unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                           unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                           samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg, SeedPool *seedPool,
                           OCC * tmpocc = NULL, DynamicUint8Array * tmpCharArray = NULL, char ** twoSamStrings = NULL, int threadId = 0 );


// output DP single-end alignment results
// the alignments of the same pair of reads have to be together
// algnNum: the number of results in "algnResult"
void outputDPSingleResult ( SingleAlgnmtResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char ** queryNames, char ** queryComments,
                            unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                            samfile_t * samOutputDPFilePtr, HSP * hsp );

void outputDPSingleResult2 ( SingleAlgnmtResult * algnResult, unsigned int algnNum,
                             unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                             unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                             samfile_t * samOutputDPFilePtr, Soap3Index * index );


// output the single-end alignment results for the pair-end reads
void outputSingleResultForPairEnds ( AllHits * allHits,
                                     unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                                     DPInfoForReads * dpInfoForReads,
                                     unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                                     samfile_t * samOutputDPFilePtr, Soap3Index * index );

#endif

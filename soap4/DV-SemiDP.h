/*
 *
 *    DV-SemiDP.h
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

#ifndef _SEMIDP_H_
#define _SEMIDP_H_

#include "PEAlgnmt.h"
#include "definitions.h"
#include "IndexHandler.h"
#include "SeedPool.h"

///////////////////////////////////////////////////////////
// The functions are for semi-global DP                  //
///////////////////////////////////////////////////////////

// To perform semi-global DP for single end DP alignment
void  semiGlobalDPForSingleEndDpAlignment ( SingleDP_Space::AlgnmtResultStream* singleDPResultStream, int insert_high, int insert_low,
                       unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                       DPInfoForReads * dpInfoForReads,
                       unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                       Soap3Index * index,
                       int alignmentType, DPParameters * dpParameters,
                       unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                       unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                       samfile_t ** samOutputDPFilePtrs,
                       BothUnalignedPairs * &unalignedReads );

// singleDPResultStream: input single end result for DP
// insert_high : maximum value of insert size
// insert_low: minimum value of insert size
// upkdQueries: the sequence of all reads. The first character of the read with ID "i"
//              is on the position "i*maxReadLength" of the array
// upkdQueryNames: the description of all reads.
// upkdReadLengths: the read lengths of all reads
// peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
// peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
// BWT: bwt structure with SA table inside
// HSP: hsp structure with packed sequence inside
// alignmentType: 1: All valid alignment; 4: random best alignment
// dpParameters: the parameters for DP
// Output:
// numDPAlignedRead: the total number of reads aligned by DP
// numDPAlignment: the total number of resulting alignments by DP
// accumReadNum: the accumulated number of reads being processed
// outputFormat: the format of output
// outputFile: the file pointer for outputing the results


#endif

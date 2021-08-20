/*
 *
 *    DV-DPForSingleReads.h
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


#ifndef __DP_FOR_SINGLE_READS_H__
#define __DP_FOR_SINGLE_READS_H__

#include "AlgnResult.h"
#include "definitions.h"
#include "IndexHandler.h"
#include "SeedPool.h"

//////////////////////////////////////////////////////////////////////////////////////
// The functions are for the DP designed for unaligned single reads                 //
//////////////////////////////////////////////////////////////////////////////////////

// To perform semi-global DP for unaligned single read
void DPForUnalignSingle2 ( UnalignedSinglesArrays * unalignedReads, SeedPool* seedPool,
                           unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                           DPInfoForReads * dpInfoForReads,
                           unsigned int maxReadLength,
                           Soap3Index * index,
                           int alignmentType, DPParameters * dpParameters,
                           unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                           unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                           SingleDP_Space::AlgnmtResultStream * &resultStream,
                           QueryIDStream * &unalignedSingleReads );

// To perform semi-global DP for unaligned single read
void DPForUnalignSingle ( UnalignedSinglesArrays * unalignedReads,
                          unsigned char * upkdQueries, char ** queryNames, char ** queryComments,
                          unsigned int * upkdReadLengths,
                          unsigned int maxReadLength,
                          Soap3Index * index,
                          int alignmentType, DPParameters * dpParameters,
                          unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                          unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr );

void outputUnpairedSingleReads ( SingleDP_Space::AlgnmtResultStream * singleDPResultStream, UnalignedSingles * unalignedSingleReads,
                                 uint * queries, uint * upkdReadLengths, uint * origReadIDs, 
                                 char ** queryNames, char ** queryComments, char * upkdQualities, int inputMaxReadLength,
                                 uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                                 Soap3Index * index,
                                 int alignmentType,
                                 uint         &        dpSAlignedRead,
                                 uint         &        dpSAlignment   );

// For deduplication
void DPForAlignSingle ( UnalignedSinglesArrays * unalignedReads,
                        unsigned int * queries, unsigned int * upkdReadLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                        DPInfoForReads * dpInfoForReads,
                        unsigned int maxReadLength,
                        Soap3Index * index,
                        unsigned int * _bwt, unsigned int * _revBwt,
                        unsigned int * _occ, unsigned int * _revOcc,
                        int alignmentType, DPParameters * dpParameters,
                        unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                        unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                        SingleDP_Space::AlgnmtResultStream * &resultStream,
                        QueryIDStream * &unalignedSingleReads );


#endif

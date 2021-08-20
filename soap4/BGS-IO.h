/*
 *
 *    BGS-IO.h
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

#ifndef __BGS_IO__
#define __BGS_IO__

#include <stdio.h>
#include <stdlib.h>
#include "PEAlgnmt.h"
#include "2bwt-flex/SRAArguments.h"
#include "SAM.h"
#include "PE.h"
#include "Release.h"

//Define the below parameter to output the alignment result(text position)
// on screen instead of writing into output file
//#define DEBUG_2BWT_OUTPUT_TO_SCREEN

//Define the below parameter to stop cache reported SA range and text position.
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_OUTPUT

//Define the below parameter to skip writing the alignment result(text position)
// into disk. This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_WRITE_FILE

//Define the below parameter to skip translating the alignment result(text position).
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_TRANSLATION

#define OCC_CONST_ALIGNMENT_HEADER 65535
#define OCC_CONST_NO_ALIGNMENT     65534

//Define below parameter to skip outputing the header for plain output
#define SKIP_PLAIN_HEADER

OCC * OCCConstruct();
void OCCReset ( OCC * occ );
void OCCFree ( OCC * occ );

void OCCWriteOutputHeader ( HSP * hsp, FILE * outFilePtr,
                                    unsigned int maxReadLength,
                                    unsigned int numOfReads,
                                    int outputFormat );

//////////////////////////
//    FOR SAM FORMAT    //
//////////////////////////

// -- start --

// For outputting DP pair-end alignment in SAM format

void pairDeepDPOutputSAMAPI ( SRAQueryInput * qInput, DeepDPAlignResult * algnResult,
                              DeepDPAlignResult * bestResult,
                              unsigned int start, unsigned int num,
                              unsigned char * query1, unsigned char * query2,
                              char * qualities1, char * qualities2,
                              int readlen1, int readlen2,
                              char * queryName1, char * queryName2,
                              char * queryComment1, char * queryComment2,
                              DynamicUint8Array * xazArray, char ** twoSamStrings = NULL,
                              int threadId = 0 );
// For outputting Deep DP pair-end alignment in SAM format

void unproperlypairDPOutputSAMAPI ( SRAQueryInput * qInput, unsigned int firstReadID, Algnmt * algn_list1, unsigned int secondReadID, Algnmt * algn_list2,
                                    int hitNum1, int hitNum2,
                                    unsigned char * query1, unsigned char * query2,
                                    char * qualities1, char * qualities2,
                                    int readlen1, int readlen2,
                                    char * queryName1, char * queryName2,
                                    char * queryComment1, char * queryComment2,
                                    DynamicUint8Array * xazArray, Algnmt ** bestAlgns );
// output the DP alignments which are not properly paired


#endif

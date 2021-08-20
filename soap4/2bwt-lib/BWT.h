/*

   BWT.h    BWT-Index

   This module contains an implementation of BWT-index for alphabet size = 4.
   The functions provided include:
    Load functions for loading BWT to memory;
    Core functions for accessing core Inverse Psi values;
    Search functions for searching patterns from text;
    Text retrieval functions for retrieving text from BWT.

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

#ifndef __BWT_H__
#define __BWT_H__

#include <stdint.h>
#include "TypeNLimit.h"
#include "MemManager.h"
#include "TextConverter.h"
#include "HSP.h"

#define BITS_PER_OCC_VALUE            16
#define OCC_VALUE_PER_WORD            2
#define OCC_INTERVAL                256
#define WORD_BETWEEN_OCC            16
#define OCC_INTERVAL_MAJOR            65536

#define SORT_ALL                    0
#define SORT_16_BIT                    1
#define SORT_NONE                    2

#define BUCKET_BIT                    16
#define NUM_BUCKET                    65536

#define MAX_APPROX_MATCH_ERROR    7
#define MAX_ARPROX_MATCH_LENGTH    64

#define BWTDP_MAX_SUBSTRING_LENGTH    512

#define ESTIMATED_OCC_DIFF            32    // 128 / 4
#define MAX_OCC_DIFF                128

static const char * ECOR1 = "GAATTCGAATTC";
const size_t BWT_HOLLOW_BETWEEN_METADATA_PAYLOAD = 3; // * sizeof(unsigned int)
const size_t OCC_HOLLOW_BETWEEN_METADATA_PAYLOAD = 3; // * sizeof(unsigned int)
const size_t SA_PAYLOAD_OFFSET = ( 1 + ALPHABET_SIZE + 1 );
const size_t BWT_PAYLOAD_OFFSET = ( 1 + ALPHABET_SIZE + 0 + BWT_HOLLOW_BETWEEN_METADATA_PAYLOAD );
const size_t OCC_PAYLOAD_OFFSET = ( 1 + ALPHABET_SIZE + 0 + OCC_HOLLOW_BETWEEN_METADATA_PAYLOAD );


typedef struct BWT {
    unsigned long long textLength;            // length of the text
    unsigned long long saInterval;            // interval between two SA values stored explicitly
    //unsigned long long inverseSaInterval;        // interval between two inverse SA stored explicitly
    unsigned long long inverseSa0;            // SA-1[0]
    unsigned long long *cumulativeFreq;        // cumulative frequency
    uint32_t *bwtCode;                // BWT code
    uint32_t *occValue;                // Occurrence values stored explicitly
    unsigned long long *occValueMajor;        // Occurrence values stored explicitly
    unsigned long long *saValue;                // SA values stored explicitly
    //unsigned long long *inverseSa;            // Inverse SA stored explicitly
    //unsigned long long *cachedSaIndex;        // Cached SA index
    //unsigned long long cachedSaIndexNumOfChar;    // Number of characters indexed in SA index range
    unsigned long long *saValueOnBoundary;    // Pre-calculated frequently referred data
    uint32_t *decodeTable;            // For decoding BWT by table lookup
    //unsigned long long decodeTableGenerated;    // == TRUE if decode table is generated on load and will be freed
    unsigned long long bwtSizeInWord;            // Temporary variable to hold the memory allocated
    unsigned long long occSizeInWord;            // Temporary variable to hold the memory allocated
    unsigned long long occMajorSizeInWord;    // Temporary variable to hold the memory allocated
    unsigned long long saValueSizeInWord;        // Temporary variable to hold the memory allocated
    //unsigned long long inverseSaSizeInWord;    // Temporary variable to hold the memory allocated
    //unsigned long long cachedSaIndexSizeInWord;    // Temporary variable to hold the memory allocated
    unsigned long long (*_bwtSaValue)(const BWT*, unsigned long long); //Function Pointer to SaValue function
} BWT;


// Load / unload functions
BWT *BWTCreate(MMPool *mmPool, const unsigned long long textLength, unsigned int *decodeTable);
BWT *BWTLoad(MMPool *mmPool, char isShareIndex, const char *bwtCodeFileName, const char *occValueFileName, const char *saValueFileName);

void BWTFree(MMPool *mmPool, BWT *bwt, char isShareIndex);

// Precalculate frequenctly accessed data
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt);

// Core functions
// The following must be customized for differenet compression schemes ***
unsigned long long BWTDecode(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, const unsigned int character);
void BWTDecodeAll(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, unsigned long long* __restrict occValue);
unsigned long long BWTOccValue(const BWT *bwt, unsigned long long index, const unsigned int character);
void BWTAllOccValue(const BWT *bwt, unsigned long long index, unsigned long long* __restrict occValue);
unsigned long long BWTOccValueOnSpot(const BWT *bwt, unsigned long long index, unsigned int* __restrict character);
unsigned long long BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned long long searchOccValue);


// Utility functions for no compression only
unsigned long long BWTResidentSizeInWord(const unsigned long long numChar);
unsigned long long BWTFileSizeInWord(const unsigned long long numChar);
void BWTClearTrailingBwtCode(BWT *bwt);

// These are generic to different compression schemes (and generic to no compression as well)
unsigned long long BWTPsiMinusValue(const BWT *bwt, const unsigned long long index);
unsigned long long BWTPsiPlusValue(const BWT *bwt, const unsigned long long index);
unsigned long long BWTSaValue(const BWT *bwt, unsigned long long index);
unsigned long long BWTFullSaValue(const BWT *bwt, unsigned long long index);
unsigned long long BWTOccIntervalMajor(const unsigned long long occInterval);
unsigned long long BWTOccValueMinorSizeInWord(const unsigned long long numChar);
unsigned long long BWTOccValueMajorSizeInWord(const unsigned long long numChar);


#endif

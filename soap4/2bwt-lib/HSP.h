/*

   HSP.h        BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

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

#ifndef __HSP_H__
#define __HSP_H__

#include "TypeNLimit.h"
#include "MemManager.h"
#include "TextConverter.h"
#include <stdint.h>

#define ALPHABET_SIZE                4 
#define BIT_PER_CHAR                2
#define CHAR_PER_128                64
#define CHAR_PER_64                    32
#define CHAR_PER_WORD                16
#define CHAR_PER_BYTE                4
#define WORD_PER_LONG               2

#define CHAR_MASK               3

#define OUTPUT_PAIRWISE                0
#define OUTPUT_WITH_IDENTITIES        1
#define OUTPUT_NO_IDENTITIES        2
#define OUTPUT_FLAT_WITH_IDENTITIES    3
#define OUTPUT_FLAT_NO_IDENTITIES    4
#define OUTPUT_BLUNT_ENDS            5
#define OUTPUT_FLAT_BLUNT_ENDS        6
#define OUTPUT_XML                    7
#define OUTPUT_TABULAR                8
#define OUTPUT_TABULAR_COMMENT        9
#define OUTPUT_ASN_TEXT                10
#define OUTPUT_ASN_BINARY            11

#define QUERY_POS_STRAND    1
#define QUERY_NEG_STRAND    2
#define QUERY_BOTH_STRAND    3

#define CONTEXT_BIT                4
#define CONTEXT_MASK            0x0FFFFFFF

#define MAX_ALIGNMENT_LENGTH    131072

//=========================================

#define GRID_SAMPLING_FACTOR 262144
#define GRID_SAMPLING_FACTOR_2_POWER 18

typedef struct Translate{
    unsigned long long startPos;
    unsigned int chrID;
    unsigned long long correction;
} Translate;

typedef struct SeqActualOffset {
    unsigned long long startPos;
    unsigned long long endPos;
} SeqActualOffset;

typedef struct Annotation {
    int gi;  
    char text[MAX_SEQ_NAME_LENGTH+1];
    char decoratedText[MAX_SEQ_NAME_LENGTH+1];
} Annotation;

typedef struct SeqOffset {
    unsigned long long startPos;
    unsigned long long endPos;
    int firstAmbiguityIndex;    // The index for the first ambiguity that starts on or after the sequence
    int lastAmbiguityIndex;        // The index for the last ambiguity that ends on or before the sequence 
} SeqOffset;

typedef struct Ambiguity {
    unsigned long long startPos;
    unsigned long long rightOfEndPos;
    int symbol;
} Ambiguity;

typedef struct HSP {
    int numOfSeq;
    int numOfAmbiguity;
    unsigned long long dnaLength;
    unsigned int minSeqLength;
    unsigned int* packedDNA;
    Annotation* annotation;
    SeqOffset* seqOffset;
    Ambiguity* ambiguity;
    unsigned long long numOfRemovedSegment;
    unsigned long long numOfGridEntry;
    unsigned int* ambiguityMap;
    Translate* translate;
    SeqActualOffset* seqActualOffset;
} HSP;


#define MAX_SEQ_NAME_LENGTH                256

#define MAX_HISTO_SIZE                    256

#define INVALID_CHAR_INDEX                15

#define DP_MATCH_MISMATCH    0    // 0000
#define DP_INSERT            2    // 0010
#define DP_DELETE            3    // 0011
#define DP_INSERT_OPEN        4    // 0100
#define DP_DELETE_OPEN        8    // 1000

#define DP_MASK                3    // 0011

#define DP_INSERT_EXTEND    0    // 0000
#define DP_DELETE_EXTEND    0    // 0000

#define DP_NEG_INFINITY                -1073741824        // use -(2^30) to leave room for decreasing the value without overflow

#define ALIGN_MATCH                    0
#define ALIGN_MISMATCH_AMBIGUITY    1
#define ALIGN_INSERT                2
#define ALIGN_DELETE                3

#define ALIGN_PER_WORD                16
#define ALIGN_BIT                    2

#define AUX_TEXT_PER_WORD            8
#define AUX_TEXT_BIT                4

#define MAX_SP_SPACE_IN_GAP            3
#define MAX_SP_MISMATCH                3
#define MAX_SP_NON_ANCHOR_LENGTH    25    // MAX_SP_NON_ANCHOR_LENGTH + 2 * MAX_SP_SPACE_IN_GAP + 1 <= 32


static const char lowercaseDnaCharIndex = 14;    // Seems that BLAST treat masked characters as 'N' (still have 1/4 chance of matching)
static const char nonMatchDnaCharIndex  = 15;
static const char dnaChar[16]            = {'A', 'C', 'G', 'T', 'M', 'R', 'S', 'V', 'W', 'Y', 'H', 'K', 'D', 'B', 'N', 'L'};
static const char dnaComplement[16]        = {'T', 'G', 'C', 'A', 'K', 'Y', 'S', 'B', 'W', 'R', 'D', 'M', 'H', 'V', 'N', 'L'};
static const char ambiguityCount[16]    = { 1 ,  1 ,  1 ,  1 ,  2 ,  2 ,  2 ,  3 ,  2 ,  2 ,  3 ,  2 ,  3 ,  3 ,  4 ,  0 };

                                                        
HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName, const char *AmbiguityFileName, const char * TranslateFileName, const unsigned int trailerBufferInWord);
void HSPFree(MMPool *mmPool, HSP *hsp, const unsigned int trailerBufferInWord);

unsigned int HSPParseFASTAToPacked(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName,const char * translateFileName,
                      const unsigned int FASTARandomSeed, const int maskLowerCase);

#endif
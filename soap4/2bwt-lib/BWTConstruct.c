/*

   BWTConstruct.c        BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "QSufSort.h"
#include "r250.h"

// Static functions
static void BWTIncConstruct(BWTInc *bwtInc, const unsigned long long numChar);

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc);
static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned long long* __restrict rank, unsigned long long* __restrict cumulativeCount, const unsigned long long numChar);
static void BWTIncBuildPackedBwt(const unsigned long long *relativeRank, unsigned int* __restrict bwt, const unsigned long long numChar,
                                 const unsigned long long *cumulativeCount, const unsigned int *packedShift);
static unsigned long long  BWTIncGetAbsoluteRank(BWT *bwt, unsigned long long* __restrict absoluteRank, unsigned long long* __restrict seq, const unsigned int *packedText, 
                                  const unsigned long long numChar, const unsigned long long* cumulativeCount, const unsigned int firstCharInLastIteration);
static void BWTIncSortKey(unsigned long long* __restrict key, unsigned long long* __restrict seq, const unsigned long long numItem);
static void BWTIncBuildRelativeRank(unsigned long long* __restrict sortedRank, unsigned long long* __restrict seq, unsigned long long* __restrict relativeRank, const unsigned long long numItem, 
                                    unsigned long long oldInverseSa0, const unsigned long long *cumulativeCount);
static void BWTIncBuildBwt(unsigned long long*  seq, const unsigned long long *relativeRank, const unsigned long long numChar, const unsigned long long *cumulativeCount);    // seq is replaced with Bwt
static void BWTIncMergeBwt(const unsigned long long *sortedRank, const unsigned int* oldBwt, const unsigned long long *insertBwt, unsigned int* __restrict mergedBwt, 
                           const unsigned long long numOldBwt, const unsigned long long numInsertBwt);


BWTInc *BWTIncCreate(MMPool *mmPool, const unsigned long long textLength, const float targetNBit,
                     const unsigned long long initialMaxBuildSize, const unsigned long long incMaxBuildSize) {

    BWTInc *bwtInc;
    unsigned long long i;

    if (targetNBit == 0) {
        fprintf(stderr, "BWTIncCreate() : targetNBit = 0!\n");
        exit(1);
    }
    
    bwtInc = (BWTInc*) MMPoolDispatch(mmPool, sizeof(BWTInc));
    
    bwtInc->numberOfIterationDone = 0;

    bwtInc->bwt = BWTCreate(mmPool, textLength, NULL);

    bwtInc->initialMaxBuildSize = initialMaxBuildSize;
    bwtInc->incMaxBuildSize = incMaxBuildSize;

    bwtInc->targetNBit = targetNBit;

    bwtInc->cumulativeCountInCurrentBuild = (unsigned long long*) MMPoolDispatch(mmPool, sizeof(unsigned long long) * (ALPHABET_SIZE + 1));
    initializeLONG(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

    // Build frequently accessed data
    bwtInc->packedShift = (unsigned int*) MMPoolDispatch(mmPool, sizeof(unsigned int) * CHAR_PER_WORD);
    for (i=0; i<CHAR_PER_WORD; i++) {
        bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
    }

    bwtInc->targetTextLength = textLength;
    bwtInc->availableWord = (unsigned long long)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
    if (bwtInc->availableWord < BWTResidentSizeInWord(textLength) + WORD_BETWEEN_OCC / 2 + BWTOccValueMinorSizeInWord(textLength) * WORD_PER_LONG) {    // + 8 words so that the 128 bits before and after an explicit occ are in the same aligned 64 byte
        fprintf(stderr, "BWTIncCreate() : targetNBit is too low!\n");
        exit(1);
    }
    if (bwtInc->availableWord < BWTINC_MIN_MEMORY_IN_WORD) {
        bwtInc->availableWord = BWTINC_MIN_MEMORY_IN_WORD;
    }
    bwtInc->workingMemory = (unsigned int*) MMUnitAllocate(bwtInc->availableWord * BYTES_IN_WORD);

    return bwtInc;

}

void BWTIncFree(MMPool *mmPool, BWTInc *bwtInc) {

    MMUnitFree(bwtInc->workingMemory, bwtInc->availableWord * BYTES_IN_WORD);
    MMPoolReturn(mmPool, bwtInc, sizeof(BWTInc));

}

BWTInc *BWTIncConstructFromPacked(MMPool *mmPool, const char *inputFileName, const unsigned int showProgress,
                                  const float targetNBit, const unsigned long long initialMaxBuildSize, const unsigned long long incMaxBuildSize) {

    FILE *packedFile;
    unsigned long long packedFileLen;
    unsigned long long totalTextLength;
    unsigned long long textToLoad, textSizeInByte;
    unsigned long long processedTextLength;
    unsigned char lastByteLength;

    BWTInc *bwtInc;

    packedFile = (FILE*)fopen64(inputFileName, "rb");

    if (packedFile == NULL) {
        fprintf(stderr, "BWTIncConstructFromPacked() : Cannot open inputFileName!\n");
        exit(1);
    }

    fseek(packedFile, -1, SEEK_END);
    packedFileLen = ftell(packedFile);
    if ((long)packedFileLen < 0) {
        fprintf(stderr, "BWTIncConstructFromPacked: Cannot determine file length!\n");
        exit(1);
    }
    fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);
    totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

    bwtInc = BWTIncCreate(mmPool, totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

    BWTIncSetBuildSizeAndTextAddr(bwtInc);

    if (bwtInc->buildSize > totalTextLength) {
        textToLoad = totalTextLength;
    } else {
        textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
    }
    textSizeInByte = textToLoad / CHAR_PER_BYTE;    // excluded the odd byte

    fseek(packedFile, -2, SEEK_CUR);
    fseek(packedFile, -((long long)textSizeInByte), SEEK_CUR);
    fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte + 1, packedFile);
    fseek(packedFile, -((long long)textSizeInByte + 1), SEEK_CUR);

    ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
    BWTIncConstruct(bwtInc, textToLoad);

    processedTextLength = textToLoad;

    while (processedTextLength < totalTextLength) {
        textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
        if (textToLoad > totalTextLength - processedTextLength) {
            textToLoad = totalTextLength - processedTextLength;
        }
        textSizeInByte = textToLoad / CHAR_PER_BYTE;
        fseek(packedFile, -((long long)textSizeInByte), SEEK_CUR);
        fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte, packedFile);
        fseek(packedFile, -((long long)textSizeInByte), SEEK_CUR);
        ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
        BWTIncConstruct(bwtInc, textToLoad);
        processedTextLength += textToLoad;
        if (showProgress && bwtInc->numberOfIterationDone % 10 == 0) {
            printf("%u iterations done. %llu characters processed.\n", bwtInc->numberOfIterationDone, processedTextLength);
        }
    }

    return bwtInc;

}

static void BWTIncConstruct(BWTInc *bwtInc, const unsigned long long numChar) {

    unsigned long long i;
    unsigned long long mergedBwtSizeInWord, mergedOccSizeInWord;
    unsigned long long firstCharInThisIteration;

    unsigned long long *relativeRank, *seq, *sortedRank;
    uint32_t *insertBwt, *mergedBwt;
    unsigned long long newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

    #ifdef DEBUG
    if (numChar > bwtInc->buildSize) {
        fprintf(stderr, "BWTIncConstruct(): numChar > buildSize!\n");
        exit(1);
    }
    #endif

    mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
    mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

    initializeLONG(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

    if (bwtInc->bwt->textLength == 0) {        // Initial build

        // Set address
        seq = (unsigned long long*)bwtInc->workingMemory;
        relativeRank = seq + bwtInc->buildSize + 1;
        mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord - WORD_BETWEEN_OCC / 2;    // build in place

        BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

        firstCharInThisIteration = relativeRank[0];
        relativeRank[numChar] = 0;

        // Sort suffix
        QSufSortSuffixSort((long long*)relativeRank, (long long*)seq, (long long)numChar, (long long)ALPHABET_SIZE - 1, 0, FALSE);
        newInverseSa0 = relativeRank[0];

        // Clear BWT area
        initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

        // Build BWT
        BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

        // so that the cumulativeCount is not deducted
        bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

    } else {        // Incremental build

#ifdef DEBUG
        if (numChar / CHAR_PER_WORD * CHAR_PER_WORD != numChar) {
            fprintf(stderr, "BWTIncConstruct(): numChar must be multiple of char per word!\n");
            exit(1);
        }
#endif

        // Set address
        sortedRank = (unsigned long long*)bwtInc->workingMemory;
        seq = sortedRank + bwtInc->buildSize + 1;
        insertBwt = (uint32_t*)seq;
        relativeRank = seq + bwtInc->buildSize + 1;

        // Store the first character of this iteration
        firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

        // Count occurrence of input text
        ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
        // Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
        bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
        bwtInc->cumulativeCountInCurrentBuild[2] += bwtInc->cumulativeCountInCurrentBuild[1];
        bwtInc->cumulativeCountInCurrentBuild[3] += bwtInc->cumulativeCountInCurrentBuild[2];
        bwtInc->cumulativeCountInCurrentBuild[4] += bwtInc->cumulativeCountInCurrentBuild[3];

        // Get rank of new suffix among processed suffix
        // The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
        oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank(bwtInc->bwt, sortedRank, seq, bwtInc->packedText, 
                                                          numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

        // Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
        for (i=0; i<ALPHABET_SIZE; i++) {
            if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
                bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
                BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
            } else {
                if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
                    BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
                }
                if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
                    BWTIncSortKey(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
                }
            }
        }

        // build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
        // the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
        BWTIncBuildRelativeRank(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
        if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
            fprintf(stderr, "BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank!\n");
            exit(1);
        }
#endif

        // Sort suffix
        QSufSortSuffixSort((long long*)relativeRank, (long long*)seq, (long long)numChar, (long long)numChar, 1, TRUE);

        newInverseSa0RelativeRank = relativeRank[0];
        newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

        sortedRank[newInverseSa0RelativeRank] = 0;    // a special value so that this is skipped in the merged bwt

        // Build BWT
        BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

        // Merge BWT
        mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord - WORD_BETWEEN_OCC / 2
                    - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
                    // minus numberOfIteration * occInterval to create a buffer for merging
        BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, seq, mergedBwt, bwtInc->bwt->textLength, numChar);

    }

    // Build auxiliary structure and update info and pointers in BWT
    bwtInc->bwt->textLength += numChar;
    bwtInc->bwt->bwtCode = mergedBwt;
    bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
    bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
    if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
        fprintf(stderr, "BWTIncConstruct() : Not enough memory allocated!\n");
        exit(1);
    }

    bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

    BWTClearTrailingBwtCode(bwtInc->bwt);
    BWTGenerateOccValueFromBwt(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
                               bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

    bwtInc->bwt->inverseSa0 = newInverseSa0;
    
    bwtInc->bwt->cumulativeFreq[1] += bwtInc->cumulativeCountInCurrentBuild[1] - (bwtInc->firstCharInLastIteration <= 0);
    bwtInc->bwt->cumulativeFreq[2] += bwtInc->cumulativeCountInCurrentBuild[2] - (bwtInc->firstCharInLastIteration <= 1);
    bwtInc->bwt->cumulativeFreq[3] += bwtInc->cumulativeCountInCurrentBuild[3] - (bwtInc->firstCharInLastIteration <= 2);
    bwtInc->bwt->cumulativeFreq[4] += bwtInc->cumulativeCountInCurrentBuild[4] - (bwtInc->firstCharInLastIteration <= 3);

    bwtInc->firstCharInLastIteration = firstCharInThisIteration;

    // Set build size and text address for the next build
    BWTIncSetBuildSizeAndTextAddr(bwtInc);
    bwtInc->numberOfIterationDone++;

}

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc) {

    unsigned long long maxBuildSize;

    if (bwtInc->bwt->textLength == 0) {
        // initial build
        // Minus 2 because n+1 entries of seq and rank needed for n char
        maxBuildSize = (bwtInc->availableWord - 2 * CHAR_PER_WORD - OCC_INTERVAL / CHAR_PER_WORD - WORD_BETWEEN_OCC / 2)
                            / (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD / 2;
        if (bwtInc->initialMaxBuildSize > 0) {
            bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
        } else {
            bwtInc->buildSize = maxBuildSize;
        }
    } else {
        // Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
        // Minus numberOfIterationDone because bwt slightly shift to left in each iteration
        maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3 * CHAR_PER_WORD
                             - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR - WORD_BETWEEN_OCC / 2) 
                             / 3 / 2;
        if (maxBuildSize < CHAR_PER_WORD) {
            fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
            exit(1);
        }
        if (bwtInc->incMaxBuildSize > 0) {
            bwtInc->buildSize = min(bwtInc->incMaxBuildSize, maxBuildSize);
        } else {
            bwtInc->buildSize = maxBuildSize;
        }
        if (bwtInc->buildSize < CHAR_PER_WORD) {
            bwtInc->buildSize = CHAR_PER_WORD;
        }
    }

    if (bwtInc->buildSize < CHAR_PER_WORD) {
        fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
        exit(1);
    }

    bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

    bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + CHAR_PER_WORD)*2;
    bwtInc->textBuffer = (unsigned char*)(bwtInc->workingMemory + bwtInc->buildSize * 2 + CHAR_PER_WORD * 2);

}

static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned long long* __restrict rank, unsigned long long* __restrict cumulativeCount, const unsigned long long numChar) {

    unsigned long long i, j;
    unsigned long long packedMask;
    unsigned long long rankIndex;
    unsigned long long lastWord, numCharInLastWord;

    unsigned char ALIGN_16 temp[CHAR_PER_WORD];

    __m128i p1, p2, mask;

    lastWord = (numChar - 1) / CHAR_PER_WORD;
    numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;
    //DEBUG printf("numCharInLastWord=%u\n",numCharInLastWord);

    packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
    rankIndex = numChar - 1;
    // Unpack word-packed text; temp[0] will be character in the least significant 2 bits
    p1 = _mm_cvtsi32_si128(packedText[lastWord]);
    p2 = _mm_srli_epi32(p1, 4);
    p1 = _mm_unpacklo_epi8(p1, p2);

    mask = _mm_set1_epi32(0x03030303);

    p2 = _mm_srli_epi32(p1, 2);
    p1 = _mm_unpacklo_epi8(p1, p2);

    p1 = _mm_and_si128(p1, mask);
    _mm_store_si128((__m128i*)temp, p1);

    for (i=CHAR_PER_WORD - numCharInLastWord; i<CHAR_PER_WORD; i++) {
    //DEBUG     printf("N temp[%u]=%u\n",i,temp[i]);
        cumulativeCount[temp[i]+1]++;
        rank[rankIndex] = temp[i];
        rankIndex--;
    }

    for (i=lastWord; i--;) {    // loop from lastWord - 1 to 0

        // Unpack word-packed text; temp[0] will be character in the least significant 2 bits
        p1 = _mm_cvtsi32_si128(packedText[i]);
        p2 = _mm_srli_epi32(p1, 4);
        p1 = _mm_unpacklo_epi8(p1, p2);

        mask = _mm_set1_epi32(0x03030303);

        p2 = _mm_srli_epi32(p1, 2);
        p1 = _mm_unpacklo_epi8(p1, p2);

        p1 = _mm_and_si128(p1, mask);
        _mm_store_si128((__m128i*)temp, p1);

        for (j=0; j<CHAR_PER_WORD; j++) {
        //DEBUG     printf("temp[%u]=%u\n",j,temp[j]);
            cumulativeCount[temp[j]+1]++;
            rank[rankIndex] = temp[j];
            rankIndex--;
        }
    }

    // Convert occurrence to cumulativeCount
    cumulativeCount[2] += cumulativeCount[1];
    cumulativeCount[3] += cumulativeCount[2];
    cumulativeCount[4] += cumulativeCount[3];
}

static void BWTIncBuildPackedBwt(const unsigned long long *relativeRank, unsigned int* __restrict bwt, const unsigned long long numChar,
                                 const unsigned long long *cumulativeCount, const unsigned int *packedShift) {

    unsigned long long i, c, r;
    unsigned long long previousRank, currentRank;
    unsigned long long wordIndex, charIndex;
    unsigned long long inverseSa0;

    inverseSa0 = previousRank = relativeRank[0];

    for (i=1; i<=numChar; i++) {
        currentRank = relativeRank[i];
        // previousRank > cumulativeCount[c] because $ is one of the char
        c = (previousRank > cumulativeCount[1]) + (previousRank > cumulativeCount[2]) 
                                                + (previousRank > cumulativeCount[3]);
        // set bwt for currentRank
        if (c > 0) {
            // c <> 'a'
            r = currentRank;
            if (r > inverseSa0) {
                // - 1 because $ at inverseSa0 is not encoded            
                r--;
            }
            wordIndex = r / CHAR_PER_WORD;
            charIndex = r - wordIndex * CHAR_PER_WORD;
            bwt[wordIndex] |= c << packedShift[charIndex];
        }
        previousRank = currentRank;
    }

}

static unsigned long long BWTIncGetAbsoluteRank(BWT *bwt, unsigned long long* __restrict absoluteRank, unsigned long long* __restrict seq, const unsigned int *packedText, 
                                  const unsigned long long numChar, const unsigned long long* cumulativeCount, const unsigned int firstCharInLastIteration) {

    unsigned long long saIndex;
    unsigned long long lastWord;
    unsigned int packedMask;
    unsigned long long i, j;
    unsigned int c, t;
    unsigned long long rankIndex;
    unsigned long long shift;
    unsigned long long seqIndexFromStart[ALPHABET_SIZE];
    unsigned long long seqIndexFromEnd[ALPHABET_SIZE];

    for (i=0; i<ALPHABET_SIZE; i++) {
        seqIndexFromStart[i] = cumulativeCount[i];
        seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
    }

    shift = BITS_IN_WORD - BIT_PER_CHAR;
    packedMask = ALL_ONE_MASK >> shift;
    saIndex = bwt->inverseSa0;
    rankIndex = numChar - 1;

    lastWord = numChar / CHAR_PER_WORD;
    for (i=lastWord; i--;) {    // loop from lastWord - 1 to 0
        t = packedText[i];
        for (j=0; j<CHAR_PER_WORD; j++) {
            c = t & packedMask;
#ifdef DEBUG
            if (c >= ALPHABET_SIZE) {
                fprintf(stderr, "BWTIncGetAbsoluteRank() : c >= ALPHABET_SIZE!\n");
                exit(1);
            }
#endif
            saIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndex, c) + 1;
#ifdef DEBUG
            if (saIndex > bwt->textLength + 1) {
                fprintf(stderr, "BWTIncGetAbsoluteRank() : saIndex > bwt->textLength + 1!\n");
                exit(1);
            }
#endif
            // A counting sort using the first character of suffix is done here
            // If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
            if (saIndex > bwt->inverseSa0) {
                seq[seqIndexFromEnd[c]] = rankIndex;
                absoluteRank[seqIndexFromEnd[c]] = saIndex;
                seqIndexFromEnd[c]--;
            } else {
                seq[seqIndexFromStart[c]] = rankIndex;
                absoluteRank[seqIndexFromStart[c]] = saIndex;
                seqIndexFromStart[c]++;
            }
            rankIndex--;
            t >>= BIT_PER_CHAR;
        }
    }

    absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;    // representing the substring of all preceding characters
    seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

#ifdef DEBUG
    if (seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]) {
        fprintf(stderr, "BWTIncGetAbsoluteRank(): seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]!\n");
    }
#endif

    return seqIndexFromStart[firstCharInLastIteration];

}

static void BWTIncSortKey(unsigned long long* __restrict key, unsigned long long* __restrict seq, const unsigned long long numItem) {

    #define EQUAL_KEY_THRESHOLD    4    // Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

    int lowIndex, highIndex, midIndex;
    int lowPartitionIndex, highPartitionIndex;
    int lowStack[32], highStack[32];
    int stackDepth;
    int i, j;
    unsigned long long tempSeq, tempKey;
    int numberOfEqualKey;

    if (numItem < 2) {
        return;
    }

    stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

    for (;;) {

        for (;;) {

            // Sort small array of data
            if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {     // Insertion sort on smallest arrays
                for (i=lowIndex+1; i<=highIndex; i++) {
                    tempSeq = seq[i];
                    tempKey = key[i];
                    for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
                        seq[j] = seq[j-1];
                        key[j] = key[j-1];
                    }
                    if (j != i) {
                        seq[j] = tempSeq;
                        key[j] = tempKey;
                    }
                }
                break;
            }

            // Choose pivot as median of the lowest, middle, and highest data; sort the three data

            midIndex = average(lowIndex, highIndex);
            if (key[lowIndex] > key[midIndex]) {
                tempSeq = seq[lowIndex];
                tempKey = key[lowIndex];
                seq[lowIndex] = seq[midIndex];
                key[lowIndex] = key[midIndex];
                seq[midIndex] = tempSeq;
                key[midIndex] = tempKey;
            }
            if (key[lowIndex] > key[highIndex]) {
                tempSeq = seq[lowIndex];
                tempKey = key[lowIndex];
                seq[lowIndex] = seq[highIndex];
                key[lowIndex] = key[highIndex];
                seq[highIndex] = tempSeq;
                key[highIndex] = tempKey;
            }
            if (key[midIndex] > key[highIndex]) {
                tempSeq = seq[midIndex];
                tempKey = key[midIndex];
                seq[midIndex] = seq[highIndex];
                key[midIndex] = key[highIndex];
                seq[highIndex] = tempSeq;
                key[highIndex] = tempKey;
            }

            // Partition data

            numberOfEqualKey = 0;

            lowPartitionIndex = lowIndex + 1;
            highPartitionIndex = highIndex - 1;

            for (;;) {
                while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
                    numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
                    lowPartitionIndex++;
                }
                while (lowPartitionIndex < highPartitionIndex) {
                    if (key[midIndex] >= key[highPartitionIndex]) {
                        numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
                        break;
                    }
                    highPartitionIndex--;
                }
                if (lowPartitionIndex >= highPartitionIndex) {
                    break;
                }
                tempSeq = seq[lowPartitionIndex];
                tempKey = key[lowPartitionIndex];
                seq[lowPartitionIndex] = seq[highPartitionIndex];
                key[lowPartitionIndex] = key[highPartitionIndex];
                seq[highPartitionIndex] = tempSeq;
                key[highPartitionIndex] = tempKey;
                if (highPartitionIndex == midIndex) {
                    // partition key has been moved
                    midIndex = lowPartitionIndex;
                }
                lowPartitionIndex++;
                highPartitionIndex--;
            }

            // Adjust the partition index
            highPartitionIndex = lowPartitionIndex;
            lowPartitionIndex--;

            // move the partition key to end of low partition
            tempSeq = seq[midIndex];
            tempKey = key[midIndex];
            seq[midIndex] = seq[lowPartitionIndex];
            key[midIndex] = key[lowPartitionIndex];
            seq[lowPartitionIndex] = tempSeq;
            key[lowPartitionIndex] = tempKey;

            if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

                // Many keys = partition key; separate the equal key data from the lower partition
        
                midIndex = lowIndex;

                for (;;) {
                    while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
                        midIndex++;
                    }
                    while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
                        lowPartitionIndex--;
                    }
                    if (midIndex >= lowPartitionIndex) {
                        break;
                    }
                    tempSeq = seq[midIndex];
                    tempKey = key[midIndex];
                    seq[midIndex] = seq[lowPartitionIndex - 1];
                    key[midIndex] = key[lowPartitionIndex - 1];
                    seq[lowPartitionIndex - 1] = tempSeq;
                    key[lowPartitionIndex - 1] = tempKey;
                    midIndex++;
                    lowPartitionIndex--;
                }

            }

            if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
                // put the larger partition to stack
                lowStack[stackDepth] = lowIndex;
                highStack[stackDepth] = lowPartitionIndex - 1;
                stackDepth++;
                // sort the smaller partition first
                lowIndex = highPartitionIndex;
            } else {
                // put the larger partition to stack
                lowStack[stackDepth] = highPartitionIndex;
                highStack[stackDepth] = highIndex;
                stackDepth++;
                // sort the smaller partition first
                if (lowPartitionIndex > lowIndex) {
                    highIndex = lowPartitionIndex - 1;
                } else {
                    // all keys in the partition equals to the partition key
                    break;
                }
            }
            continue;

        }

        // Pop a range from stack
        if (stackDepth > 0) {
            stackDepth--;
            lowIndex = lowStack[stackDepth];
            highIndex = highStack[stackDepth];
            continue;
        } else {
            return;
        }

    }


}

static void BWTIncBuildRelativeRank(unsigned long long* __restrict sortedRank, unsigned long long* __restrict seq, unsigned long long* __restrict relativeRank, const unsigned long long numItem, 
                                    unsigned long long oldInverseSa0, const unsigned long long *cumulativeCount) {

    unsigned long long i, c;
    unsigned long long s, r;
    unsigned long long lastRank, lastIndex;
    unsigned long long oldInverseSa0RelativeRank = 0;
    unsigned long long freq;

    lastIndex = numItem;
    lastRank = sortedRank[numItem];
    if (lastRank > oldInverseSa0) {
        sortedRank[numItem]--;    // to prepare for merging; $ is not encoded in bwt
    }
    s = seq[numItem];
    relativeRank[s] = numItem;
    if (lastRank == oldInverseSa0) {
        oldInverseSa0RelativeRank = numItem;
        oldInverseSa0++;    // so that this segment of code is not run again
        lastRank++;            // so that oldInverseSa0 become a sorted group with 1 item
    }

    c = ALPHABET_SIZE - 1;
    freq = cumulativeCount[c];

    for (i=numItem; i--;) {    // from numItem - 1 to 0
        r = sortedRank[i];
        if (r > oldInverseSa0) {
            sortedRank[i]--;    // to prepare for merging; $ is not encoded in bwt
        }
        s = seq[i];
        if (i < freq) {
            if (lastIndex >= freq) {
                lastRank++;    // to trigger the group across alphabet boundary to be split
            }
            c--;
            freq = cumulativeCount[c];
        }
        if (r == lastRank) {
            relativeRank[s] = lastIndex;
        } else {
            if (i == lastIndex - 1) {
                if (lastIndex < numItem && (long long)seq[lastIndex + 1] < 0) {
                    seq[lastIndex] = seq[lastIndex + 1] - 1;
                } else {
                    seq[lastIndex] = (unsigned long long)-1;
                }
            }
            lastIndex = i;
            lastRank = r;
            relativeRank[s] = i;
            if (r == oldInverseSa0) {
                oldInverseSa0RelativeRank = i;
                oldInverseSa0++;    // so that this segment of code is not run again
                lastRank++;            // so that oldInverseSa0 become a sorted group with 1 item
            }
        }
    }

}

static void BWTIncBuildBwt(unsigned long long*  seq, const unsigned long long *relativeRank, const unsigned long long numChar, const unsigned long long *cumulativeCount) {

    unsigned long long i, c;
    unsigned long long previousRank, currentRank;

    previousRank = relativeRank[0];

    for (i=1; i<=numChar; i++) {
        currentRank = relativeRank[i];
        c = (previousRank >= cumulativeCount[1]) + (previousRank >= cumulativeCount[2])
                                                   + (previousRank >= cumulativeCount[3]);
        seq[currentRank] = c;
        previousRank = currentRank;
    }

}

static void BWTIncMergeBwt(const unsigned long long *sortedRank, const unsigned int* oldBwt, const unsigned long long *insertBwt, unsigned int* __restrict mergedBwt, 
                           const unsigned long long numOldBwt, const unsigned long long numInsertBwt) {

    unsigned long long bitsInWordMinusBitPerChar;
    unsigned long long leftShift, rightShift;
    unsigned long long o;
    unsigned long long oIndex, iIndex, mIndex;
    unsigned long long mWord, mChar, oWord, oChar;
    unsigned long long  numInsert;
    unsigned int insertBwtConv;

    bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

    oIndex = 0;
    iIndex = 0;
    mIndex = 0;

    mWord = 0;
    mChar = 0;

    mergedBwt[0] = 0;    // this can be cleared as merged Bwt slightly shift to the left in each iteration

    while (oIndex < numOldBwt) {

        // copy from insertBwt
        while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
            if (sortedRank[iIndex] != 0) {    // special value to indicate that this is for new inverseSa0
                insertBwtConv = insertBwt[iIndex];
                mergedBwt[mWord] |= insertBwtConv << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
                mIndex++;
                mChar++;
                if (mChar == CHAR_PER_WORD) {
                    mChar = 0;
                    mWord++;
                    mergedBwt[mWord] = 0;    // no need to worry about crossing mergedBwt boundary
                }
            }
            iIndex++;
        }

        // Copy from oldBwt to mergedBwt
        if (iIndex <= numInsertBwt) {
            o = sortedRank[iIndex];
        } else {
            o = numOldBwt;
        }
        numInsert = o - oIndex;

        oWord = oIndex / CHAR_PER_WORD;
        oChar = oIndex - oWord * CHAR_PER_WORD;
        if (oChar > mChar) {
            leftShift = (oChar - mChar) * BIT_PER_CHAR;
            rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
            mergedBwt[mWord] = mergedBwt[mWord]
                                | (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
                                | (oldBwt[oWord+1] >> rightShift);
            oIndex += min(numInsert, CHAR_PER_WORD - mChar);
            while (o > oIndex) {
                oWord++;
                mWord++;
                mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
                oIndex += CHAR_PER_WORD;
            }
        } else if (oChar < mChar) {
            rightShift = (mChar - oChar) * BIT_PER_CHAR;
            leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
            mergedBwt[mWord] = mergedBwt[mWord] 
                                | (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
            oIndex += min(numInsert, CHAR_PER_WORD - mChar);
            while (o > oIndex) {
                oWord++;
                mWord++;
                mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
                oIndex += CHAR_PER_WORD;
            }
        } else { // oChar == mChar
            mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
            oIndex += min(numInsert, CHAR_PER_WORD - mChar);
            while (o > oIndex) {
                oWord++;
                mWord++;
                mergedBwt[mWord] = oldBwt[oWord];
                oIndex += CHAR_PER_WORD;
            }
        }
        oIndex = o;
        mIndex += numInsert;

        // Clear the trailing garbage in mergedBwt
        mWord = mIndex / CHAR_PER_WORD;
        mChar = mIndex - mWord * CHAR_PER_WORD;
        if (mChar == 0) {
            mergedBwt[mWord] = 0;
        } else {
            mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
        }

    }

    // copy from insertBwt
    while (iIndex <= numInsertBwt) {
        if (sortedRank[iIndex] != 0) {
            insertBwtConv = insertBwt[iIndex];
            mergedBwt[mWord] |= insertBwtConv << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
            mIndex++;
            mChar++;
            if (mChar == CHAR_PER_WORD) {
                mChar = 0;
                mWord++;
                mergedBwt[mWord] = 0;    // no need to worry about crossing mergedBwt boundary
            }
        }
        iIndex++;
    }

}

unsigned long long BWTGenerateOccValueToFileFromBwt(const char *bwtFileName, const char *occValueFileName, unsigned int*  decodeTable) {

    FILE *bwtFile, *occValueFile;
    unsigned int *bwt;
    unsigned long long bwtFileSizeInWord, bwtResidentSizeInWord;
    unsigned long long textLength;
    unsigned long long inverseSa0;
    unsigned int *occValue;
    unsigned long long *occValueMajor;
    unsigned long long occSizeInWord, occMajorSizeInWord;
    unsigned long long i;
    unsigned long long cumulativeFreq[ALPHABET_SIZE];

    bwtFile = (FILE*)fopen64(bwtFileName, "rb");
    if (bwtFile == NULL) {
        fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open BWT file!\n");
        exit(1);
    }

    fread(&inverseSa0, sizeof(unsigned long long), 1, bwtFile);

    fread(cumulativeFreq, sizeof(unsigned long long), ALPHABET_SIZE, bwtFile);
    textLength = cumulativeFreq[ALPHABET_SIZE - 1];

    bwtResidentSizeInWord = BWTResidentSizeInWord(textLength);
    bwt = (unsigned int*) MMUnitAllocate(bwtResidentSizeInWord * sizeof(unsigned int));
    bwtFileSizeInWord = BWTFileSizeInWord(textLength);
    fread(bwt, sizeof(unsigned int), bwtFileSizeInWord, bwtFile);
    fclose(bwtFile);
    for (i=bwtFileSizeInWord; i<bwtResidentSizeInWord; i++) {
        bwt[i] = 0;
    }

    // occValue File

    occValueFile = (FILE*)fopen64(occValueFileName, "wb");
    if (occValueFile == NULL) {
        fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open occ value file!\n");
        exit(1);
    }

    fwrite(&inverseSa0, sizeof(unsigned long long), 1, occValueFile);
    fwrite(cumulativeFreq, sizeof(unsigned long long), ALPHABET_SIZE, occValueFile);

    occSizeInWord = BWTOccValueMinorSizeInWord(textLength);
    occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
    occValue = (unsigned int*) MMUnitAllocate(occSizeInWord * sizeof(unsigned int));
    occValueMajor = (unsigned long long*) MMUnitAllocate(occMajorSizeInWord * sizeof(unsigned long long));

    if (decodeTable == NULL) {
        decodeTable = (unsigned int*) MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
        GenerateDNAOccCountTable(decodeTable);
    }

    BWTGenerateOccValueFromBwt(bwt, occValue, occValueMajor, textLength, decodeTable);

    fwrite(occValue, sizeof(unsigned int), occSizeInWord, occValueFile);
    fwrite(occValueMajor, sizeof(unsigned long long), occMajorSizeInWord, occValueFile);
    fclose(occValueFile);

    MMUnitFree(occValue, occSizeInWord * sizeof(unsigned int));
    MMUnitFree(occValueMajor, occMajorSizeInWord * sizeof(unsigned long long));
    MMUnitFree(bwt, bwtResidentSizeInWord * sizeof(unsigned int));
    MMUnitFree(decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(uint32_t));

    return textLength;

}

void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned long long* __restrict occValueMajor,
                                const unsigned long long textLength, const unsigned int*  decodeTable) {

    unsigned long long numberOfOccValueMajor, numberOfOccValue;
    unsigned long long wordBetweenOccValue;
    unsigned long long numberOfOccIntervalPerMajor;
    unsigned long long c;
    unsigned long long i, j;
    unsigned long long occMajorIndex;
    unsigned long long occIndex, bwtIndex;
    unsigned long long sum;
    unsigned int tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

    wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

    // Calculate occValue

    numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;                // Value at both end for bi-directional encoding
    numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
    numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

    tempOccValue0[0] = 0;
    tempOccValue0[1] = 0;
    tempOccValue0[2] = 0;
    tempOccValue0[3] = 0;
    occValueMajor[0] = 0;
    occValueMajor[1] = 0;
    occValueMajor[2] = 0;
    occValueMajor[3] = 0;

    occIndex = 0;
    bwtIndex = 0;
    for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

        for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

            sum = 0;
            tempOccValue1[0] = tempOccValue0[0];
            tempOccValue1[1] = tempOccValue0[1];
            tempOccValue1[2] = tempOccValue0[2];
            tempOccValue1[3] = tempOccValue0[3];

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum += decodeTable[c >> 16];
                sum += decodeTable[c & 0x0000FFFF];
                bwtIndex++;
            }
            if (!DNA_OCC_SUM_EXCEPTION(sum)) {
                tempOccValue1[0] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue1[1] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue1[2] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue1[3] += sum;
            } else {
                if (sum == 0x00000100) {
                    tempOccValue1[0] += 256;
                } else if (sum == 0x00010000) {
                    tempOccValue1[1] += 256;
                } else if (sum == 0x01000000) {
                    tempOccValue1[2] += 256;
                } else {
                    tempOccValue1[3] += 256;
                }
            }
            occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
            occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
            occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
            occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
            tempOccValue0[0] = tempOccValue1[0];
            tempOccValue0[1] = tempOccValue1[1];
            tempOccValue0[2] = tempOccValue1[2];
            tempOccValue0[3] = tempOccValue1[3];
            sum = 0;

            occIndex++;

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum += decodeTable[c >> 16];
                sum += decodeTable[c & 0x0000FFFF];
                bwtIndex++;
            }
            if (!DNA_OCC_SUM_EXCEPTION(sum)) {
                tempOccValue0[0] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue0[1] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue0[2] += (sum & 0x000000FF);    sum >>= 8;
                tempOccValue0[3] += sum;
            } else {
                if (sum == 0x00000100) {
                    tempOccValue0[0] += 256;
                } else if (sum == 0x00010000) {
                    tempOccValue0[1] += 256;
                } else if (sum == 0x01000000) {
                    tempOccValue0[2] += 256;
                } else {
                    tempOccValue0[3] += 256;
                }
            }
        }

        occValueMajor[occMajorIndex * 4 + 0] = occValueMajor[(occMajorIndex - 1) * 4 + 0] + tempOccValue0[0];
        occValueMajor[occMajorIndex * 4 + 1] = occValueMajor[(occMajorIndex - 1) * 4 + 1] + tempOccValue0[1];
        occValueMajor[occMajorIndex * 4 + 2] = occValueMajor[(occMajorIndex - 1) * 4 + 2] + tempOccValue0[2];
        occValueMajor[occMajorIndex * 4 + 3] = occValueMajor[(occMajorIndex - 1) * 4 + 3] + tempOccValue0[3];
        tempOccValue0[0] = 0;
        tempOccValue0[1] = 0;
        tempOccValue0[2] = 0;
        tempOccValue0[3] = 0;

    }

    while (occIndex < (numberOfOccValue-1)/2) {
        sum = 0;
        tempOccValue1[0] = tempOccValue0[0];
        tempOccValue1[1] = tempOccValue0[1];
        tempOccValue1[2] = tempOccValue0[2];
        tempOccValue1[3] = tempOccValue0[3];
        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum += decodeTable[c >> 16];
            sum += decodeTable[c & 0x0000FFFF];
            bwtIndex++;
        }
        if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            tempOccValue1[0] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[1] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[2] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[3] += sum;
        } else {
            if (sum == 0x00000100) {
                tempOccValue1[0] += 256;
            } else if (sum == 0x00010000) {
                tempOccValue1[1] += 256;
            } else if (sum == 0x01000000) {
                tempOccValue1[2] += 256;
            } else {
                tempOccValue1[3] += 256;
            }
        }
        occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
        occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
        occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
        occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
        tempOccValue0[0] = tempOccValue1[0];
        tempOccValue0[1] = tempOccValue1[1];
        tempOccValue0[2] = tempOccValue1[2];
        tempOccValue0[3] = tempOccValue1[3];
        sum = 0;
        occIndex++;

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum += decodeTable[c >> 16];
            sum += decodeTable[c & 0x0000FFFF];
            bwtIndex++;
        }
        if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            tempOccValue0[0] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue0[1] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue0[2] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue0[3] += sum;
        } else {
            if (sum == 0x00000100) {
                tempOccValue0[0] += 256;
            } else if (sum == 0x00010000) {
                tempOccValue0[1] += 256;
            } else if (sum == 0x01000000) {
                tempOccValue0[2] += 256;
            } else {
                tempOccValue0[3] += 256;
            }
        }
    }

    sum = 0;
    tempOccValue1[0] = tempOccValue0[0];
    tempOccValue1[1] = tempOccValue0[1];
    tempOccValue1[2] = tempOccValue0[2];
    tempOccValue1[3] = tempOccValue0[3];

    if (occIndex * 2 < numberOfOccValue - 1) {
        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum += decodeTable[c >> 16];
            sum += decodeTable[c & 0x0000FFFF];
            bwtIndex++;
        }
        if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            tempOccValue1[0] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[1] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[2] += (sum & 0x000000FF);    sum >>= 8;
            tempOccValue1[3] += sum;
        } else {
            if (sum == 0x00000100) {
                tempOccValue1[0] += 256;
            } else if (sum == 0x00010000) {
                tempOccValue1[1] += 256;
            } else if (sum == 0x01000000) {
                tempOccValue1[2] += 256;
            } else {
                tempOccValue1[3] += 256;
            }
        }
    }

    occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
    occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
    occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
    occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];

}

void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName) {

    FILE *bwtFile, *occValueFile;
    unsigned long long bwtLength;

    bwtFile = (FILE*)fopen64(bwtFileName, "wb");
    if (bwtFile == NULL) {
        fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open BWT code file!\n");
        exit(1);
    }

    fwrite(&bwt->inverseSa0, sizeof(unsigned long long), 1, bwtFile);
    fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned long long), ALPHABET_SIZE, bwtFile);
    bwtLength = BWTFileSizeInWord(bwt->textLength);
    fwrite(bwt->bwtCode, sizeof(unsigned int), bwtLength, bwtFile);
    fclose(bwtFile);

    occValueFile = (FILE*)fopen64(occValueFileName, "wb");
    if (occValueFile == NULL) {
        fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
        exit(1);
    }

    fwrite(&bwt->inverseSa0, sizeof(unsigned long long), 1, occValueFile);
    fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned long long), ALPHABET_SIZE, occValueFile);

    fwrite(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
    fwrite(bwt->occValueMajor, sizeof(unsigned long long), bwt->occMajorSizeInWord, occValueFile);
    fclose(occValueFile);

}

void BWTGenerateSaValue(BWT *bwt, const unsigned long long saValueFreq, unsigned int showProgressInterval) {

    unsigned long long saValue;
    unsigned long long saIndex;
    unsigned long long numberOfSaValue;
    unsigned long long numberOfSaValueGenerated;
    unsigned long long progressInterval;
    unsigned long long i;

    if (bwt->saInterval != ALL_ONE_MASK) {
        fprintf(stderr, "BWTGenerateSaValue() : saValue already exist!\n");
        exit(1);
    }

    //if (saValueFreq == 1) {
        //BWTGenerateFullSaValue(bwt);
        //return;
    //}
    saValue = bwt->textLength;
    saIndex = 0;

    numberOfSaValue = (bwt->textLength + saValueFreq) / saValueFreq;
    bwt->saValueSizeInWord = numberOfSaValue;
    bwt->saInterval = saValueFreq;

    bwt->saValue = (unsigned long long*) MMUnitAllocate(bwt->saValueSizeInWord * sizeof(unsigned long long));

    #ifdef DEBUG
    initializeVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK);
    #endif

    numberOfSaValueGenerated = 0;

    if (showProgressInterval == 0 || showProgressInterval > numberOfSaValue) {
        progressInterval = numberOfSaValue;
    } else {
        progressInterval = showProgressInterval;
    }

    while (numberOfSaValueGenerated < numberOfSaValue) {
        progressInterval = min(numberOfSaValue - numberOfSaValueGenerated, progressInterval);
        for (i=0; i<progressInterval; i++) {
            while (saIndex % saValueFreq != 0) {
                #ifdef DEBUG
                if (saValue == 0) {
                    fprintf(stderr, "BWTGenerateSaValue(): saValue < 0!\n");
                    exit(1);
                }
                #endif
                saValue--;
                saIndex = BWTPsiMinusValue(bwt, saIndex);
                #ifdef DEBUG
                if (saIndex > bwt->textLength) {
                    fprintf(stderr, "BWTGenerateSaValue() : saIndex > textLength!\n");
                    exit(1);
                }
                #endif
            }
            bwt->saValue[saIndex/saValueFreq] = saValue;
            saValue--;
            saIndex = BWTPsiMinusValue(bwt, saIndex);
        }
        numberOfSaValueGenerated += progressInterval;
        if (showProgressInterval > 0) {
            printf("SA Value generated : %llu\n", numberOfSaValueGenerated);
        }
    }

    #ifdef DEBUG
    if (numberOfMatchInVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK) > 0) {
        fprintf(stderr, "BWTGenerateSaValue() : some saValue is not filled!\n");
        exit(1);
    }
    #endif

    bwt->saValue[0] = (unsigned long long)-1;    // Special handling

}

void BWTGenerateFullSaValue(BWT *bwt) {

    unsigned long long i, j;
    unsigned long long c, t, n;
    unsigned long long freq[ALPHABET_SIZE];

    if (bwt->saInterval != ALL_ONE_MASK) {
        fprintf(stderr, "BWTGenerateSaValue() : saValue already exist!\n");
        exit(1);
    }

    bwt->saValueSizeInWord = bwt->textLength + 1;
    bwt->saInterval = 1;

    bwt->saValue = (unsigned long long*) MMUnitAllocate(bwt->saValueSizeInWord * sizeof(unsigned long long));

    for (i=0; i<ALPHABET_SIZE; i++) {
        freq[i] = 0;
    }

    n = 0;
    for (i=0; i<bwt->textLength / CHAR_PER_WORD; i++) {
        t = bwt->bwtCode[i];
        for (j=0; j<CHAR_PER_WORD; j++) {
            c = t >> (BITS_IN_WORD - BIT_PER_CHAR);
            t <<= BIT_PER_CHAR;
            freq[c]++;
            bwt->saValue[n + (n>=bwt->inverseSa0)] = freq[c] + bwt->cumulativeFreq[c];
            n++;
        }
    }
    t = bwt->bwtCode[i];
    for (j=0; j<bwt->textLength - i * CHAR_PER_WORD; j++) {
        c = t >> (BITS_IN_WORD - BIT_PER_CHAR);
        t <<= BIT_PER_CHAR;
        freq[c]++;
        bwt->saValue[n + (n>=bwt->inverseSa0)] = freq[c] + bwt->cumulativeFreq[c];
        n++;
    }

    t = 0;
    for (i=bwt->textLength; i>0; i--) {
        c = t;
        t = bwt->saValue[c];
        bwt->saValue[c] = i;
    }

    bwt->saValue[bwt->inverseSa0] = 0;

    bwt->saValue[0] = (unsigned long long)-1;    // Special handling

}

void BWTSaveSaValue(const BWT *bwt, const char *saValueFileName) {

    FILE *saValueFile;

    saValueFile = (FILE*)fopen64(saValueFileName, "wb");
    if (saValueFile == NULL) {
        fprintf(stderr, "BWTSaveSaValue() : cannot open saValueFile!\n");
        exit(1);
    }

    fwrite(&bwt->inverseSa0, sizeof(unsigned long long), 1, saValueFile);
    fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned long long), ALPHABET_SIZE, saValueFile);
    fwrite(&bwt->saInterval, sizeof(unsigned long long), 1, saValueFile);

    fwrite(&bwt->textLength, sizeof(unsigned long long), 1, saValueFile);    // Save SA values without special handling on SA[0]
    fwrite(bwt->saValue + 1, sizeof(unsigned long long), bwt->saValueSizeInWord - 1, saValueFile);

    fclose(saValueFile);

}





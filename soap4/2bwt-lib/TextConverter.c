/*

   TextConverter.c        Text Converter

   This module contains miscellaneous text conversion functions.

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
#include <sys/mman.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "r250.h"


unsigned int BitPerWordPackedChar(const unsigned int alphabetSize) {

    #ifdef DEBUG
    if (alphabetSize < 2) {
        fprintf(stderr, "BitPerWordPackedChar() : alphabetSize < 2!\n");
        exit(1);
    }
    #endif

    return ceilLog2(alphabetSize);

}


unsigned long long TextLengthFromWordPacked(unsigned long long wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength) {

    return (wordPackedLength - 1) * (BITS_IN_WORD / bitPerChar) + lastWordLength;

}



unsigned int BitPerBytePackedChar(const unsigned int alphabetSize) {

    unsigned int bitPerChar;

    #ifdef DEBUG
    if (alphabetSize < 2) {
        fprintf(stderr, "BitPerBytePackedChar() : alphabetSize < 2!\n");
        exit(1);
    }
    #endif

    bitPerChar = ceilLog2(alphabetSize);

    #ifdef DEBUG
    if (bitPerChar > BITS_IN_BYTE) {
        fprintf(stderr, "BitPerBytePackedChar() : bitPerChar > BITS_IN_BYTE!\n");
        exit(1);
    }
    #endif

    // Return the largest number of bit that does not affect packing efficiency
    if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar) {
        bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
    }
    return bitPerChar;
}

unsigned long long TextLengthFromBytePacked(unsigned long long bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength) {

    if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
        fprintf(stderr, "TextLengthFromBytePacked(): text length > 2^32!\n");
    }
    return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;

}

void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned int bitPerChar, charPerWord;
    unsigned long long i, j, k;
    unsigned int c;
    unsigned long long charValue;

    bitPerChar = BitPerWordPackedChar(alphabetSize);
    charPerWord = BITS_IN_WORD / bitPerChar;

    for (i=0; i<textLength/charPerWord; i++) {
        c = 0;
        j = i * charPerWord;
        for (k=0; k<charPerWord; k++) {
            charValue = charMap[input[j+k]];
            if (charValue >= alphabetSize) {
                charValue = 0;
            }
            c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
        }
        output[i] = c;
    }
    if (i * charPerWord < textLength) {
        c = 0;
        j = i * charPerWord;
        for (k=0; j+k < textLength; k++) {
            charValue = charMap[input[j+k]];
            if (charValue >= alphabetSize) {
                charValue = 0;
            }
            c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
        }
        output[i] = c;
    }

}

void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned int bitPerChar, charPerByte;
    unsigned long long i, j, k;
    unsigned char c;

    bitPerChar = BitPerBytePackedChar(alphabetSize);
    charPerByte = BITS_IN_BYTE / bitPerChar;

    for (i=0; i<textLength/charPerByte; i++) {
        c = 0;
        j = i * charPerByte;
        for (k=0; k<charPerByte; k++) {
            c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
        }
        output[i] = c;
    }
    if (i * charPerByte < textLength) {
        c = 0;
        j = i * charPerByte;
        for (k=0; j+k < textLength; k++) {
            c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
        }
        output[i] = c;
    }

}

void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned int bitPerChar, charPerWord;
    unsigned long long i, j, k;
    unsigned int c;

    bitPerChar = BitPerWordPackedChar(alphabetSize);
    charPerWord = BITS_IN_WORD / bitPerChar;

    for (i=0; i<textLength/charPerWord; i++) {
        c = input[i];
        j = i * charPerWord;
        for (k=0; k<charPerWord; k++) {
            output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
            c <<= bitPerChar;
        }
    }
    if (i * charPerWord < textLength) {
        c = input[i];
        j = i * charPerWord;
        for (k=0; j+k<textLength; k++) {
            output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
            c <<= bitPerChar;
        }
    }

}

void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned int bitPerChar, charPerByte;
    unsigned long long i, j, k;
    unsigned char c;

    bitPerChar = BitPerBytePackedChar(alphabetSize);
    charPerByte = BITS_IN_BYTE / bitPerChar;

    for (i=0; i<textLength/charPerByte; i++) {
        c = input[i];
        j = i * charPerByte;
        for (k=0; k<charPerByte; k++) {
            output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
            c <<= bitPerChar;
        }
    }
    if (i * charPerByte < textLength) {
        c = input[i];
        j = i * charPerByte;
        for (k=0; j+k<textLength; k++) {
            output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
            c <<= bitPerChar;
        }
    }

}

void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned int bitPerChar, charPerByte;
    unsigned long long i, j, k;
    unsigned char c;

    bitPerChar = BitPerBytePackedChar(alphabetSize);
    charPerByte = BITS_IN_BYTE / bitPerChar;

    for (i=0; i<textLength/charPerByte; i++) {
        c = input[i];
        j = i * charPerByte;
        for (k=0; k<charPerByte; k++) {
            output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
            c <<= bitPerChar;
        }
    }
    if (i * charPerByte < textLength) {
        c = input[i];
        j = i * charPerByte;
        for (k=0; j+k<textLength; k++) {
            output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
            c <<= bitPerChar;
        }
    }

}

void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned long long i, j, k;
    unsigned int c;
    unsigned long long bitPerBytePackedChar;
    unsigned long long bitPerWordPackedChar;
    unsigned int charPerWord;
    unsigned int charPerByte;
    unsigned long long bytePerIteration;
    unsigned long long byteProcessed = 0;
    unsigned long long wordProcessed = 0;
    unsigned long long mask, shift;
    
    unsigned long long buffer[BITS_IN_WORD];

    bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
    bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
    charPerWord = BITS_IN_WORD / bitPerBytePackedChar;
    charPerByte = BITS_IN_BYTE / bitPerWordPackedChar;

    bytePerIteration = charPerWord / charPerByte;
    mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
    shift = BITS_IN_WORD - bitPerWordPackedChar;

    while ((wordProcessed + 1) * charPerWord < textLength) {

        c = input[wordProcessed];
        for (i=0; i<charPerWord; i++) {
            buffer[i] = c >> shift;
            c <<= bitPerWordPackedChar;
        }
        wordProcessed++;

        k = 0;
        for (i=0; i<bytePerIteration; i++) {
            c = 0;
            for (j=0; j<charPerByte; j++) {
                c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
                k++;
            }
            output[byteProcessed] = (unsigned char)c;
            byteProcessed++;
        }

    }

    c = input[wordProcessed];
    for (i=0; i < textLength - wordProcessed * charPerWord; i++) {
        buffer[i] = c >> shift;
        c <<= bitPerWordPackedChar;
    }

    k = 0;
    while (byteProcessed * charPerByte < textLength) {
        c = 0;
        for (j=0; j < textLength - wordProcessed * charPerWord; j++) {
            c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
            k++;
        }
        output[byteProcessed] = (unsigned char)c;
        byteProcessed++;
    }

}

void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned long long textLength) {

    unsigned long long i, j, k;
    unsigned int c;
    unsigned long long bitPerBytePackedChar;
    unsigned long long bitPerWordPackedChar;
    unsigned int charPerWord;
    unsigned int charPerByte;
    unsigned long long bytePerIteration;
    unsigned long long byteProcessed = 0;
    unsigned long long wordProcessed = 0;
    unsigned long long mask, shift;
    
    unsigned long long buffer[BITS_IN_WORD];

    bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
    bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
    charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
    charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

    bytePerIteration = charPerWord / charPerByte;
    mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
    shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

    while ((wordProcessed + 1) * charPerWord < textLength) {

        k = 0;
        for (i=0; i<bytePerIteration; i++) {
            c = (unsigned int)input[byteProcessed] << shift;
            for (j=0; j<charPerByte; j++) {
                buffer[k] = c & mask;
                c <<= bitPerBytePackedChar;
                k++;
            }
            byteProcessed++;
        }

        c = 0;
        for (i=0; i<charPerWord; i++) {
            c |= buffer[i] >> bitPerWordPackedChar * i;
        }
        output[wordProcessed] = c;
        wordProcessed++;

    }

    k = 0;
    for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
        c = (unsigned int)input[byteProcessed] << shift;
        for (j=0; j<charPerByte; j++) {
            buffer[k] = c & mask;
            c <<= bitPerBytePackedChar;
            k++;
        }
        byteProcessed++;
    }

    c = 0;
    for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
        c |= buffer[i] >> bitPerWordPackedChar * i;
    }
    output[wordProcessed] = c;

}

void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned long long textLength) {

    unsigned long long i;

    for (i=0; i< textLength; i++) {
        output[i] = charMap[input[i]];
    }

}

void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned long long textLength) {

    unsigned long long i;

    for (i=0; i< textLength; i++) {
        output[i] = reverseCharMap[input[i]];
    }

}



void * DNAMapPacked ( const char * inputFileName, unsigned long long * textLength, const unsigned int trailerBufferInWord )
{
    unsigned int * packedText;
    FILE * inputFile = ( FILE * ) fopen64 ( inputFileName, "rb" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "DNALoadPacked() : Cannot open inputFileName!\n" );
        exit ( 1 );
    }

    fread ( textLength, sizeof ( unsigned long long ), 1, inputFile );
    unsigned long long wordToProcess = ( *textLength + 64 - 1 ) / 64 * 4 + ( ( ( trailerBufferInWord + 3 ) / 4 * 4 ) * 4 );
    fprintf ( stderr, "[DNAMapPacked] Using shared index in host memory\n" );
    int inputFileFD = fileno ( inputFile );
    packedText = ( unsigned int * ) mmap ( NULL, ( sizeof ( unsigned int ) * ( wordToProcess + DNAPACK_PAYLOAD_OFFSET_BYWORD ) ), PROT_READ, MAP_SHARED | MAP_NORESERVE, inputFileFD, 0 );
    int rtVal;
    rtVal = mlock ( packedText, ( sizeof ( unsigned int ) * ( wordToProcess + DNAPACK_PAYLOAD_OFFSET_BYWORD ) ) );

    if ( ( void * ) rtVal == MAP_FAILED )
    {
        fprintf ( stderr, "\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:" );
        perror ( "" );
        exit ( EXIT_FAILURE );
    }

    packedText += DNAPACK_PAYLOAD_OFFSET_BYWORD;
    fclose ( inputFile );
    return ( void * ) packedText;
}



// Alphabet size of DNA must be 4
void *DNALoadPacked(const char *inputFileName, unsigned long long *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord) {

    FILE *inputFile;
    unsigned char tempChar[4];
    unsigned int *packedText;
    unsigned long long packedFileLen;
    unsigned char lastByteLength;
    unsigned long long wordToProcess;
    unsigned long long i;
    unsigned long long trailerBufferIn128;

    trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

    inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");

    if (inputFile == NULL) {
        fprintf(stderr, "DNALoadPacked() : Cannot open inputFileName!\n");
        exit(1);
    }

    fseek(inputFile, -1, SEEK_END);
    packedFileLen = ftell(inputFile);
    if ((long)packedFileLen < 0) {
        fprintf(stderr, "DNALoadPacked(): Cannot determine file length!\n");
        exit(1);
    }
    fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

    *textLength = (packedFileLen - 1) * 4 + lastByteLength;

    wordToProcess = (*textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4;        // allocate multiple of 128 bit + trailer buffer

    packedText = (unsigned int*) MMUnitAllocate(wordToProcess * sizeof(unsigned int));
    for (i=(*textLength)/16; i<wordToProcess; i++) {
        packedText[i] = 0;
    }

    fseek(inputFile, 0, SEEK_SET);
    fread(packedText, 1, packedFileLen, inputFile);
    fclose(inputFile);

    if (convertToWordPacked) {

        for (i=0; i<wordToProcess; i++) {
    
            *(unsigned int*)tempChar = packedText[i];
            packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

        }

    }

    return (void*)packedText;

}

void DNAFreePacked(void* packedDNA, const unsigned long long textLength, const unsigned int trailerBufferInWord, char isShareIndex) {

    unsigned int trailerBufferIn128;

    trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

    if ( isShareIndex )
    {
        packedDNA = ( ( unsigned int * ) packedDNA - DNAPACK_PAYLOAD_OFFSET_BYWORD );
        munmap ( packedDNA, ( ( ( textLength + 64 - 1 ) / 64 * 4 + trailerBufferIn128 * 4 ) * sizeof ( unsigned int ) ) );
    }
    else
    {
    MMUnitFree(packedDNA, ((textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4) * sizeof(unsigned int));
    }

}

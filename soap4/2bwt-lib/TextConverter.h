/*

   TextConverter.h        Text Converter

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

#ifndef __TEXTCONVERTOR_H__
#define __TEXTCONVERTOR_H__

#include <stdint.h>
#include "TypeNLimit.h"
#include "MemManager.h"

#define INVALID_CHAR 0xFF
#define CHAR_MAP_SIZE 256
#define PACKED_BUFFER_SIZE            (PACKED_BUFFER_SIZE_IN_WORD * BYTES_IN_WORD)
#define PACKED_BUFFER_SIZE_IN_WORD    65536
#define MAX_SEQ_NAME_LENGTH            256
#define RANDOM_SUBSTITUTE            'R'

#define DNAPACK_HOLLOW_BETWEEN_METADATA_PAYLOAD 12
#define DNAPACK_PAYLOAD_OFFSET_BYWORD ((DNAPACK_HOLLOW_BETWEEN_METADATA_PAYLOAD / 4) + 1)


// Word packed text functions
unsigned int BitPerWordPackedChar(const unsigned int alphabetSize);
unsigned long long TextLengthFromWordPacked(unsigned long long wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength);
// Byte packed text functions
unsigned int BitPerBytePackedChar(const unsigned int alphabetSize);
unsigned long long TextLengthFromBytePacked(unsigned long long bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength);

// Conversion functions
void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int textLength);
void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int textLength);

// Full load function
void *DNALoadPacked(const char *inputFileName, unsigned long long *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord);
void DNAFreePacked(void* packedDna, const unsigned long long textLength, const unsigned int trailerBufferInWord, char isShareIndex );
void * DNAMapPacked ( const char * inputFileName, unsigned long long * textLength, const unsigned int trailerBufferInWord );

#endif

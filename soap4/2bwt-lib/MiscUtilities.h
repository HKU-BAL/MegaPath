/*

   MiscUtilities.h		Miscellaneous Utilities

   This module contains miscellaneous utility functions.


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

#ifndef __MISC_UTILITIES_H__
#define __MISC_UTILITIES_H__

#include "TypeNLimit.h"
#include "stdint.h"
#include "stdio.h"

#define init(variable)							variable = 0;	// this is for avoiding compiler warning
																// disable it if compiler becomes smarter!

#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
// alignBoundary must be power of 2
#define nextAlignedBoundary(offset, alignBoundary)	( ((offset) + (alignBoundary) - 1) & (- (alignBoundary)) )
#define lastAlignedBoundary(offset, alignBoundary)		( (offset) & (- (alignBoundary)) )
#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define med3Index(key, ia, ib, ic)				( key[ia]<key[ib] ? (key[ib]<key[ic] ? ib : key[ia]<key[ic] ? ic : ia) : (key[ib]>key[ic] ? ib : key[ia]>key[ic] ? ic : ia))
#define swap(a, b, t);							t = a; a = b; b = t;

void Dust(const unsigned int len, unsigned char *pattern, const unsigned int level, const unsigned int window, const unsigned int word);

int QSortUnsignedIntOrder(const void *data, const long long index1, const long long index2);
void QSort(void* __restrict data, const long long numData, const long long dataWidth, int (*QSortComp)(const void*, const long long, const long long) );

unsigned int checkDuplicate(int *input, const unsigned int numItem, const int minValue, const int maxValue, char* text);
unsigned int leadingZero(const unsigned int input);
unsigned int ceilLog2(const unsigned int input);
unsigned int floorLog2(const unsigned int input);
unsigned int power(const unsigned int base, const unsigned int power);
void formatVALAsBinary(const unsigned int input, char* output, unsigned int bitGroup);
unsigned int getRandomSeed();

void ConvertBytePackedDNAToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int textLength);


unsigned int reverseBit(unsigned int x);
void initializeLONG(unsigned long long *startAddr, const unsigned long long length, const unsigned long long initValue);
void initializeVAL(unsigned int *startAddr, const unsigned long long length, const unsigned int initValue);
void initializeCHAR(unsigned char *startAddr, const unsigned long long length, const unsigned char initValue);
unsigned long long numberOfMatchInVAL(unsigned int *startAddr, const unsigned long long length, const unsigned int searchValue);
unsigned long long numberOfMatchInVAL(unsigned long long *startAddr, const unsigned long long length, const unsigned long long searchValue);
unsigned long long numberOfMatchInCHAR(unsigned char *startAddr, const unsigned long long length, const unsigned char searchValue);

void bitCopyNoDestOffset(unsigned int *destinationAddress, const unsigned int *sourceAddress,
							int sourceBitOffset, int copyLengthInBit);
void bitCopyNoDestBitOffset(unsigned int *destinationAddress, int destinationWordOffset,
							const unsigned int *sourceAddress, int sourceWordOffset,
							int sourceBitOffset, int copyLengthInBit);
unsigned int bitCopy(unsigned int *destinationAddress, int destinationWordOffset, int destinationBitOffset,
			 const unsigned int *sourceAddress, int sourceBitOffset, int copyLengthInBit);

unsigned int nextPrime(const unsigned int number);
unsigned int popCount(const unsigned int bitVector);

#endif


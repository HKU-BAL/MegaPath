/*
 *
 *    DPfunctions.cu
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

#include "CPU_DPfunctions.h"
#include "OutputDPResult.h"
#include "CPU_DP.h"
#include <assert.h>

#ifdef USE_CPU_DP

#include <algorithm>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// SemiGlobalAligner  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

uint SemiGlobalAligner::estimateThreadSize ( int maxReadLength, int maxDNALength )
{
	uint tableSizeOfThread = 2 * maxDNALength * maxReadLength * sizeof ( short );
	uint otherSizeOfThread = MC_CeilDivide16 ( maxDNALength ) * sizeof ( uint ) + //_packedDNASequence
	                         MC_CeilDivide16 ( maxReadLength ) * sizeof ( uint ) + //_packedReadSequence
	                         sizeof ( uint ) + //_DNALengths
	                         sizeof ( uint ) + //_readLengths
	                         sizeof ( int ) * 2 + //_scores, _cutoffThresholds
	                         sizeof ( uint ) + //_startOffsets
	                         sizeof ( uint ) * 2 + //_clipLtSizes, Rt
	                         sizeof ( uint ) * 2 + //_anchorLeftLocs, Right
	                         sizeof ( uint ) + //_hitLocs
	                         sizeof ( uint ) + //_maxScoreCounts
	                         ( maxReadLength + maxDNALength ) * sizeof ( uchar ); //_pattern
	return tableSizeOfThread + otherSizeOfThread + 16; // 16 is padding size
}

SemiGlobalAligner::SemiGlobalAligner()
{
	int i = 0;
    blockConf[i] = 256;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 192;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 128;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 64;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 48;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 32;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 16;
    coefConf[i] = 3;
    ++i;
    blockConf[i] = 8;
    coefConf[i] = 2;
    ++i;
    blockConf[i] = 2;
    coefConf[i] = 1.5; // <- desperate

    n_conf = i;
}
/*
SemiGlobalAligner::SemiGlobalAligner()
{
	n_conf = 6;

	for ( int i = 0; i < 4; i++ )
	{
		blockConf[i] = 64 - 16 * i;
		coefConf[i] = 3;
	}

	blockConf[4] = 8;
	coefConf[4] = 2;
	blockConf[5] = 2;
	coefConf[5] = 1.5; // <- desperate
}
*/
int SemiGlobalAligner::tryAlloc ( size_t estimatedThreadSize, size_t numOfBlocks )
{
    // Assume CPU memory is always enough 
    return 0;
}

void SemiGlobalAligner::decideConfiguration (
    int maxReadLength, int maxDNALength,
    int & maxDPTableLength, int & numOfBlocks,
    int & patternLength, DPParameters & dpPara
)
{
#define DP_Effective_Region(l)  (l + l/2 + 8)
	size_t avail, total;
	// cudaMemGetInfo ( &avail, &total );
    avail = 6442450944; // 6GB
    total = 6442450944; // 6GB
	size_t availableMemory = avail * 88 / 100; // 88%
	size_t estimatedThreadSize      = estimateThreadSize ( maxReadLength, maxDNALength );
	size_t availableBlocks          = availableMemory / ( estimatedThreadSize * DP_THREADS_PER_BLOCK );
	size_t estimatedThreadSize_2    = estimateThreadSize ( maxReadLength, DP_Effective_Region ( maxReadLength ) );
	size_t availableBlocks_2        = availableMemory / ( estimatedThreadSize_2 * DP_THREADS_PER_BLOCK );
	int successFlag = 0;

	for ( int i = 0; i < n_conf; i++ )
	{
		if ( availableBlocks >= blockConf[i] )
		{
			//              printf("[%d] scheme 1, try allocate\n", i);
			if ( 0 == tryAlloc ( estimatedThreadSize, blockConf[i] ) )
			{
				successFlag = 1;
				alignmentScheme = 1;
				numOfBlocks = blockConf[i];
				maxDPTableLength = maxDNALength;
				break;
			}
		}

		if ( availableBlocks_2 >= blockConf[i] &&
		        maxDNALength >= ( int ) ( maxReadLength * coefConf[i] ) )
		{
			//              printf("[%d] scheme 2, try allocate\n", i);
			if ( 0 == tryAlloc ( estimatedThreadSize_2, blockConf[i] ) )
			{
				successFlag = 1;
				alignmentScheme = 2;
				numOfBlocks = blockConf[i];
				maxDPTableLength = DP_Effective_Region ( maxReadLength );
				break;
			}
		}
	}

	if ( !successFlag )
	{
		// configuration failed
		printf ( "[DPfunc] error: insufficient GPU memory, cannot perform DP\n" );
		exit ( -1 );
	}

	// check invalid configuration
	if ( numOfBlocks < 32 )
	{
		printf ( "[DPfunc] warning: insufficient GPU memory, performance might degrade\n" );
	}

	if ( dpPara.matchScore > 30 )
	{
		printf ( "[DPfunc] warning: MatchScore (set to %d) should not exceed 30\n", dpPara.matchScore );
	}

	fflush ( stdout );
	patternLength = maxDNALength + maxReadLength;
}

void SemiGlobalAligner::decideConfigurationForIndelRA (
    int maxReadLength, int maxDNALength,
    int & maxDPTableLength, int & numOfBlocks,
    int & patternLength, DPParameters & dpPara
)
{
#define DP_Effective_Region_RA(l)  (l * 9)
	size_t avail, total;
    //	cudaMemGetInfo ( &avail, &total );
    avail = 6442450944; // 6GB
    total = 6442450944; // 6GB
	size_t availableMemory = avail * 88 / 100; // 88%
	size_t estimatedThreadSize      = estimateThreadSize ( maxReadLength, maxDNALength );
	size_t availableBlocks          = availableMemory / ( estimatedThreadSize * DP_THREADS_PER_BLOCK );
	size_t estimatedThreadSize_2    = estimateThreadSize ( maxReadLength, DP_Effective_Region_RA ( maxReadLength ) );
	size_t availableBlocks_2        = availableMemory / ( estimatedThreadSize_2 * DP_THREADS_PER_BLOCK );
	int successFlag = 0;

	for ( int i = 0; i < n_conf; i++ )
	{
		if ( availableBlocks >= blockConf[i] )
		{
			//              printf("[%d] scheme 1, try allocate\n", i);
			if ( 0 == tryAlloc ( estimatedThreadSize, blockConf[i] ) )
			{
				successFlag = 1;
				alignmentScheme = 1;
				numOfBlocks = blockConf[i];
				maxDPTableLength = maxDNALength;
				break;
			}
		}

		if ( availableBlocks_2 >= blockConf[i] &&
		        maxDNALength >= ( int ) ( maxReadLength * coefConf[i] ) )
		{
			//              printf("[%d] scheme 2, try allocate\n", i);
			if ( 0 == tryAlloc ( estimatedThreadSize_2, blockConf[i] ) )
			{
				successFlag = 1;
				alignmentScheme = 2;
				numOfBlocks = blockConf[i];
				maxDPTableLength = DP_Effective_Region_RA ( maxReadLength );
				break;
			}
		}
	}

	if ( !successFlag )
	{
		// configuration failed
		printf ( "[DPfunc] error: insufficient GPU memory, cannot perform DP\n" );
		exit ( -1 );
	}

	// check invalid configuration
	if ( numOfBlocks < 32 )
	{
		printf ( "[DPfunc] warning: insufficient GPU memory, performance might degrade\n" );
	}

	if ( dpPara.matchScore > 30 )
	{
		printf ( "[DPfunc] warning: MatchScore (set to %d) should not exceed 30\n", dpPara.matchScore );
	}

	fflush ( stdout );
	patternLength = maxReadLength + maxDNALength;
}

void SemiGlobalAligner::init (
    int batchSize,
    int maxReadLength, int maxDNALength, int maxDPTableLength,
    DPParameters & dpPara
)
{
	MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, dpPara );
#ifdef __AVX2__
	_DPTable =  _mm_malloc ( ( maxDNALength + 2 ) * ( maxReadLength + 1 ) * sizeof ( __m256i ), sizeof(__m256i)  );
	_DPTempRow = (__m256i*) _mm_malloc ( ( maxReadLength + 1 ) * sizeof ( __m256i ), sizeof(__m256i)  );
	//_DPTempRow = (__m256i*) _DPTable + (maxDNALength+1) * (maxReadLength+1);
#else
	_DPTable =  _mm_malloc ( ( maxDNALength + 2 ) * ( maxReadLength + 1 ) * sizeof ( __m128i ), sizeof(__m128i)  );
	_DPTempRow =  (__m128i*)_mm_malloc ( ( maxReadLength + 1 ) * sizeof ( __m128i ), sizeof(__m128i)  );
	//_DPTempRow = (__m128i*) _DPTable + (maxDNALength+1) * (maxReadLength+1);
#endif
	// use last row of _DPTable as temp storage
}

void SemiGlobalAligner::performAlignment (
    uint * packedDNASequence, uint * DNALengths,
    uint * packedReadSequence, uint * readLengths,
    int * cutoffThresholds, int * scores, uint * hitLocs,
    uint * maxScoreCounts,
    uchar * pattern, int numOfThreads,
    uint * clipLtSizes, uint * clipRtSizes,
    uint * anchorLeftLocs, uint * anchorRightLocs )
{
    // double startTime = setStartTime ();
    // init inputs
    _packedDNASequence = packedDNASequence;
    _packedReadSequence = packedReadSequence;
    _DNALengths = DNALengths;
    _readLengths = readLengths;
    _cutoffThresholds = cutoffThresholds;
    _clipLtSizes = clipLtSizes;
    _clipRtSizes = clipRtSizes;
    _anchorLeftLocs = anchorLeftLocs;
    _anchorRightLocs = anchorRightLocs;

    // init output locations
    _scores = scores;
    _hitLocs = hitLocs;
    _pattern = pattern;
    _maxScoreCounts = maxScoreCounts;
    
    // triger threads
    callDP ( _packedDNASequence, _DNALengths, maxDNALength,
             _packedReadSequence, _readLengths, maxReadLength,
             _clipLtSizes[0], _clipRtSizes[0],
             _anchorLeftLocs, _anchorRightLocs,
             dpPara.matchScore, dpPara.mismatchScore,
             dpPara.openGapScore, dpPara.extendGapScore,
             _cutoffThresholds,
             _DPTable, _DPTempRow, numOfThreads,
             _scores, _hitLocs,
             _maxScoreCounts, _pattern );

    // printf ( "%f\n", getElapsedTime ( startTime ) );

    return ;
}

void SemiGlobalAligner::freeMemory()
{
	_mm_free ( _DPTable );
	_mm_free ( _DPTempRow );
	// no need to free _DPTempRow again
}

#endif

/*
 *
 *    AlgnResult.h
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

#include "AlgnResult.h"

///////////////////////////////////////////////////////////
// Structures for storing DP Infos of a read             //
///////////////////////////////////////////////////////////
DPInfoForReads * constructDPInfoForReads ( unsigned long long size )
{
    DPInfoForReads * dpInfoForReads = ( DPInfoForReads *  ) malloc ( sizeof ( DPInfoForReads ) );
    dpInfoForReads->startPositions = ( unsigned long long * ) malloc ( size * sizeof ( unsigned long long ) );
    dpInfoForReads->strand_dpLengths = ( unsigned short * ) malloc ( size * sizeof ( unsigned short ) );
    dpInfoForReads->peLeftAnchors = ( uint * ) malloc ( size * sizeof ( uint ) );
    dpInfoForReads->peRightAnchors = ( uint * ) malloc ( size * sizeof ( uint ) );
    dpInfoForReads->size = size;

    memset ( dpInfoForReads->startPositions, -1, dpInfoForReads->size * sizeof ( unsigned long long ) );
    memset ( dpInfoForReads->strand_dpLengths, -1, dpInfoForReads->size * sizeof ( unsigned short ) );
    memset ( dpInfoForReads->peLeftAnchors, -1, dpInfoForReads->size * sizeof ( uint ) );
    memset ( dpInfoForReads->peRightAnchors, -1, dpInfoForReads->size * sizeof ( uint ) );
    return dpInfoForReads;
}

void freeDPInfoForReads ( DPInfoForReads * dpInfoForReads )
{
    free ( dpInfoForReads->startPositions );
    free ( dpInfoForReads->strand_dpLengths );
    free ( dpInfoForReads->peLeftAnchors );
    free ( dpInfoForReads->peRightAnchors );
    free ( dpInfoForReads );
}


///////////////////////////////////////////////////////////
// Structures for storing single-end alignment result    //
///////////////////////////////////////////////////////////

BothUnalignedPairs * constructBothUnalignedPairs ()
{
    BothUnalignedPairs * bothUnalignedPairs = ( BothUnalignedPairs * ) malloc ( sizeof ( BothUnalignedPairs ) );
    bothUnalignedPairs->readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN );
    bothUnalignedPairs->totalNum = 0;
    bothUnalignedPairs->size = INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN;
    return bothUnalignedPairs;
}

void freeBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs )
{
    free ( bothUnalignedPairs->readIDs );
    free ( bothUnalignedPairs );
}

void resetBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArrays )
{
    for ( unsigned int i = 0; i < bothUnalignedPairsArrays->arrayNum; i++ )
    {
        bothUnalignedPairsArrays->array[i]->totalNum = 0;
    }
}

void addReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int readID )
{
    if ( bothUnalignedPairs->totalNum == bothUnalignedPairs->size )
    {
        // enlarge the arrays "readIDs" by at least double
        unsigned int new_size = bothUnalignedPairs->size * 2;
        unsigned int * new_readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        memcpy ( new_readIDs, bothUnalignedPairs->readIDs, sizeof ( unsigned int ) *bothUnalignedPairs->totalNum );
        free ( bothUnalignedPairs->readIDs );
        bothUnalignedPairs->size = new_size;
        bothUnalignedPairs->readIDs = new_readIDs;
    }
    bothUnalignedPairs->readIDs[bothUnalignedPairs->totalNum] = readID;
    bothUnalignedPairs->totalNum++;
}

// Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
// assume totalReadNum is an even number
void addAllFirstReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int totalReadNum )
{
    unsigned int newReadNum = totalReadNum / 2;

    if ( bothUnalignedPairs->totalNum + newReadNum > bothUnalignedPairs->size )
    {
        // enlarge the arrays such that it is enough for storing all the new read IDs
        unsigned int new_size = ( bothUnalignedPairs->totalNum + newReadNum ) * 1.2;
        unsigned int * new_readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        memcpy ( new_readIDs, bothUnalignedPairs->readIDs, sizeof ( unsigned int ) *bothUnalignedPairs->totalNum );
        free ( bothUnalignedPairs->readIDs );
        bothUnalignedPairs->readIDs = new_readIDs;
    }

    unsigned int i;

    for ( i = 0; i <= totalReadNum - 2; i += 2 )
    {
        bothUnalignedPairs->readIDs[bothUnalignedPairs->totalNum] = i;
        bothUnalignedPairs->totalNum++;
    }
}

BothUnalignedPairsArrays * constructBothUnalignedPairsArrays ( unsigned int numCPUThreads )
{
    BothUnalignedPairsArrays * newArray;
    newArray = ( BothUnalignedPairsArrays * ) malloc ( sizeof ( BothUnalignedPairsArrays ) );
    newArray->arrayNum = numCPUThreads;
    newArray->array = ( BothUnalignedPairs ** ) malloc ( sizeof ( BothUnalignedPairs * ) * numCPUThreads );

    for ( unsigned int i = 0; i < numCPUThreads; i++ )
    {
        newArray->array[i] = constructBothUnalignedPairs ();
    }

    return newArray;
}

void freeBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArray )
{
    if ( bothUnalignedPairsArray == NULL ) return;
    for ( unsigned int i = 0; i < bothUnalignedPairsArray->arrayNum; i++ )
    {
        freeBothUnalignedPairs ( bothUnalignedPairsArray->array[i] );
    }

    free ( bothUnalignedPairsArray->array );

    free ( bothUnalignedPairsArray );
}


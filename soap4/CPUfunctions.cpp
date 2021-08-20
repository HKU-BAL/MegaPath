/*
 *
 *    CPUfunctions.cpp
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

#include "CPUfunctions.h"

int getDefaultMismatchNum ( uint read_length )
{
    // for user does not specify # of mismatches and dp is disabled
    // if the read length < 50, then return DEFAULT_NUM_MISMATCH_FOR_SHORT_READ
    // else return DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ
    if ( read_length > 50 )
    { return DEFAULT_NUM_MISMATCH_NO_DP_NORMAL_READ; }
    else
    { return DEFAULT_NUM_MISMATCH_NO_DP_SHORT_READ; }
}

int getMaxHitNumForDefaultDP ( uint read_length )
{
    // get the max hit # for default DP
    if ( read_length > 50 )
    { return MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { return MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }
}

void getParameterForAllDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params )
{
    // get the parameters for All DP
    // not related to read length
    dp_params.matchScore = ini_params.Ini_MatchScore;
    dp_params.mismatchScore = ini_params.Ini_MismatchScore;
    dp_params.openGapScore = ini_params.Ini_GapOpenScore;
    dp_params.extendGapScore = ini_params.Ini_GapExtendScore;
    dp_params.numOfCPUThreads = ini_params.Ini_NumOfCpuThreads;
    dp_params.numOfCPUForSeeding = ini_params.Ini_NumOfCpuThreads;
}


void getParameterForDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 )
{
    int max_read_length = (read_length > read_length2) ? (read_length) : (read_length2);
    // get the parameters for Default DP
    // for paired-end reads

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 50 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }

    if ( read_length2 > 50 )
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }

    dp_params.isExtensionDP = 0;
    // the following parameters are not needed for default DP
    // seedLength, sampleDist, tailTrimLen, singleDPSeedNum and singleDPSeedPos[]
}

void getParameterForNewDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 )
{
    // get the parameters for Default DP
    // for paired-end reads

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;
    dp_params.isExtensionDP = 0;
}


void getParameterForDeepDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2, int maxReadLength, int roundNum )
{
    int max_read_length = (read_length > read_length2) ? (read_length) : (read_length2);
    max_read_length = ( maxReadLength > max_read_length ? maxReadLength : max_read_length );
    SeedingProperties seedingProperties;

    if ( max_read_length > 120 )
        { seedingProperties = ini_params.seedingPropertiesForLongReads[roundNum-1]; }
    else
        { seedingProperties = ini_params.seedingPropertiesForShortReads[roundNum-1]; }

    // get the parameters for Deep DP

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 50 )
    { dp_params.paramRead[0].maxHitNum = seedingProperties.maxNumberOfSeedHit; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_SHORT_READ; }

    if ( read_length2 > 50 )
    { dp_params.paramRead[1].maxHitNum = seedingProperties.maxNumberOfSeedHit; }
    else
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_SHORT_READ; }

    dp_params.tailTrimLen = TAIL_TRIM_SEEDING_DEEP_DP;
    dp_params.seedProperties = seedingProperties;
    dp_params.roundNum = roundNum-1;
    if ( roundNum == 1 && read_length > 120 )
    {
        dp_params.isExtensionDP = 0;
    }
    else 
    {
        dp_params.isExtensionDP = 0;
    }

    dp_params.mmpProperties = ini_params.mmpProperties;
    fprintf(stderr, "[MMConfig] nbm: %u restricted: %u\n", dp_params.seedProperties.numNBMismatch, dp_params.seedProperties.restrictedNBM );
    // the following parameters are not needed for deep DP
    // singleDPSeedNum and singleDPSeedPos[]
}

void getParameterForSingleDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length )
{
    int max_read_length = read_length;

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 300 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_V_LONG_READ; }
    else if ( read_length > 80 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_SHORT_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_V_SHORT_READ; }

    // update: for reads longer than 100, one more seed for every extra 100 bases
    if ( read_length > 100 )
    {
        dp_params.singleDPSeedNum = SEED_NUM_SINGLE_DP + read_length / 100;
    }
    else
    {
        dp_params.singleDPSeedNum = SEED_NUM_SINGLE_DP;
    }

    // the followings are obsoleted
    dp_params.singleDPSeedPos[0] = 0;
    int X;

    if ( read_length > 80 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_SHORT_READ; }
    else
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_V_SHORT_READ; }

    dp_params.singleDPSeedPos[2] = ( read_length - X ) * 0.5 - 1;

    /*
    if ( dp_params.singleDPSeedPos[2] > read_length - dp_params.paramRead[0].seedLength )
    { dp_params.singleDPSeedPos[2] = read_length - dp_params.paramRead[0].seedLength; }
    */

    dp_params.isExtensionDP = 0;
    dp_params.singleDPSeedPos[1] = ( dp_params.singleDPSeedPos[0] + dp_params.singleDPSeedPos[2] ) / 2 - 1;
    dp_params.tailTrimLen = 0;
    // the following parameters are not needed for single DP
    // sampleDist and tailTrimLen

    if ( ini_params.Ini_skipSingleEndDP ) {
        dp_params.singleDPSeedNum = 0; // i.e. skip all
    }
}



// pack the reads with no alignment together
// return number of reads with no alignment
uint packReads ( uint * queries, uint * readIDs, uint * readLengths, uint * seedLengths, unsigned char * noAlignment,
                 uint wordPerQuery, uint numQueries )
{
    ullint total_no_alignment = 0;

    for ( ullint readId = 0; readId < numQueries; ++readId )
    {
        if ( noAlignment[readId] == 1 )
        {
            // the read has no alignment
            if ( total_no_alignment < readId )
            {
                ullint srcQueryOffset = ( readId ) / 32 * 32 * wordPerQuery + readId % 32;
                ullint srcNewQueryOffset = ( total_no_alignment ) / 32 * 32 * wordPerQuery + total_no_alignment % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    queries[srcNewQueryOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                readLengths[total_no_alignment] = readLengths[readId];

                if ( seedLengths != NULL )
                { seedLengths[total_no_alignment] = seedLengths[readId]; }

                readIDs[total_no_alignment] = readIDs[readId];
            }

            total_no_alignment++;
        }
    }

    return total_no_alignment;
}

// pack the reads with no alignment together
// return readIDS of the unaligned reads
uint * packReads2 ( uint * queries, uint * readLengths, unsigned char * noAlignment,
                    uint wordPerQuery, uint numQueries, uint & numUnAligned )
{
    uint * readIDs = ( uint * ) malloc ( sizeof ( unsigned int ) * numQueries );
    ullint total_no_alignment = 0;

    for ( ullint readId = 0; readId < numQueries; ++readId )
    {
        if ( noAlignment[readId] == 1 )
        {
            // the read has no alignment
            if ( total_no_alignment < readId )
            {
                ullint srcQueryOffset = ( readId ) / 32 * 32 * wordPerQuery + readId % 32;
                ullint srcNewQueryOffset = ( total_no_alignment ) / 32 * 32 * wordPerQuery + total_no_alignment % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    queries[srcNewQueryOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                readLengths[total_no_alignment] = readLengths[readId];
            }

            readIDs[total_no_alignment] = readId + 1;
            total_no_alignment++;
            // printf("total_no_alignment = %u\n", total_no_alignment);
        }
    }

    numUnAligned = total_no_alignment;
    return readIDs;
}



// pack the reads which are unpaired
void packUnPairedReads ( uint * queries, uint * readIDs, uint * readLengths, uint * unAlignedPair,
                         uint wordPerQuery, uint numOfUnPaired, ullint maxBatchSize )
{
    ullint roundUp = ( numOfUnPaired + 31 ) / 32 * 32;
    ullint totalQueryLength = roundUp * wordPerQuery;
    uint * newqueries;
    uint * newreadLengths;
    uint * newreadIDs;

    if ( numOfUnPaired > 0 )
    {
        // copy the content to the new array
        newqueries = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
        newreadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        newreadIDs = ( uint * ) malloc ( roundUp * sizeof ( uint ) );

        for ( ullint i = 0; i < numOfUnPaired; i++ )
        {
            ullint readId = unAlignedPair[i];
            ullint srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            ullint srcNewQueryOffset = i / 32 * 32 * wordPerQuery + i % 32;

            for ( ullint j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[i] = readLengths[readId];
            newreadIDs[i] = readIDs[readId];
        }
    }

    // clear the content of the old array
    ullint maxRoundUp = ( maxBatchSize + 31 ) / 32 * 32;
    ullint maxTotalQueryLength = maxRoundUp * wordPerQuery;
    memset ( queries, 0, maxTotalQueryLength * sizeof ( uint ) );
    memset ( readLengths, 0, maxRoundUp * sizeof ( uint ) );
    memset ( readIDs, 0, maxRoundUp * sizeof ( uint ) );

    if ( numOfUnPaired > 0 )
    {
        // copy the content of the new array to the old array
        memcpy ( queries, newqueries, totalQueryLength * sizeof ( uint ) );
        memcpy ( readLengths, newreadLengths, roundUp * sizeof ( uint ) );
        memcpy ( readIDs, newreadIDs, roundUp * sizeof ( uint ) );
        // free the memory
        free ( newqueries );
        free ( newreadLengths );
        free ( newreadIDs );
    }
}



// repack the reads
// no read will be removed, but
// the reads which need to be processed in next-round by soap3 will be duplicated
// to the front of the list. The readIDs are stored inside the array called "needProcessPair"
// the corresponding readIDs inside "readInputForDP", "readInputForNewDP" and
// "bothUnalignedPairs" need to be updated correspondingly.
void repackUnPairedReads ( uint ** orgqueries, uint ** orgreadIDs, uint ** orgreadLengths, uint * needProcessPair,
                           uint wordPerQuery, uint numOfReadsToProcess, ullint numOfTotalReads,
                           ReadInputForDPArrays * readInputForDPall,
                           ReadInputForDPArrays * readInputForNewDPall,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays )
{
    uint * queries = ( *orgqueries );
    uint * readIDs = ( *orgreadIDs );
    uint * readLengths = ( *orgreadLengths );

    if ( numOfReadsToProcess > 0 )
    {
        ullint i, j;
        ullint roundUp = ( numOfReadsToProcess + numOfTotalReads + 31 ) / 32 * 32;
        ullint totalQueryLength = roundUp * wordPerQuery;
        ullint readId, srcQueryOffset, srcNewQueryOffset, newReadId;
        // copy the content to the new array
        uint * newqueries = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) );
        uint * newreadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        uint * newreadIDs = ( uint * ) malloc ( roundUp * sizeof ( uint ) );

        for ( i = 0; i < numOfReadsToProcess; i++ )
        {
            readId = needProcessPair[i];
            srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            srcNewQueryOffset = i / 32 * 32 * wordPerQuery + i % 32;

            for ( j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[i] = readLengths[readId];
            newreadIDs[i] = readIDs[readId];
        }

        // append the content of the old array to the new array
        newReadId = numOfReadsToProcess;

        for ( readId = 0; readId < numOfTotalReads; readId++ )
        {
            srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            srcNewQueryOffset = newReadId / 32 * 32 * wordPerQuery + newReadId % 32;

            for ( j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[newReadId] = readLengths[readId];
            newreadIDs[newReadId] = readIDs[readId];
            newReadId++;
        }

        // update the read IDs in the readInputForDPall
        for ( i = 0; i < readInputForDPall->numArrays; i++ )
        {
            ReadInputForDP * currArray = readInputForDPall->inputArrays[i];

            // sa list
            for ( j = 0; j < currArray->saRangeTotalNum; j++ )
            {
                ( currArray->sa_list[j] ).readID += numOfReadsToProcess;
            }

            // occ list
            for ( j = 0; j < currArray->occTotalNum; j++ )
            {
                ( currArray->occ_list[j] ).readID += numOfReadsToProcess;
            }
        }

        // update the read IDs in the readInputForNewDPall
        for ( i = 0; i < readInputForNewDPall->numArrays; i++ )
        {
            ReadInputForDP * currArray = readInputForNewDPall->inputArrays[i];

            // sa list
            for ( j = 0; j < currArray->saRangeTotalNum; j++ )
            {
                ( currArray->sa_list[j] ).readID += numOfReadsToProcess;
            }

            // occ list
            for ( j = 0; j < currArray->occTotalNum; j++ )
            {
                ( currArray->occ_list[j] ).readID += numOfReadsToProcess;
            }
        }

        // update the read IDs in the bothUnalignedPairsArrays
        for ( i = 0; i < bothUnalignedPairsArrays->arrayNum; i++ )
        {
            BothUnalignedPairs * currArray = bothUnalignedPairsArrays->array[i];

            for ( j = 0; j < currArray->totalNum; j++ )
            {
                ( currArray->readIDs ) [j] += numOfReadsToProcess;
            }
        }

        // free the memory of the original array
        free ( queries );
        free ( readLengths );
        free ( readIDs );
        ( *orgqueries ) = newqueries;
        ( *orgreadLengths ) = newreadLengths;
        ( *orgreadIDs ) = newreadIDs;
    }
}

#define lastAlignedBoundary(offset, alignBoundary)              ( (offset) & (- (alignBoundary)) )
void retrieve2BWTQueryFormat ( unsigned int * packedPatterns,
                               unsigned int queryIdx, unsigned int wordPerQuery, unsigned char * unpackedPattern )
{
    unsigned int i;
    int j;
    int k = 0;
    unsigned int warpSkipSize = lastAlignedBoundary ( queryIdx, BGS_DEVICE_WARP_SIZE ) * wordPerQuery;
    packedPatterns += ( warpSkipSize + queryIdx % BGS_DEVICE_WARP_SIZE );

    for ( i = 0; i < wordPerQuery; i++ )
    {
        unsigned int currentWord = * ( packedPatterns );

        for ( j = 0; j < CHAR_PER_WORD; j++ )
        {
            char bit = currentWord & CHAR_MASK;
            currentWord >>= BIT_PER_CHAR;
            unpackedPattern[k++] = bit;
        }

        packedPatterns += BGS_DEVICE_WARP_SIZE;
    }
}

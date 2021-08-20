/*
 *
 *    OutputDPResult.c
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

#include "OutputDPResult.h"
#include "omp.h"
#include "assert.h"

void getReadInfo ( SRAQueryInput * qInput, unsigned int readID, unsigned int displayID,
                   unsigned char * upkdQueries, char ** queryNames, char ** queryComments, unsigned int * upkdReadLengths,
                   unsigned int maxReadLength, int maxReadNameLength )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    qInfo->ReadName = queryNames[readID];
    qInfo->ReadCode = upkdQueries + ( unsigned long long ) readID * maxReadLength;
    qInfo->ReadLength = upkdReadLengths[readID];
    qInfo->ReadId = displayID;
}

void getReadInfoPackedSeq ( SRAQueryInput * qInput, unsigned int readID, unsigned int displayID,
                            unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char ** queryNames,
char ** queryComments,                             unsigned int word_per_query, unsigned char * unpackedQueries, char * upkdQualities,
                            unsigned int maxReadLength, int maxReadNameLength )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    qInfo->ReadName = queryNames[upkdReadIDs[readID] - 1];
    retrieve2BWTQueryFormat ( queries, readID, word_per_query, unpackedQueries );
    qInfo->ReadCode = unpackedQueries;
    qInfo->ReadLength = readLengths[readID];
    qInfo->ReadId = displayID;
    qInfo->ReadQuality = upkdQualities + ( unsigned long long ) ( upkdReadIDs[readID] - 1 ) * maxReadLength;
}

void getAlgnInfo ( PEPairs & pePair, AlgnmtDPResult & dpResult, int peStrandLeftLeg, int peStrandRightLeg )
{
    pePair.algnmt_1 = dpResult.algnmt_1;
    pePair.strand_1 = dpResult.strand_1;
    pePair.mismatch_1 = dpResult.score_1;
    pePair.algnmt_2 = dpResult.algnmt_2;
    pePair.strand_2 = dpResult.strand_2;
    pePair.mismatch_2 = dpResult.score_2;
}

// output best alignment result pointer
// pointer is pointing to data of algnResult, no free to the output
DeepDPAlignResult * outputDeepDPResult2 ( DeepDPAlignResult * algnResult, unsigned int algnNum,
                           unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                           unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
                           SeedPool *seedPool, OCC * tmpocc, DynamicUint8Array * tmpCharArray, char ** twoSamStrings, int threadId )
{
    if ( algnNum == 0 )
    {
        return NULL;
    }
    BWT * bwt = index->sraIndex->bwt;
    HSPAux * hspaux = index->sraIndex->hspaux;
    
    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy(&(aIndex),index->sraIndex,sizeof(SRAIndex));
    
    SRASetting qSetting;
    qSetting.OutFilePtr = NULL;
    if ( tmpocc != NULL )
        { qSetting.occ = tmpocc; }
    else
        { qSetting.occ = OCCConstruct(); } 
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int upkdReadID;
    unsigned int displayID;
    DeepDPAlignResult * bestAlgn = NULL;
    int maxScore;
    unsigned char upkdPosQuery[MAX_READ_LENGTH];
    unsigned char upkdPosQuery2[MAX_READ_LENGTH];
    char * qualities1, *qualities2;
    unsigned int readLen1, readLen2;
    char * queryName1, *queryName2;
    char * queryComment1, * queryComment2;
    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    unsigned int startIndex = 0;
    DynamicUint8Array * charArray = NULL;
    
    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { 
        if ( tmpCharArray == NULL)
            { charArray = DynamicUint8ArrayConstruct(); }
        else 
            { charArray = tmpCharArray; }
    }
    double startTime = setStartTime();
    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        int isFirstUnaligned = ( algnResult[i].algnmt_1 == NOT_ALIGNED );
        int isSecondUnaligned = ( algnResult[i].algnmt_2 == NOT_ALIGNED );
        int properlyPaired = ( ( !isFirstUnaligned ) && ( !isSecondUnaligned ) );
        if ( properlyPaired && seedPool != NULL)
        {
            freeSeedAlignmentsOfReadWithReadId( seedPool, algnResult[i].readID );
            freeSeedAlignmentsOfReadWithReadId( seedPool, algnResult[i].readID+1 );
        }
        if ( algnResult[i].readID != preReadID )
        {
            if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
            {
                if ( preReadID != 0xFFFFFFFF )
                {
                    // ==========================================================
                    // for those not properly paired, needs to perform deep-dp for unaligned reads,
                    // ==========================================================
                    if ( bestAlgn->algnmt_1 == NOT_ALIGNED && bestAlgn->algnmt_2 == NOT_ALIGNED && ( ( AllHits * ) hspaux->allHits )->readNum < MAX_UNALIGN_READS_NUM_FOR_DEEP_DP )
                    {
                        int readID = preReadID;
                        addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID );
                        addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID + 1 );
                    }
                    else
                    {
                        pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                                 startIndex, i - startIndex,
                                                 upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                                 readLen1, readLen2, queryName1, queryName2, queryComment1, queryComment2,
                                                 charArray, twoSamStrings, threadId );
                    }

                    // ==========================================================
                }

                // reset the bestAlgn
                bestAlgn = & ( algnResult[i] );

                if ( isFirstUnaligned && ( !isSecondUnaligned ) )
                {
                    maxScore = bestAlgn->score_2;
                }
                else if ( isSecondUnaligned && ( !isFirstUnaligned ) )
                {
                    maxScore = bestAlgn->score_1;
                }
                else if ( properlyPaired )
                {
                    maxScore = bestAlgn->score_1 + bestAlgn->score_2;
                }
                else
                {
                    maxScore = -127;
                }
            }

            // get the information of both reads
            readLen1 = readLengths[algnResult[i].readID]; // length of the first read
            readLen2 = readLengths[algnResult[i].readID + 1]; // length of the second read
            // name of the first read
            queryName1 = queryNames[upkdReadIDs[algnResult[i].readID] - 1];
            queryComment1 = queryComments[upkdReadIDs[algnResult[i].readID] - 1];
            // name of the second read
            queryName2 = queryNames[upkdReadIDs[algnResult[i].readID + 1] - 1];
            queryComment2 = queryComments[upkdReadIDs[algnResult[i].readID + 1] - 1];
            // sequence the first read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID, word_per_query, upkdPosQuery );
            // sequence the second read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID + 1, word_per_query, upkdPosQuery2 );
            // base qualities of the first read
            qualities1 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID] - 1 ) * maxReadLength;
            // base qualities of the second read
            qualities2 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID + 1] - 1 ) * maxReadLength;
            // set qInput for the second read
            qInfo.ReadName = queryName2;
            qInfo.ReadCode = upkdPosQuery2;
            qInfo.ReadLength = readLen2;
            qInfo.ReadQuality = qualities2;
            upkdReadID = upkdReadIDs[algnResult[i].readID + 1];
            displayID = ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) ? ( upkdReadID + accumReadNum ) / 2 : upkdReadID + accumReadNum;
            qInfo.ReadId = displayID;
            startIndex = i;

            preReadID = algnResult[i].readID;
        }
        else if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            // SAM format
            // update the bestAlgn
            if ( isFirstUnaligned && !isSecondUnaligned )
            {
                if ( algnResult[i].score_2 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_2;
                }
            }
            else if ( isSecondUnaligned && !isFirstUnaligned )
            {
                if ( algnResult[i].score_1 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_1;
                }
            }
            else if ( properlyPaired )
            {
                if ( algnResult[i].score_1 + algnResult[i].score_2 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_1 + algnResult[i].score_2;
                }
            }
        }
    }

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && preReadID != 0xFFFFFFFF )
    {
        // ==========================================================
        // for those not properly paired, needs to perform deep-dp for unaligned reads,
        // ==========================================================
        if ( bestAlgn->algnmt_1 == NOT_ALIGNED && bestAlgn->algnmt_2 == NOT_ALIGNED && ( ( AllHits * ) hspaux->allHits )->readNum < MAX_UNALIGN_READS_NUM_FOR_DEEP_DP )
        {
            int readID = preReadID;
            addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID );
            addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID + 1 );
        }
        else
        {
            pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                     startIndex, algnNum - startIndex,
                                     upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                     readLen1, readLen2, queryName1, queryName2, queryComment1, queryComment2,
                                     charArray, twoSamStrings, threadId );            
        }

        // ==========================================================
    }
    
    if ( tmpocc == NULL ) 
        { OCCFree ( qSetting.occ ); }
    
    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && tmpCharArray == NULL )
    { DynamicUint8ArrayFree ( charArray ); }
    return bestAlgn;
}

void outputDPSingleResult2 ( SingleAlgnmtResult * algnResult, unsigned int algnNum,
                             unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                             unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, Soap3Index * index )
{
    if ( algnNum == 0 ) { return; }
    
    HSPAux * hspaux = index->sraIndex->hspaux;
    assert ( hspaux->readType == PAIR_END_READ );

    // preforming deep dp for the unaligned pair-end reads
    // save all the alignment results to array
    inputAlgnmtsToArray ( ( AllHits * ) hspaux->allHits, algnResult, algnNum );
}


// output the single-end alignment results for the pair-end reads
// only for SAM format
void outputSingleResultForPairEnds ( AllHits * allHits,
                                     unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                                     DPInfoForReads * dpInfoForReads,
                                     unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                                     samfile_t * samOutputUnpairFilePtr, Soap3Index * index )
{
    if ( allHits->readNum == 0 )
    {
        return;
    }

    HSPAux * hspaux = index->sraIndex->hspaux;
    
    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy(&(aIndex),index->sraIndex,sizeof(SRAIndex));
    
    SRASetting qSetting;
    qSetting.occ = OCCConstruct();
    qSetting.SAMOutFilePtr = samOutputUnpairFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;

    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    vector<DynamicUint8Array*> charArrays(omp_get_max_threads());

    for (size_t i = 0; i < charArrays.size(); ++i) {
        charArrays[i] = DynamicUint8ArrayConstruct();
    }

    if (hspaux->outputBAM) {
        omp_set_num_threads(1);
    }

    #pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < allHits->readNum - 1; i += 2 )
    {   
        unsigned int firstReadID;
        unsigned int secondReadID;
        unsigned int upkdReadID;
        unsigned int displayID;
        unsigned char upkdPosQuery[MAX_READ_LENGTH];
        unsigned char upkdPosQuery2[MAX_READ_LENGTH];
        char * qualities1, *qualities2;
        unsigned int readLen1, readLen2;
        char * queryName1, *queryName2;
        char * queryComment1, *queryComment2;

        Algnmt * algn_list1;
        Algnmt * algn_list2;
        Algnmt * bestAlgns[2] = {0};
        DynamicUint8Array * charArray = charArrays[omp_get_thread_num()];


        firstReadID = allHits->readPtrArray[i].readID;
        secondReadID = allHits->readPtrArray[i + 1].readID;
        if ( ( firstReadID + 1 ) != secondReadID )
        {
            // Error! It should not happen
            fprintf ( stderr, "%u %u\n",firstReadID, secondReadID);
            fprintf ( stderr, "[Inside the function: outputSingleResultForPairEnds] Error appears. Reads are not paired-end.\n" );
            exit(-1);
        }

        algn_list1 = & ( allHits->hitArray[ ( allHits->readPtrArray[i].startIndex )] );
        algn_list2 = & ( allHits->hitArray[ ( allHits->readPtrArray[i + 1].startIndex )] );
        // get the information of both reads
        readLen1 = readLengths[firstReadID]; // length of the first read
        readLen2 = readLengths[secondReadID]; // length of the second read
        // name of the first read
        queryName1 = queryNames[upkdReadIDs[firstReadID] - 1];
        queryComment1 = queryComments[upkdReadIDs[firstReadID] - 1];
        // name of the second read
        queryName2 = queryNames[upkdReadIDs[secondReadID] - 1];
        queryComment2 = queryComments[upkdReadIDs[secondReadID] - 1];
        // sequence the first read
        retrieve2BWTQueryFormat ( queries, firstReadID, word_per_query, upkdPosQuery );
        // sequence the second read
        retrieve2BWTQueryFormat ( queries, secondReadID, word_per_query, upkdPosQuery2 );
        // base qualities of the first read
        qualities1 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[firstReadID] - 1 ) * maxReadLength;
        // base qualities of the second read
        qualities2 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[secondReadID] - 1 ) * maxReadLength;
        // set qInput for the second read
        qInfo.ReadName = queryName2;
        qInfo.ReadCode = upkdPosQuery2;
        qInfo.ReadLength = readLen2;
        qInfo.ReadQuality = qualities2;
        upkdReadID = upkdReadIDs[secondReadID];
        displayID = upkdReadID + accumReadNum;
        qInfo.ReadId = displayID;
        unproperlypairDPOutputSAMAPI ( &qInput, firstReadID, algn_list1, secondReadID, algn_list2,
                                       allHits->readPtrArray[i].numAlgnmt, allHits->readPtrArray[i + 1].numAlgnmt,
                                       upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                       readLen1, readLen2, queryName1, queryName2, queryComment1, queryComment2,
                                       charArray, bestAlgns );
        if ( bestAlgns[0] != NULL )
        {
            dpInfoForReads->startPositions[firstReadID] = bestAlgns[0]->startPosition;
            dpInfoForReads->strand_dpLengths[firstReadID] = bestAlgns[0]->refDpLength | ( ( bestAlgns[0]->strand - 1 ) << 15 );
            dpInfoForReads->peLeftAnchors[firstReadID] = bestAlgns[0]->peLeftAnchor;
            dpInfoForReads->peRightAnchors[firstReadID] = bestAlgns[0]->peRightAnchor;
        }
        if ( bestAlgns[1] != NULL )
        {
            dpInfoForReads->startPositions[secondReadID] = bestAlgns[1]->startPosition;
            dpInfoForReads->strand_dpLengths[secondReadID] = bestAlgns[1]->refDpLength | ( ( bestAlgns[1]->strand - 1 ) << 15 );
            dpInfoForReads->peLeftAnchors[secondReadID] = bestAlgns[1]->peLeftAnchor;
            dpInfoForReads->peRightAnchors[secondReadID] = bestAlgns[1]->peRightAnchor;
        }
    }

    for (size_t i = 0; i < charArrays.size(); ++i) {
        DynamicUint8ArrayFree ( charArrays[i] );
    }
    OCCFree ( qSetting.occ );
}

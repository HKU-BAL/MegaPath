/*
 *
 *    alignment.cu
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
#include <iostream>
#include <assert.h>
#include "alignment.h"
#include "DV-DPfunctions.h"

// Perform SOAP3-DP Paired-End Alignment
void soap3_dp_pair_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                           ullint maxBatchSize, uint numQueries, uint accumReadNum,
                           Soap3Index * index,
                           IniParams ini_params, InputOptions input_options,
                           uint maxReadLength, uint detected_read_length, uint detected_read_length2,
                           char * upkdQualities,
                           uint * readIDs, char ** queryNames, char ** queryComments,
                           DPInfoForReads * dpInfoForReads, samfile_t ** currSamOutputFilePtr, samfile_t * samOutputDPFilePtr,
                           samfile_t * samOutputUnpairFilePtr,
                           unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays, SeedPool * seedPool,
                           double startTime, double & lastEventTime, double & totalAlignmentTime)
{
    double alignmentTime, copyTime;
    unsigned int numDPAlignedPair = 0;
    unsigned int numDPAlignment = 0;
    DPParameters dpParameters;
    int orig_align_type;
    HSPAux * hspaux = index->sraIndex->hspaux;
    SingleDP_Space::AlgnmtResultStream* singleDPResultStream = NULL;

    if ( ini_params.Ini_skipSOAP3Alignment == 1 )
    {
        // Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
        addAllFirstReadIDToBothUnalignedPairs ( bothUnalignedPairsArrays->array[0], numQueries );
        // For DP module, if SAM format and all-best are both selected, then
        // output format is needed to set to all-valid.
        orig_align_type = input_options.alignmentType;

        if ( input_options.alignmentType == OUTPUT_ALL_BEST)
        {
            input_options.alignmentType = OUTPUT_ALL_VALID;
        }

        // Parameters for DP
        getParameterForAllDP ( dpParameters, input_options, ini_params );
    } else {
        assert(false);
    }

#ifdef PERFORM_DEEP_DP
    ////////////////////////////////////////////////////////////
    // PERFORM DP FOR BOTH ENDS UNALIGNED                     //
    ////////////////////////////////////////////////////////////
    numDPAlignedPair = 0;
    numDPAlignment = 0;

    int numberOfRoundOfDeepDp = 0;
    SeedingProperties * seedingProperties=NULL;
    if ( maxReadLength > 120 )
    {
        numberOfRoundOfDeepDp = ini_params.Ini_numberOfRoundOfDeepDpForLongReads;
        seedingProperties = ini_params.seedingPropertiesForLongReads;
    }
    else
    {
        numberOfRoundOfDeepDp = ini_params.Ini_numberOfRoundOfDeepDpForShortReads;
        seedingProperties = ini_params.seedingPropertiesForShortReads;
    }

    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
    {
        for (int roundNumber = 1; roundNumber <= numberOfRoundOfDeepDp; ++roundNumber )
        {
            hspaux->dpStageId = roundNumber;
            // DP Parameters for Deep DP

            double ini_params_cutoff_threshold_bak = ini_params.Ini_DPScoreThreshold;
            if (seedingProperties[roundNumber-1].dpScoreThreshold >= 0) {
                ini_params.Ini_DPScoreThreshold = seedingProperties[roundNumber-1].dpScoreThreshold;
            }
            getParameterForDeepDP ( dpParameters, input_options, ini_params, detected_read_length, detected_read_length2, maxReadLength, roundNumber );
            // printDPParameters ( dpParameters );
            hspaux->ProceedDPForTooManyHits = seedingProperties[roundNumber-1].proceedDPForTooManyHits;
            ini_params.Ini_DPScoreThreshold = ini_params_cutoff_threshold_bak;

            unsigned int totalReadsProceedToDP = 0;

            for ( int threadId = 0; threadId <= ini_params.Ini_NumOfCpuThreads; ++threadId )
            {
                totalReadsProceedToDP += bothUnalignedPairsArrays->array[threadId]->totalNum;
            }

            if ( totalReadsProceedToDP > 0 )
            {
                fprintf (stderr, "[Main] %u pairs of reads are proceeded to deep DP Round %d.\n", totalReadsProceedToDP, roundNumber );
                DPForUnalignPairs2 ( bothUnalignedPairsArrays, seedPool, input_options.insert_high, input_options.insert_low,
                                     queries, readLengths, readIDs, queryNames, queryComments, upkdQualities,
                                     maxReadLength, ini_params.Ini_PEStrandLeftLeg, ini_params.Ini_PEStrandRightLeg,
                                     dpInfoForReads,
                                     index,
                                     input_options.alignmentType, &dpParameters,
                                     numDPAlignedPair, numDPAlignment,
                                     accumReadNum, SRA_OUTPUT_FORMAT_SAM_API, samOutputDPFilePtr, currSamOutputFilePtr);
                // fprintf("stderr,Finished DP for both-end unaligned reads\n");
                // fprintf("stderr,[Main] Number of pairs aligned by DP: %u (number of alignments: %u)\n", numDPAlignedPair, numDPAlignment);
                fprintf (stderr, "[Main] Number of pairs aligned by DP: %u\n", numDPAlignedPair );
                fprintf (stderr, "[Main] Number of alignments aligned by DP: %u\n", numDPAlignment );
                alignmentTime = getElapsedTime ( startTime );
                fprintf (stderr, "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
                totalAlignmentTime += alignmentTime - lastEventTime;
                lastEventTime = alignmentTime;
                numOfAlignedRead +=  numDPAlignedPair * 2;
                numOfAnswer += numDPAlignment;
                // fprintf("stderr,[Main] Total Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
                fprintf (stderr, "[Main] Total Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
                fprintf (stderr, "\n" );
            }
        }
    }
#endif

    ////////////////////////////////////////////////////////////
    // TO PERFORM SINGLE DEEP-DP FOR UNALIGNED READS          //
    ////////////////////////////////////////////////////////////
    unsigned int numSingleDPAligned = 0;
    unsigned int numSingleDPAlignment = 0;
    QueryIDStream * unalignedReadsFromSingleEndDP = new QueryIDStream;
    UnalignedSinglesArrays * unalignedSingleEndArrays = NULL;
    if ( bothUnalignedPairsArrays->array[0]->totalNum != 0 )
    {
        unalignedSingleEndArrays = (UnalignedSinglesArrays *) malloc(sizeof (UnalignedSinglesArrays));
        unalignedSingleEndArrays->arrayNum = 1;
        unalignedSingleEndArrays->array = (UnalignedSingles **) malloc( sizeof (UnalignedSingles*) );
        unalignedSingleEndArrays->array[0] = (UnalignedSingles *) malloc( sizeof (UnalignedSingles) );
        int new_size = bothUnalignedPairsArrays->array[0]->totalNum * 2;
        unalignedSingleEndArrays->array[0]->readIDs = (uint *) malloc ( new_size * sizeof (uint) );
        unalignedSingleEndArrays->array[0]->size = new_size;
        unalignedSingleEndArrays->array[0]->totalNum= new_size;
        for (int i=0;i<new_size;i+=2)
        {
            unalignedSingleEndArrays->array[0]->readIDs[i] = bothUnalignedPairsArrays->array[0]->readIDs[i/2];
            unalignedSingleEndArrays->array[0]->readIDs[i+1] = bothUnalignedPairsArrays->array[0]->readIDs[i/2]+1;
        }
    }

    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ && numberOfRoundOfDeepDp >= 1 && unalignedSingleEndArrays != NULL )
    {
        // hspaux->dpStageId++;
        // reset the readType to single read
        input_options.readType = SINGLE_READ;
        // DP Parameters for single-end alignment
        getParameterForSingleDP ( dpParameters, input_options, ini_params, detected_read_length );
        hspaux->singleDPcutoffRatio = DP_SCORE_THRESHOLD_RATIO;
        hspaux->singleDPcutoffLB = DP_SCORE_THRESHOLD_LOWER_BOUND;
        // printDPParameters(dpParameters);
        unsigned int totalReadsProceedToDP = 0;

        totalReadsProceedToDP += unalignedSingleEndArrays->array[0]->totalNum;
        fprintf (stderr, "[Main] %u unaligned reads are proceeded to DP.\n", totalReadsProceedToDP );

        if ( totalReadsProceedToDP > 0 )
        {
            DPForUnalignSingle2 ( unalignedSingleEndArrays, seedPool,
                                  queries, readLengths, readIDs, queryNames, queryComments, upkdQualities,
                                  dpInfoForReads,
                                  maxReadLength,
                                  index,
                                  input_options.alignmentType, &dpParameters,
                                  numSingleDPAligned, numSingleDPAlignment,
                                  accumReadNum, SRA_OUTPUT_FORMAT_SAM_API, 
                                  NULL, singleDPResultStream,
                                  unalignedReadsFromSingleEndDP );

            int sum=0,cnt=0,pre_readID=-1, ac=0;
            for (int i=0; i < singleDPResultStream->dpSResult.size(); ++i)
            {
                SingleDP_Space::SingleDPResultBatch *batch = singleDPResultStream->dpSResult[i];
                for (int j=0;j<batch->size();++j )
                {
                    ++cnt;
                    if ( (*batch)[j].readID != pre_readID ){
                        ac += cnt;
                        cnt = 0;
                        pre_readID = (*batch)[j].readID;
                    }
                }
                sum += batch->size();
            }
        }
        // fprintf("stderr,Finished DP for single unaligned reads\n");
        // fprintf("stderr,[Main] Number of reads aligned by DP: %u (number of alignments: %u)\n", numDPAlignedSingle, numDPAlignment);
        fprintf (stderr, "[Main] Number of reads aligned by single-end DP: %u\n", numSingleDPAligned );
        fprintf (stderr, "[Main] Number of alignments aligned by DP: %u\n", numSingleDPAlignment );
        alignmentTime = getElapsedTime ( startTime );
        fprintf (stderr, "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
        totalAlignmentTime += alignmentTime - lastEventTime;
        lastEventTime = alignmentTime;

        fprintf (stderr, "\n" );
        input_options.readType = PAIR_END_READ;
    }

    ////////////////////////////////////////////////////////////
    // PERFORM SEMI-GLOBAL DP IF NECESSARY                    //
    ////////////////////////////////////////////////////////////
#ifndef SKIP_DEFAULT_DP
    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ && unalignedSingleEndArrays != NULL && ini_params.Ini_skipDefaultDP == 0 )
    {
        hspaux->dpStageId++;
        // DP Parameters for half-aligned reads
        getParameterForDefaultDP ( dpParameters, input_options, ini_params, detected_read_length, detected_read_length2 );
        // printDPParameters ( dpParameters );

        unsigned int totalReadsProceedToDP = 0;
        uint preReadId = 0xffffffff;
        for ( int arrayIndex = 0; arrayIndex < singleDPResultStream->dpSResult.size(); ++arrayIndex )
        {
            for (int i=0;i<singleDPResultStream->dpSResult[arrayIndex]->size(); ++i)
            {
                if ( preReadId != (*(singleDPResultStream->dpSResult[arrayIndex]))[i].readID >> 1)
                {
                    preReadId = (*(singleDPResultStream->dpSResult[arrayIndex]))[i].readID >> 1;
                    ++totalReadsProceedToDP;
                }
            }
        }

        if ( totalReadsProceedToDP &&  singleDPResultStream != NULL && singleDPResultStream->numOut != 0 )
        {
            fprintf (stderr, "[Main] %u half-aligned pairs of reads are proceeded to DP.\n", totalReadsProceedToDP );
            semiGlobalDPForSingleEndDpAlignment ( singleDPResultStream, input_options.insert_high, input_options.insert_low,
                            queries, readLengths, readIDs, queryNames, queryComments, upkdQualities,
                            dpInfoForReads,
                            maxReadLength, ini_params.Ini_PEStrandLeftLeg, ini_params.Ini_PEStrandRightLeg,
                            index,
                            input_options.alignmentType, &dpParameters,
                            numDPAlignedPair, numDPAlignment,
                            accumReadNum, SRA_OUTPUT_FORMAT_SAM_API,  samOutputDPFilePtr, currSamOutputFilePtr,
                            bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads] );

            // fprintf("stderr,Finished semi-global DP\n");
            // fprintf("stderr,[Main] Number of pairs aligned by DP: %u (number of alignments: %u)\n", numDPAlignedPair, numDPAlignment);
            fprintf (stderr, "[Main] Number of pairs aligned by DP: %u\n", numDPAlignedPair );
            fprintf (stderr, "[Main] Number of alignments aligned by DP: %u\n", numDPAlignment );
            // fprintf("stderr, Number of reads passed to deep-dp: %u\n", bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads]->totalNum);
            alignmentTime = getElapsedTime ( startTime );
            fprintf (stderr, "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
            totalAlignmentTime += alignmentTime - lastEventTime;
            lastEventTime = alignmentTime;
            numOfAlignedRead +=  numDPAlignedPair * 2;
            numOfAnswer += numDPAlignment;
            // fprintf("stderr,[Main] Total Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
            fprintf (stderr, "[Main] Total Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
            fprintf (stderr, "\n" );
        }
    }
    else if ( singleDPResultStream != NULL )
    {
        BothUnalignedPairs * bothUnalignedPairs = bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads];
        uint lastReadID=0xffffffff;
        for (int i=0; i < singleDPResultStream->dpSResult.size(); ++i)
        {
            SingleDP_Space::SingleDPResultBatch *batch = singleDPResultStream->dpSResult[i];
            for (int j=0;j<batch->size();++j )
            {   
                if ( lastReadID != ( (*batch)[j].readID & 0xfffffffe) )
                {
                    addReadIDToBothUnalignedPairs( bothUnalignedPairs, ( (*batch)[j].readID & 0xfffffffe) );
                    lastReadID = ( (*batch)[j].readID & 0xfffffffe);
                }
            }
        }
    }
    else
    {
        singleDPResultStream = new SingleDP_Space::AlgnmtResultStream();
    }
#endif

    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
    {
        BothUnalignedPairs * bothUnalignedPairs = bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads];
        UnalignedSingles * unalignedSingleReads = ( UnalignedSingles * ) malloc ( sizeof ( UnalignedSingles ) );
        
        filterOutUnpairedSingleReads ( unalignedSingleReads, bothUnalignedPairs, unalignedReadsFromSingleEndDP->data );
        
        //Output file for SAM output is handled by SAM API
        //=> OutFilePtr is not used for SAM API output format.
        hspaux->dpStageId++;
        uint dpSAlignedRead = 0, dpSAlignment = 0; 
        // this is to collect all hits in hspaux->allHits
        outputUnpairedSingleReads ( singleDPResultStream, unalignedSingleReads,
                                     queries, readLengths, readIDs, 
                                     queryNames, queryComments, upkdQualities, maxReadLength,
                                     accumReadNum, SRA_OUTPUT_FORMAT_SAM_API,  samOutputDPFilePtr,
                                     index,
                                     input_options.alignmentType,
                                     dpSAlignedRead, dpSAlignment );
        
        // sort the alignment results
        sortReadPtrs ( ( AllHits * ) hspaux->allHits );
        // output the alignment results
        outputSingleResultForPairEnds ( ( AllHits * ) hspaux->allHits,
                                        queries, readLengths, readIDs, queryNames, queryComments, upkdQualities,
                                        dpInfoForReads,
                                        maxReadLength, accumReadNum, SRA_OUTPUT_FORMAT_SAM_API, 
                                        samOutputUnpairFilePtr,
                                        index );
        
        // testing
        int preReadID = 0;
        for ( int i=0;i< singleDPResultStream->dpSResult.size(); ++i )
        {
            SingleDP_Space::SingleDPResultBatch * singleDPResultBatch = singleDPResultStream->dpSResult[i];
            for (int j=0; j < singleDPResultBatch->size(); ++j)
            {
                if ( (*singleDPResultBatch)[j].readID < preReadID )
                    fprintf( stderr, "Single-end DP result not in order %u %u\n", (*singleDPResultBatch)[j].readID, preReadID);
                preReadID = (*singleDPResultBatch)[j].readID ;
            }
        }
        
        fprintf( stderr,"[Main] Number of aligned reads: %d\n", dpSAlignedRead);
        fprintf( stderr,"[Main] Number of aligned DP: %d\n", dpSAlignment);
        alignmentTime = getElapsedTime ( startTime );
        fprintf (stderr, "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
        totalAlignmentTime += alignmentTime - lastEventTime;
        lastEventTime = alignmentTime;
        freeBothUnalignedPairs ( unalignedSingleReads);
        freeBothUnalignedPairsArrays ( unalignedSingleEndArrays );
        input_options.readType = PAIR_END_READ;
    }
    delete unalignedReadsFromSingleEndDP;
    delete singleDPResultStream;
    input_options.alignmentType = orig_align_type;
}

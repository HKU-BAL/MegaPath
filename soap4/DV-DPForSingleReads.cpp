/*
 *
 *    DV-DPForSingleReads.cu
 *    Soap3(gpu)
 *
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

#include "DV-DPForSingleReads.h"
#include "OutputDPResult.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <semaphore.h>

#include <functional>
#include <vector>
#include <parallel/algorithm>
#include "kxsort.h"
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace SingleDP_Space;
class SingleDPWrapper
{
        UnalignedSinglesArrays * unalignedReads;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char *upkdQualities;
        char ** queryNames;
        char ** queryComments;
        DPInfoForReads * dpInfoForReads;
        uint maxReadLength;
        Soap3Index * index;
        int alignmentType;
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;

        DPParameters * dpParameters;
        SOAP3Wrapper<void> * soap3Wrapper;
        SeedPool * seedPool;
        
    public:
        AlgnmtResultStream * resultStream;
        uint numDPAlignedRead, numDPAlignment;

        SingleDPWrapper ( UnalignedSinglesArrays * unalignedReads, SeedPool* seedPool,
                          uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                          char ** queryNames, char ** queryComments, char * upkdQualities, uint maxReadLength,
                          DPInfoForReads * dpInfoForReads,
                          Soap3Index * index,
                          int alignmentType, DPParameters * dpParameters,
                          uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr )
        {
            MC_MemberCopy3 ( this->, , unalignedReads, dpParameters, seedPool );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy4 ( this->, , queryNames, queryComments, upkdQualities, maxReadLength );
            MC_MemberCopy ( this->, , dpInfoForReads );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy2 ( this->, , outputFormat, samOutputDPFilePtr );
            soap3Wrapper =
                new SOAP3Wrapper<void> (index);
        }
        ~SingleDPWrapper()
        {
            delete soap3Wrapper;
        }

        void alignment ( CandidateStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            SingleEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                queryNames, queryComments, origReadIDs, upkdQualities,
                dpInfoForReads,
                index,
                alignmentType,
                accumReadNum, outputFormat, samOutputDPFilePtr,
                /* output */
                unalignedIDStream, alignedRead, alignment, resultStream );
            numDPAlignedRead += alignedRead;
            numDPAlignment += alignment;
        }

        struct RadixTraitsCandidateInfo {
            bool operator()(const CandidateInfo &x, const CandidateInfo &y) {
                if (x.readID < y.readID) return true;
                else if (x.readID > y.readID) return false;
                else if (x.strand < y.strand) return true;
                else if (x.strand > y.strand) return false;
                else if (x.pos < y.pos) return true;
                else if (x.pos > y.pos) return false;
                else return x.seedAlignmentLength < y.seedAlignmentLength;
            }
        };

        void transferSeed( SeedPool * seedPool, QueryIDStream * inputStream, CandidateStream * canStream, QueryIDStream * unseededIDStream)
        {
            QueryIDStream * tmpUnseededIDStream = new QueryIDStream;
            #define MC_AppendPos(posIter, id, seedStrand, ePos, alignmentLength) { \
                posIter->readID = id; \
                posIter->pos = ePos; \
                posIter->strand = seedStrand; \
                posIter->seedAlignmentLength = alignmentLength; \
                posIter->peLeftAnchor = 0; \
                posIter->peRightAnchor = 0; \
                ++posIter; \
            }
            int readId;
            SeedAlignmentsOfRead * seedAlignments;
            CandidateInfo *candidate;
            int sum=0,numberCnt=0;
            for (int i=0;i<inputStream->data->size();++i)
            {
                readId = (*(inputStream->data))[i];
                for (int threadId=0; threadId<seedPool->numberOfThread;++threadId)
                {
                    if ( seedPool->seedAlignmentsOfReadArrays[threadId][readId] != NULL)
                    {
                        seedAlignments = seedPool->seedAlignmentsOfReadArrays[threadId][readId];
                        sum += seedAlignments->numberOfAlignments;
                    }
                }
            }
            fprintf(stderr, "[Single-end DP] Number of Seed Alignments: %d\n", sum);
            MC_CheckMalloc ( candidate,     SeedPos,    sum + 1 );
            CandidateInfo *iter_candidateInfo = candidate;
            int cnt=0;
            fprintf(stderr, "[Single-end DP] Total Reads: %u\n", inputStream->data->size());
            for (int i=0;i<inputStream->data->size();++i)
            {
                bool unseededReads=true;
                readId = (*(inputStream->data))[i];
                for (int threadId=0; threadId<seedPool->numberOfThread;++threadId)
                {
                    if ( seedPool->seedAlignmentsOfReadArrays[threadId][readId] != NULL)
                    {
                        unseededReads = false;
                        seedAlignments = seedPool->seedAlignmentsOfReadArrays[threadId][readId];
                        cnt += seedAlignments->numberOfAlignments;
                        for ( int j=0; j<seedAlignments->occList->size; ++j )
                        {
                            OccRecord occList = seedAlignments->occList->get(j);
                            // fprintf(stderr, "P: %u\n",occList.pos);
                            MC_AppendPos( iter_candidateInfo, readId, occList.strand, occList.pos, occList.seedLength );
                            // fprintf(stderr, "[%s] Seed Length %d\n", __func__, occList.seedLength);
                            ++numberCnt;
                        }
                    }
                }
                if ( unseededReads )
                    { tmpUnseededIDStream->data->push_back( readId ); }
                else
                    { numberCnt+=0; }
                freeSeedAlignmentsOfReadWithReadId(seedPool, readId);
            }
            MC_AppendPos ( iter_candidateInfo, 0x7FFFFFFF, 2, 0xFFFFFFFF, 0xFFFFFFFF );
            omp_set_num_threads(dpParameters->numOfCPUThreads);
            __gnu_parallel::sort(candidate, candidate + sum + 1, RadixTraitsCandidateInfo());
            SingleEndSeedingEngine* engine = new SingleEndSeedingEngine();

            QueryIDStream * eliminatedReadIDStream = new QueryIDStream;
            vector<CandidateInfo> * tmp = engine->singleMerge( candidate, eliminatedReadIDStream );
            fprintf(stderr, "[Single-end DP] Total Number of unseeded reads: %u\n", eliminatedReadIDStream->data->size());
            int j;
            numberCnt = 0;
            for (int i=0;i<tmp->size(); i=j)
            {
                ++numberCnt;
                j=i+1;
                while ( j<tmp->size() && (((*tmp)[j]).readID) == (((*tmp)[i]).readID) )
                    { ++j; }
                if ( j-i ) // max hit of a seed
                {
                    int l=0;
                    for ( int k=i;k<j && k <i+200;++k) // WARNING HARDCODE
                    {
                        canStream->data.push_back((*tmp)[k]);
                        ++l;
                    }
                }
            }
            
            free ( candidate );
            delete tmp;
            sort ( eliminatedReadIDStream->data->begin(), eliminatedReadIDStream->data->end() );
            mergeTwoSortedQueryIDStream( tmpUnseededIDStream, eliminatedReadIDStream, unseededIDStream );
            delete tmpUnseededIDStream;
            delete eliminatedReadIDStream;
            delete engine;
        }

        void formCandidate( QueryIDStream * inputStream, CandidateStream * canStream, QueryIDStream * unalignedIDStream )
        {
            #define MC_AppendPosForPreparedRegion(posIter, id, seedStrand, ePos, alignmentLength, leftAnchor, rightAnchor) { \
                posIter->readID = id; \
                posIter->pos = ePos; \
                posIter->strand = seedStrand; \
                posIter->seedAlignmentLength = alignmentLength; \
                posIter->peLeftAnchor = leftAnchor; \
                posIter->peRightAnchor = rightAnchor; \
                ++posIter; \
            }
            uint readID;
            uint strand;
            uint position;
            uint length;
            uint peLeftAnchor;
            uint peRightAnchor;
            uint size = inputStream->data->size();
            CandidateInfo *candidate, *auxCadidate;
            MC_CheckMalloc ( candidate,     SeedPos,    size + 1 );
            MC_CheckMalloc ( auxCadidate,  SeedPos,    size + 1 );
            CandidateInfo *iter_candidateInfo = candidate;
            uint cnt= 0;
            for (int i=0;i<size;++i)
            {
                readID = (*(inputStream->data))[i];
                position = dpInfoForReads->startPositions[readID];
                strand = ( ( dpInfoForReads->strand_dpLengths[readID] >> 15 ) & 1 ) + 1;
                length = ( ( dpInfoForReads->strand_dpLengths[readID] ) & 0x7fff );
                peLeftAnchor = dpInfoForReads->peLeftAnchors[readID];
                peRightAnchor = dpInfoForReads->peRightAnchors[readID];
                if ( position != 0xffffffff )
                { 
                    // fprintf (stderr, "%u %u\n", position, length);
                    MC_AppendPosForPreparedRegion ( iter_candidateInfo, readID, strand, position, length, peLeftAnchor, peRightAnchor );
                    ++cnt;
                }
                else 
                    { unalignedIDStream->data->push_back( readID ); }
            }
            MC_AppendPos ( iter_candidateInfo, 0x7FFFFFFF, 2, 0xFFFFFFFF, 0xFFFFFFFF );
            // MC_RadixSort_32_16 ( candidate, seedAlignmentLength, auxCadidate, cnt+1 );
            
            for (int i=0;i<cnt; i++)
            {
                canStream->data.push_back( candidate[i] );
            }
            free ( auxCadidate );
            free ( candidate );
        }

        
        void singleDPOneRound ( QueryIDStream * inputStream,
                                QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            CandidateStream * canStream = new CandidateStream;
            // transfer seed to candidate from seedPool
            transferSeed( seedPool, inputStream, canStream, unseededIDStream);
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }


        // Assume there is only one result for each readID
        /*
        void sortResultStreamByReadID ( AlgnmtResultStream * resultStream )
        { 
            uint size = resultStream->dpSResult.size();
            SingleAlgnmtResult * tmpArray = ( SingleAlgnmtResult * ) malloc ( size * sizeof ( SingleAlgnmtResult ) );
            SingleAlgnmtResult * bufferArray = ( SingleAlgnmtResult * ) malloc ( size * sizeof ( SingleAlgnmtResult ) );
            for ( int i=0;i<size;++i )
                { tmpArray[i]. }
            
        }
        */
        
        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DPSOutputUnalignedReads (
                unalignedIDStream,
                queries, upkdReadLengths, maxReadLength,
                index, queryNames, queryComments, origReadIDs, upkdQualities,
                accumReadNum, outputFormat, samOutputDPFilePtr
            );
        }

        QueryIDStream * run()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            QueryIDStream * input = new QueryIDStream ( unalignedReads );
            QueryIDStream * unaligned_round1 = new QueryIDStream;
            QueryIDStream * unseededIDStream = new QueryIDStream;
            singleDPOneRound ( input, unseededIDStream, unaligned_round1 );
            delete input;
            // outputUnaligned ( unaligned_round1 );
            // delete unaligned_round1;
            QueryIDStream * unalignedIDStream = new QueryIDStream;
            mergeTwoSortedQueryIDStream ( unseededIDStream, unaligned_round1, unalignedIDStream );
            delete unseededIDStream;
            delete unaligned_round1;
            return unalignedIDStream ;
        }
};

void DPForUnalignSingle2 ( UnalignedSinglesArrays * unalignedReads, SeedPool* seedPool,
                           unsigned int * queries, unsigned int * upkdReadLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                           DPInfoForReads * dpInfoForReads,
                           unsigned int maxReadLength,
                           Soap3Index * index,
                           int alignmentType, DPParameters * dpParameters,
                           unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                           unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                           SingleDP_Space::AlgnmtResultStream * &resultStream,
                           QueryIDStream * &unalignedSingleReads )
{
    using namespace SingleDP_Space;
    SingleDPWrapper
    singleDPWrapper (
        unalignedReads, seedPool,
        queries, upkdReadLengths, origReadIDs,
        queryNames, queryComments, upkdQualities, maxReadLength,
        dpInfoForReads,
        index, 
        alignmentType, dpParameters,
        accumReadNum, outputFormat, samOutputDPFilePtr
    );
    unalignedSingleReads = singleDPWrapper.run();
    numDPAlignedRead = singleDPWrapper.numDPAlignedRead;
    numDPAlignment = singleDPWrapper.numDPAlignment;
    resultStream = singleDPWrapper.resultStream;
}

void outputUnpairedSingleReads ( SingleDP_Space::AlgnmtResultStream * singleDPResultStream, UnalignedSingles * unalignedSingleReads,
                                 uint * queries, uint * upkdReadLengths, uint * origReadIDs, 
                                 char ** queryNames, char ** queryComments, char * upkdQualities, int inputMaxReadLength,
                                 uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                                 Soap3Index * index,
                                 int alignmentType,
                                 uint         &        dpSAlignedRead,
                                 uint         &        dpSAlignment   )
{
    SingleDP_Space::DPSOutputUnpairedAlignment( singleDPResultStream, unalignedSingleReads,
                               queries, upkdReadLengths, origReadIDs, 
                               queryNames, queryComments, upkdQualities, inputMaxReadLength,
                               accumReadNum, outputFormat, samOutputDPFilePtr,
                               index,
                               alignmentType,
                               dpSAlignedRead, dpSAlignment );    
}

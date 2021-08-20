/*
 *
 *    DV-DPForBothUnalign.cu
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


#include "DV-DPForBothUnalign.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include "kxsort.h"
using namespace std;

#define DP2_SEEDING_BATCH_SIZE 128 * 1024 
// #define DP2_MARGIN(l) ((l>100) ? 30 : 25)

#define MC_Max(x,y) (x > y ? x : y)

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace DeepDP_Space;
class DeepDPWrapper
{
        BothUnalignedPairsArrays * unalignedReads;
        SeedPool *seedPool;
        int insert_high, insert_low;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char *upkdQualities;
        char ** queryNames;
        char ** queryComments;
        uint maxReadLength;
        int peStrandLeftLeg, peStrandRightLeg;
        DPInfoForReads * dpInfoForReads;
        Soap3Index * index;
        int alignmentType;
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;
        samfile_t ** samOutputDPFilePtrs;

        bool shortDnaLength;
        DPParameters * dpParameters;
        SOAP3Wrapper<void> * soap3Wrapper;

    public:
        uint numDPAlignedRead, numDPAlignment;

        DeepDPWrapper ( BothUnalignedPairsArrays * unalignedReads, SeedPool* seedPool, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char ** queryNames, char ** queryComments, char * upkdQualities, uint maxReadLength,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        DPInfoForReads * dpInfoForReads,
                        Soap3Index * index,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                        samfile_t ** samOutputDPFilePtrs,
                        bool shortDnaLength )
        {
            MC_MemberCopy4 ( this->, , unalignedReads, seedPool, insert_high, insert_low );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy4 ( this->, , queryNames, queryComments, upkdQualities, maxReadLength );
            MC_MemberCopy3 ( this->, , dpParameters, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy ( this->, , dpInfoForReads );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy3 ( this->, , outputFormat, samOutputDPFilePtrs, samOutputDPFilePtr );
            MC_MemberCopy ( this->, , shortDnaLength );
            soap3Wrapper =
                new SOAP3Wrapper<void> ( index);
        }
        ~DeepDPWrapper()
        {
            delete soap3Wrapper;
        }

        void seeding_mmp ( QueryIDStream * inputStream,
                           CandidateStream * canStream, QueryIDStream * unseededIDStream )
        {
            double startTime = setStartTime();
            PairEndSeedingEngine::
            performMmpSeeding (
                /* input */
                inputStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                soap3Wrapper, index, seedPool,
                /* output */
                canStream,
                unseededIDStream );
            fprintf(stderr, "[%s] Seeding Time: %.3f\n", __func__, getElapsedTime(startTime));
            //fprintf(stderr, "Unseeded: %d\n", unseededIDStream->data->size());
        }

        void alignment ( CandidateStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            PairEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters, shortDnaLength,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                queryNames, queryComments, origReadIDs, upkdQualities,
                dpInfoForReads,
                index,
                alignmentType,
                accumReadNum, outputFormat, samOutputDPFilePtr, seedPool, samOutputDPFilePtrs,
                /* output */
                unalignedIDStream, alignedRead, alignment );
            numDPAlignedRead += alignedRead;
            numDPAlignment += alignment;
        }

        void deepDP_MmpSeeding ( QueryIDStream * inputStream,
                                 QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            CandidateStream * canStream = new CandidateStream;
            seeding_mmp ( inputStream, canStream, unseededIDStream );
            //fprintf(stderr, "Number of DP instances: %d\n", (canStream->data).size());
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void DPExtension ( QueryIDStream * inputStream,
                          QueryIDStream * unalignedIDStream, PairInputForDPExtensionArrays * alignedShortPairs )
        {
            CandidateStream * canStream = new CandidateStream;
            vector<CandidateInfo> * candidateRegion = new vector<CandidateInfo>; 
            AlgnmtFlags * algnmtFlags = new AlgnmtFlags;
            
            for ( int i=0; i < alignedShortPairs->numArrays; ++i)
            {
                for ( uint j = 0; j < alignedShortPairs->inputArrays[i]->canRegionsSize; ++j )
                {
                    CandidateInfo tmp;
                    tmp.readIDLeft = alignedShortPairs->inputArrays[i]->canRegions[j].readIDLeft;
                    tmp.pos[0] = alignedShortPairs->inputArrays[i]->canRegions[j].pos[0];
                    tmp.pos[1] = alignedShortPairs->inputArrays[i]->canRegions[j].pos[1];
                    candidateRegion->push_back ( tmp ) ;
                }
            }
            canStream->append( candidateRegion, algnmtFlags );
            delete algnmtFlags;
            delete candidateRegion;
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DP2OutputUnalignedReads (
                unalignedIDStream, seedPool,
                queries, upkdReadLengths, maxReadLength,
                index, peStrandLeftLeg, peStrandRightLeg,
                queryNames, queryComments, origReadIDs, upkdQualities,
                accumReadNum, outputFormat, samOutputDPFilePtr
            );
        }
        void outputInvalidPair (QueryIDStream * unalignedIDStream, 
                                                    unsigned int maxHitNumForDP, unsigned int maxHitNumForDP2,
                                                    ReadInputForDP * dpInput, ReadInputForDP * dpInputForNewDefault,
                                                    ReadInputForDP * otherSoap3Result, BothUnalignedPairs * bothUnalignedPairs )
        {
            DPExtensionPassUnalignedReads (
                unalignedIDStream,
                queries, upkdReadLengths, maxReadLength,
                index, peStrandLeftLeg, peStrandRightLeg,
                queryNames, queryComments, origReadIDs, upkdQualities,
                accumReadNum, outputFormat, samOutputDPFilePtr,
                maxHitNumForDP, maxHitNumForDP2,
                dpInput, dpInputForNewDefault,
                otherSoap3Result, bothUnalignedPairs
            );
        }

        void runForMmpSeeding()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            QueryIDStream * input = new QueryIDStream ( unalignedReads );
            QueryIDStream * unaligned_round1 = new QueryIDStream;
            deepDP_MmpSeeding ( input, unaligned_round1, unaligned_round1 );
            delete input;
            // fprintf(stderr, "total unaligned: %d\n", (int) unaligned_round1->data->size());

            // prepare for next round
            resetBothUnalignedPairsArrays ( unalignedReads );
            free (unalignedReads->array[0]->readIDs);
            int new_size = unaligned_round1->data->size();
            unalignedReads->array[0]->readIDs = (uint *) malloc ( new_size * sizeof (uint)); 
            for ( int i=0;i< unaligned_round1->data->size(); ++i)
            {
                unalignedReads->array[0]->readIDs[i] = (*(unaligned_round1->data))[i];
            }
            unalignedReads->array[0]->size = new_size;
            unalignedReads->array[0]->totalNum = unaligned_round1->data->size();
            // uint * auxArray = ( uint * ) malloc ( new_size * sizeof(uint) );
            kx::radix_sort(unalignedReads->array[0]->readIDs, unalignedReads->array[0]->readIDs+unalignedReads->array[0]->totalNum);
            // MC_RadixSortForArray_32_16( unalignedReads->array[0]->readIDs, auxArray, new_size );
            // free ( auxArray );
            //for (int i=1;i<new_size;++i)
            //    if (unalignedReads->array[0]->readIDs[i] -2 != unalignedReads->array[0]->readIDs[i-1])
            //        fprintf(stderr, "%d\n",unalignedReads->array[0]->readIDs[i]-2);
            // end of preparation
            // outputUnaligned ( unaligned_round1 );
            
            delete unaligned_round1;
        }


};

//////////////////////////////////////////////////////////////////////////////////////
// The functions are for the DP designed for pairs of which both ends are unaligned //
//////////////////////////////////////////////////////////////////////////////////////

void DPForUnalignPairs2 ( BothUnalignedPairsArrays * unalignedReads, SeedPool* seedPool, int insert_high, int insert_low,
                          unsigned int * queries, unsigned int * upkdReadLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                          unsigned int maxReadLength,
                          int peStrandLeftLeg, int peStrandRightLeg,
                          DPInfoForReads * dpInfoForReads,
                          Soap3Index * index,
                          int alignmentType, DPParameters * dpParameters,
                          unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                          unsigned int accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, samfile_t ** samOutputDPFilePtrs)
{
    DeepDPWrapper
    deepDPWrapper (
        unalignedReads, seedPool, insert_high, insert_low,
        queries, upkdReadLengths, origReadIDs,
        queryNames, queryComments, upkdQualities, maxReadLength,
        peStrandLeftLeg, peStrandRightLeg,
        dpInfoForReads,
        index,
        alignmentType, dpParameters,
        accumReadNum, outputFormat, samOutputDPFilePtr, samOutputDPFilePtrs, false );
    if ( dpParameters->seedProperties.mmpSeedingScheme )
        { deepDPWrapper.runForMmpSeeding(); }
    else {
        assert(false);
    }
    numDPAlignedRead = deepDPWrapper.numDPAlignedRead;
    numDPAlignment = deepDPWrapper.numDPAlignment;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The functions are for the DP designed for short read pairs of which both ends are aligned in SOAP_SEED_LEN but fail in popcount              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




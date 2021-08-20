/*
 *
 *    DV-SemiDP.cu
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

#include "DV-SemiDP.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <stdlib.h>
#include <parallel/algorithm>
#include <vector>
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

bool compStartPos(const SingleAlgnmtResult & a, const SingleAlgnmtResult &b) {
    if (a.readID != b.readID) 
        return a.readID < b.readID;
    else if (a.score != b.score) 
        return a.score > b.score;
    else
        return a.startPos < b.startPos;
}


using namespace DP_Space;
class SemiDPWrapper
{
        Soap3Index * index;
        ReadInputForDPArrays * readInputArrays;
        SingleDP_Space::AlgnmtResultStream* singleDPResultStream;
        BothUnalignedPairs * unalignedReads;
        int insert_high, insert_low;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char * upkdQualities;
        char ** queryNames;
        char ** queryComments;
        DPInfoForReads * dpInfoForReads;
        uint maxReadLength;
        int peStrandLeftLeg, peStrandRightLeg;
        int alignmentType;
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;
        samfile_t ** samOutputDPFilePtrs;

        DPParameters * dpParameters;
        SOAP3Wrapper<void> * soap3Wrapper;

    public:
        uint numDPAlignedRead, numDPAlignment;

        SemiDPWrapper ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char ** queryNames, char ** queryComments, char * upkdQualities, uint maxReadLength,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        Soap3Index * index,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                        samfile_t ** samOutputDPFilePtrs,
                        BothUnalignedPairs * unalignedReads )
        {
            MC_MemberCopy3 ( this->, , readInputArrays, dpParameters, unalignedReads );
            MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy4 ( this->, , queryNames, queryComments, upkdQualities, maxReadLength );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy2 ( this->, , outputFormat, samOutputDPFilePtr );
            MC_MemberCopy  ( this->, , samOutputDPFilePtrs );
            soap3Wrapper =
                new SOAP3Wrapper<void> ( index);
        }
        SemiDPWrapper ( SingleDP_Space::AlgnmtResultStream* singleDPResultStream, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char ** queryNames, char ** queryComments, char * upkdQualities, uint maxReadLength,
                        DPInfoForReads * dpInfoForReads,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        Soap3Index * index,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                        samfile_t ** samOutputDPFilePtrs,
                        BothUnalignedPairs * unalignedReads )
        {
            readInputArrays = NULL;
            MC_MemberCopy3 ( this->, , singleDPResultStream, dpParameters, unalignedReads );
            MC_MemberCopy ( this->, , dpInfoForReads );
            MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy4 ( this->, , queryNames, queryComments, upkdQualities, maxReadLength );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy3 ( this->, , outputFormat, samOutputDPFilePtrs, samOutputDPFilePtr );
            soap3Wrapper = NULL;
        }
        ~SemiDPWrapper()
        {
            if ( soap3Wrapper )
            { delete soap3Wrapper; }
        }

        void alignment ( HalfEndOccStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, numAlignment = 0;
            HalfEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                dpInfoForReads,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                queryNames, queryComments, origReadIDs, upkdQualities,
                index,
                alignmentType,
                accumReadNum, outputFormat, samOutputDPFilePtr, samOutputDPFilePtrs,
                /* output */
                unalignedIDStream, alignedRead, numAlignment );
            numDPAlignedRead += alignedRead;
            numDPAlignment += numAlignment;
        }

        void passUnalignedToDeepDP ( QueryIDStream * unalignedIDStream )
        {
            for ( int i = 0; i < unalignedIDStream->data->size(); i++ )
            {
                addReadIDToBothUnalignedPairs ( unalignedReads, ( * ( unalignedIDStream->data ) ) [i] );
            }
        }


        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DeepDP_Space::DP2OutputUnalignedReads (
                unalignedIDStream, NULL,
                queries, upkdReadLengths, maxReadLength,
                index, peStrandLeftLeg, peStrandRightLeg,
                queryNames, queryComments, origReadIDs, upkdQualities,
                accumReadNum, outputFormat, samOutputDPFilePtr
            );
        }

        void runForSingleDPResultStream()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            omp_set_num_threads(dpParameters->numOfCPUThreads);
            for (int i=0;i<singleDPResultStream->dpSResult.size();i++) __gnu_parallel::sort(singleDPResultStream->dpSResult[i][0].begin(), singleDPResultStream->dpSResult[i][0].end(), compStartPos);
            HalfEndOccStream * canStream = new HalfEndOccStream ( singleDPResultStream );
            QueryIDStream * unalignedIDStream = new QueryIDStream;
            // printf("Go Aligning!!\n");
            alignment ( canStream, unalignedIDStream );
            fprintf(stderr, "[Default DP] Total unaligned reads: %d\n", (int) unalignedIDStream->data->size());
            passUnalignedToDeepDP ( unalignedIDStream );
            delete canStream;
            delete unalignedIDStream;
        }

};

// TODO
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN start here ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

struct CompareIDAndScore {
    bool operator() (const SingleAlgnmtResult &a, const SingleAlgnmtResult &b) {
        if (a.readID < b.readID) return true;
        else if (a.readID > b.readID) return false;
        return a.score > b.score;
    }
};

// To perform semi-global DP for single end DP alignment
void  semiGlobalDPForSingleEndDpAlignment ( SingleDP_Space::AlgnmtResultStream* singleDPResultStream , int insert_high, int insert_low,
                       unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char ** queryNames, char ** queryComments, char * upkdQualities,
                       DPInfoForReads * dpInfoForReads,
                       unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                       Soap3Index * index,
                       int alignmentType, DPParameters * dpParameters,
                       unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                       unsigned int accumReadNum, int outputFormat,
                       samfile_t * samOutputDPFilePtr, samfile_t ** samOutputDPFilePtrs,
                       BothUnalignedPairs * &unalignedReads )
{
    // singleDPResultStream: input single end result for DP
    // insert_high : maximum value of insert size
    // insert_low: minimum value of insert size
    // upkdQueries: the sequence of all reads. The first character of the read with ID "i"
    //              is on the position "i*maxReadLength" of the array
    // readLengths: the read length of all reads
    // peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
    // peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
    // BWT: bwt structure with SA table inside
    // HSP: hsp structure with packed sequence inside
    // alignmentType: 1: All valid alignment; 2: all best alignment
    // numDPResult: the total number of alignment results outputted by DP

    SingleDP_Space::AlgnmtResultStream * sortedSingleDPResultStream = new SingleDP_Space::AlgnmtResultStream;
    SingleDP_Space::SingleDPResultBatch * tmpResultBatch = new SingleDP_Space::SingleDPResultBatch;
    SingleDP_Space::SingleDPResultBatch * auxResultBatch = new SingleDP_Space::SingleDPResultBatch;

    // The following part is the old code for sorting Single DP results first by first ReadID and second the score
    // It is very slow in some cases and I (Dinghua had refactored) this part, but this should be reviewed by Chun

    // for ( int i=0; i<singleDPResultStream->dpSResult.size(); ++i )
    // {
    //     SingleDP_Space::SingleDPResultBatch & tmp = (*(singleDPResultStream->dpSResult[i]));
    //     for ( int j=0; j<tmp.size(); ++j)
    //     {
    //         tmpResultBatch->push_back ( tmp[j] );
    //     }
    // }
    // int dpResultLen = tmpResultBatch->size();
    // auxResultBatch->resize ( dpResultLen );
    // MC_RadixSort_32_16 ( tmpResultBatch->begin(), score, auxResultBatch->begin(), dpResultLen );
    // for ( int i=0;i<dpResultLen>>1;++i )
    // { 
    //     SingleAlgnmtResult tmpSingle = (*tmpResultBatch)[i];
    //     (*tmpResultBatch)[i] = (*tmpResultBatch)[dpResultLen-i-1];
    //     (*tmpResultBatch)[dpResultLen-i-1] = tmpSingle;
    // }
    // MC_RadixSort_32_16 ( tmpResultBatch->begin(), readID, auxResultBatch->begin(), dpResultLen );
    // sortedSingleDPResultStream->dpSResult.push_back( tmpResultBatch );

    for ( int i=0; i<singleDPResultStream->dpSResult.size(); ++i ) {
        SingleDP_Space::SingleDPResultBatch & tmp = (*(singleDPResultStream->dpSResult[i]));
        for ( int j=0; j<tmp.size(); ++j)
        {
            tmpResultBatch->push_back ( tmp[j] );
        }
    }

    omp_set_num_threads(dpParameters->numOfCPUThreads);
    __gnu_parallel::sort(tmpResultBatch->begin(), tmpResultBatch->end(), CompareIDAndScore());
    sortedSingleDPResultStream->dpSResult.push_back( tmpResultBatch );

    // ---- end of refactoring ----

    SemiDPWrapper
    semiDPWrapper (
        sortedSingleDPResultStream , insert_high, insert_low,
        queries, readLengths, origReadIDs,
        queryNames, queryComments, upkdQualities, maxReadLength,
        dpInfoForReads,
        peStrandLeftLeg, peStrandRightLeg,
        index,
        alignmentType, dpParameters,
        accumReadNum, outputFormat, samOutputDPFilePtr, samOutputDPFilePtrs,
        unalignedReads );
    semiDPWrapper.runForSingleDPResultStream();
    numDPAlignedRead = semiDPWrapper.numDPAlignedRead;
    numDPAlignment = semiDPWrapper.numDPAlignment;
    delete auxResultBatch;
    for ( int i = 0; i < sortedSingleDPResultStream->dpSResult.size(); i++ )
    {
        // Belwo comment is a remind that the cigar string should not be cleaned
        // because "sortedSingleDPResultStream" is a copy of singleDPResultStream
        // freeing the cigar string would be done when deleting singleDPResultStream
        //
        // Please do not amend this comment
        /**************************************************************************************************
        SingleDP_Space::SingleDPResultBatch & resultBatch = * ( sortedSingleDPResultStream->dpSResult[i] );
        for ( int j = 0; j < resultBatch.size(); j++ )
        {
            free ( resultBatch[j].cigarString );
        }
        ****************************************************************************************************/
        delete sortedSingleDPResultStream->dpSResult[i];
    }
    sortedSingleDPResultStream->dpSResult.clear();
    delete sortedSingleDPResultStream;
    return;
}




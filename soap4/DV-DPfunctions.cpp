/*
 *
 *    DV-DPfunctions.cu
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
#include <queue>
#include <omp.h>
#include "DV-DPfunctions.h"
#include "OutputDPResult.h"
#include <assert.h>
#include <pthread.h>
#include <algorithm>
#include <parallel/algorithm>
#include "kxsort.h"
#include <tuple>
using namespace std;

const unsigned long long kMaxULL = 0xFFFFFFFFFFFFFFFFULL;

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// For output ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

AlgnmtFlags::AlgnmtFlags ( uint range )
{
    size = ( range + 31 ) / 32;
    MC_CheckMalloc ( flags,   uint,   size );

    for ( int i = 0; i < 32; i++ )
    {
        MASK[i] = 1 << i;
    }

    pthread_mutex_init ( &occupy_mutex, NULL );
    clear();
}

void AlgnmtFlags::clear()
{
    memset ( flags, 0, size * sizeof ( uint ) );
}

void AlgnmtFlags::increaseSize ( uint newSize )
{
    uint * oldFlags = flags;
    uint oldSize = size;
    size  = newSize;
    MC_CheckMalloc ( flags,   uint,   size );
    memcpy ( flags, oldFlags, oldSize * sizeof ( uint ) );
    memset ( flags + oldSize, 0, ( newSize - oldSize ) * sizeof ( uint ) );
    free ( oldFlags );
}

void AlgnmtFlags::set ( int readID )
{
    pthread_mutex_lock ( &occupy_mutex );
    uint offset = readID >> 5;

    if ( offset >= size )
    {
        uint newSize = size * 2;

        while ( offset >= newSize )
        { newSize *= 2; }

        increaseSize ( newSize );
    }

    flags[offset] |= MASK[readID & 0x1F];
    pthread_mutex_unlock ( &occupy_mutex );
}

#define AlgnmtFlags_Get(in_flag, int32Offset, diff) { \
        uint flag = in_flag; \
        if (flag != 0) { \
            int offset = int32Offset << 5; \
            for (int j = 0; flag != 0; j++) { \
                if (flag & 1) { \
                    diff->push_back(offset + j); \
                } \
                flag >>= 1; \
            } \
        } \
    }

void AlgnmtFlags::get ( vector<int> * diff )
{
    for ( int i = 0; i < size; i++ )
    {
        AlgnmtFlags_Get ( flags[i], i, diff );
    }
}

inline void AlgnmtFlags::reserveSize ( AlgnmtFlags * algnFlags )
{
    if ( size < algnFlags->size )
    {
        this->increaseSize ( algnFlags->size );
    }
    else if ( size > algnFlags->size )
    {
        algnFlags->increaseSize ( size );
    }
}

void AlgnmtFlags::getXOR ( AlgnmtFlags * algnFlags, vector<int> * diff )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        AlgnmtFlags_Get ( flags[i] ^ algnFlags->flags[i], i, diff );
    }
}

void AlgnmtFlags::AND ( AlgnmtFlags * algnFlags )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        flags[i] &= algnFlags->flags[i];
    }
}

void AlgnmtFlags::XOR ( AlgnmtFlags * algnFlags )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        flags[i] ^= algnFlags->flags[i];
    }
}

AlgnmtFlags::~AlgnmtFlags()
{
    free ( flags );
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Alignment modules //////////////////////////////////////
/////////////////// The following code better be placed in a seperate file ///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// standard space ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


template <>
int isValid ( AlgnmtDPResult & a )
{
    return ( a.whichFromDP < 2 );
}

template <>
int ScoreCompare ( AlgnmtDPResult & a, AlgnmtDPResult & b )
{
#define MC_DPScoreCompare_SetValue(result, aligned, mismatch, score) { \
        if (result.whichFromDP == 0) { \
            aligned = 1; \
            mismatch = result.score_2; \
            score = result.score_1; \
        } \
        else \
            if (result.whichFromDP == 1) { \
                aligned = 1; \
                mismatch = result.score_1; \
                score = result.score_2; \
            } \
            else { \
                aligned = 0; \
                if (result.algnmt_1 != kMaxULL) \
                    mismatch = result.score_1; \
                else \
                    mismatch = result.score_2; \
                score = 0; \
            } \
    }
    uint aligned_a, mismatch_a, score_a;
    uint aligned_b, mismatch_b, score_b;
    MC_DPScoreCompare_SetValue ( a, aligned_a, mismatch_a, score_a );
    MC_DPScoreCompare_SetValue ( b, aligned_b, mismatch_b, score_b );
    uint64 value_a = ( ( uint64 ) aligned_a << 63 ) | ( ( uint64 ) ( 0x1FFFFFFF - mismatch_a ) << 32 ) | ( score_a + 0x1FFFFFFF );
    uint64 value_b = ( ( uint64 ) aligned_b << 63 ) | ( ( uint64 ) ( 0x1FFFFFFF - mismatch_b ) << 32 ) | ( score_b + 0x1FFFFFFF );

    if ( value_a > value_b )
    { return 1; }
    else if ( value_a < value_b )
    { return -1; }
    else
    { return 0; }
}

template <>
int ScoreCompare ( SingleAlgnmtResult & a, SingleAlgnmtResult & b )
{
    if ( a.score > b.score )
    { return 1; }
    else if ( a.score < b.score )
    { return -1; }
    else
    { return 0; }
}

template <>
int ScoreCompare ( DeepDPAlignResult & a, DeepDPAlignResult & b )
{
    int score_a = a.score_1 + a.score_2;
    int score_b = b.score_1 + b.score_2;

    if ( score_a > score_b )
    { return 1; }
    else if ( score_a < score_b )
    { return -1; }
    else
    { return 0; }
}

template <>
bool ResultCompare ( const AlgnmtDPResult & a, const AlgnmtDPResult & b )
{
    return make_tuple(a.algnmt_1, a.algnmt_2, a.score_1, a.score_2) < make_tuple(b.algnmt_1, b.algnmt_2, b.score_1, b.score_2);
}
template <>
bool ResultCompare ( const SingleAlgnmtResult & a, const SingleAlgnmtResult & b )
{
    return make_tuple(a.algnmt, a.score) < make_tuple(b.algnmt, b.score);
}
template <>
bool ResultCompare ( const DeepDPAlignResult & a, const DeepDPAlignResult & b )
{
    return make_tuple(a.algnmt_1, a.algnmt_2, a.score_1, a.score_2) < make_tuple(b.algnmt_1, b.algnmt_2, b.score_1, b.score_2);

    //return ( a.algnmt_1 == b.algnmt_1 ? ( a.score_1 == b.score_1 ? ( a.algnmt_2 == b.algnmt_2 ? a.score_2 > b.score_2 : a.algnmt_2 < b.algnmt_2 ) : a.score_1 > b.score_1 ) : a.algnmt_1 < b.algnmt_1 );
}

//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// single-dp space //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

using namespace SingleDP_Space;
#define DPS_SEEDING_BATCH_SIZE 256 * 1024
#define DPS_MARGIN(l) ((l>100) ? 30 : 25)

CandidateStream::CandidateStream()
{
    pthread_mutex_init ( &occupy_mutex, NULL );
}
void CandidateStream::append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags )
{
    pthread_mutex_lock ( &occupy_mutex );

    for ( vector<CandidateInfo>::iterator it = canInfo->begin();
            it < canInfo->end(); ++it )
    {
        data.push_back ( *it );
        alignFlags->set ( it->readID );
    }

    pthread_mutex_unlock ( &occupy_mutex );
}

SingleEndSeedingEngine::SingleEndSeedingEngine() {}

struct LongerSeedLength {
    bool operator() (const SeedPos &a, const SeedPos &b) {
        return a.seedAlignmentLength > b.seedAlignmentLength;
    }
};

vector<CandidateInfo> * SingleEndSeedingEngine::singleMerge (
    SingleDP_Space::SeedPos * readPos, QueryIDStream * eliminatedReadIDs
)
{
    vector<CandidateInfo> * canInfo = new vector<CandidateInfo>();
    SingleDP_Space::SeedPos * p = readPos;
    uint seedReadCnt = 0;
    uint lastReadID = 0x7FFFFFFF;
// TODO: omp this part; may be significant
    while ( p->readID != 0x7FFFFFFF )
    {
        uint readID = p->readID;
        bool keepReads = false;
        size_t oldSize = canInfo->size();

        for (; p->readID == readID; p++ )
        {
            if (p->seedAlignmentLength < 17) {
                continue;
            }

            keepReads = true;
            canInfo->push_back(*p);

            while ((p+1)->readID == readID) {
                if ((p+1)->pos < canInfo->back().pos + DPS_DIVIDE_GAP && canInfo->back().strand == (p+1)->strand){
                    if ((p+1)->seedAlignmentLength > canInfo->back().seedAlignmentLength)
                        canInfo->back() = *(p+1);
                } else {
                    break;
                }
                ++p;
            }
        }
        if (!keepReads && eliminatedReadIDs != NULL) {
            eliminatedReadIDs->data->push_back(readID);
        }
        ++seedReadCnt;

        std::sort(canInfo->begin() + oldSize, canInfo->end(), LongerSeedLength());
        if (oldSize < canInfo->size()) {
            while (canInfo->back().seedAlignmentLength < (*canInfo)[oldSize].seedAlignmentLength * 0.6) {
                canInfo->pop_back();
            }
        }
    }        
    return canInfo;
}
// ****
SingleEndAlignmentEngine::SingleEndAlgnBatch::SingleEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths,
    DPInfoForReads * dpInfoForReads
)
{
    MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy3 ( this->, , queries, inputMaxReadLength, upkdLengths );
    MC_MemberCopy ( this->, , dpInfoForReads );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->wordPerOldQuery   = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery      = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA        = MC_CeilDivide16 ( maxDNALength );
    this->packedDNA         = index->sraIndex->hsp->packedDNA;
    this->fullDNALength     = index->sraIndex->hsp->dnaLength;
    this->index             = index;
    MC_CheckMalloc ( canInfos,            CandidateInfo,  batchSize );
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( packedDNASeq,        uint,           batchSize * MC_CeilDivide16 ( maxDNALength ) );
    MC_CheckMalloc ( packedReadSeq,       uint,           batchSize * MC_CeilDivide16 ( maxReadLength ) );
    MC_CheckMalloc ( scores,              int,            batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );
    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( peLeftAnchorLocs,    uint,           batchSize );
    MC_CheckMalloc ( peRightAnchorLocs,   uint,           batchSize );
    MC_CheckMalloc ( hitLocs,             uint,           batchSize );
    MC_CheckMalloc ( pattern,             uchar,          batchSize * patternLength );
    MC_CheckMalloc ( maxScoreCounts,      uint,           batchSize );
    clear();
}

SingleEndAlignmentEngine::SingleEndAlgnBatch::~SingleEndAlgnBatch()
{
    free ( canInfos );
    free ( DNALengths );
    free ( lengths );
    free ( packedDNASeq );
    free ( packedReadSeq );
    free ( scores );
    free ( cutoffThresholds );
    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( peLeftAnchorLocs );
    free ( peRightAnchorLocs );
    free ( hitLocs );
    free ( pattern );
    free ( maxScoreCounts );
}

void SingleEndAlignmentEngine::SingleEndAlgnBatch::clear()
{
    numOfThreads = 0;
}

int SingleEndAlignmentEngine::SingleEndAlgnBatch::pack (
    CandidateInfo & canInfo
)
{
    if ( numOfThreads >= batchSize )
    {
        return 0;
    }

    uint readID = canInfo.readID;
    uint readLength = upkdLengths[readID];
    int margin = DPS_MARGIN ( readLength );
    unsigned long long DNAStart = canInfo.pos - margin;

    if ( DNAStart >= fullDNALength )
    {
        DNAStart = 0;
    }

    uint DNALength = readLength + margin * 2;

    if ( DNAStart + DNALength > fullDNALength )
    {
        DNALength = fullDNALength - DNAStart;
    }

    packRead ( packedReadSeq, numOfThreads,
               readID, readLength,
               canInfo.strand );
    repackDNA ( packedDNASeq, numOfThreads,
                packedDNA, DNAStart, DNALength );
    softClipLtSizes[numOfThreads] = ( canInfo.strand == 1 ) ?
                                    softClipLeft : softClipRight;
    softClipRtSizes[numOfThreads] = ( canInfo.strand == 1 ) ?
                                    softClipRight : softClipLeft;
    peLeftAnchorLocs[numOfThreads] = maxDNALength;
    peRightAnchorLocs[numOfThreads] = 0;
    DNALengths[numOfThreads] = DNALength;
    lengths[numOfThreads] = readLength;
    cutoffThresholds[numOfThreads] = std::max(DP_SCORE_THRESHOLD_RATIO * readLength, DP_SCORE_THRESHOLD_LOWER_BOUND);
    canInfo.pos = DNAStart;
    canInfos[numOfThreads] = canInfo;
    ++numOfThreads;
    return 1;
}

inline void SingleEndAlignmentEngine::SingleEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}

inline void SingleEndAlignmentEngine::SingleEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, unsigned long long start, uint length
)
{
    size_t dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    uint *src = seq + start / 16;
    uint src_ofs = start % 16;
    uint *dst = packedSeq + dnaTPARA;
    uint dst_ofs = 1;
    *dst = 0;

    // we need a 0 at the beginning to reserve the first column of DP table
    while (length > 0) {
        uint len = min(min(16 - dst_ofs, 16 - src_ofs), length);
        *dst |= *src << (src_ofs * 2) >> (32 - len*2) << (32-(dst_ofs+len)*2);
        length -= len;
        dst_ofs += len;
        src_ofs += len;

        if (src_ofs == 16) { ++src; src_ofs = 0; }
        if (dst_ofs == 16) { dst += 32; dst_ofs = 0; *dst = 0; }
    }
}

// ****
void SingleEndAlignmentEngine::SingleEndAlgnThreadContext::init ( SingleEndAlgnBatch * batch )
{
    int batchSize = engine->DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                             engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}

void SingleEndAlignmentEngine::SingleEndAlgnThreadContext::freeMemory()
{
    semiGlobalAligner.freeMemory();
    delete batch;
}

// ****
SingleDP_Space::AlgnmtResultStream::AlgnmtResultStream()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}

SingleDP_Space::AlgnmtResultStream::~AlgnmtResultStream()
{
    for ( int i = 0; i < dpSResult.size(); i++ )
    {
        SingleDPResultBatch & resultBatch = * ( dpSResult[i] );

        for ( int j = 0; j < resultBatch.size(); j++ )
        {
            free ( resultBatch[j].cigarString );
        }

        delete dpSResult[i];
    }

    dpSResult.clear();
}

// ****
void SingleEndAlignmentEngine::performAlignment (
    uint & numDPAlignedRead, uint & numDPAlignment
)
{
    /* initialize */
    algnBatchCount = 0;
    dpSAlignedRead = 0;
    dpSAlignment = 0;
    lastReadID = -1;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    outputBuf = new OutputBuffer<SingleAlgnmtResult>();
    outputBuf->setAlignmentType ( alignmentType );
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = maxReadLength + 2 * DPS_MARGIN ( inputMaxReadLength ) + 8;
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DPS_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );
    algnSwapBatch =
        new SingleEndAlgnBatch ( DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                 maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                 index, queries, inputMaxReadLength, upkdReadLengths,
                                 dpInfoForReads );
    algnThreadContext = new SingleEndAlgnThreadContext[dpPara->numOfCPUThreads];

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        SingleEndAlgnBatch * batch =
            new SingleEndAlgnBatch ( DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                     maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                     index, queries, inputMaxReadLength, upkdReadLengths,
                                     dpInfoForReads );
        algnThreadContext[i].init ( batch );
    }

    outputThreadDelegator.init ( 1, DPSOutputThread,
                                 NULL, DPSOutputThreadFinalize );
    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, algnmtCPUThread );
    /* perform alignment */
    int threadId;
    void * empty;

    for ( uint i = 0; i < canStream->data.size(); i++ )
    {
        CandidateInfo & info = canStream->data[i];
        inputFlags->set ( info.readID );

        if ( !algnSwapBatch->pack ( info ) )
        {
            // launch one batch
            threadId = algnmtCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
            algnSwapBatch->clear();
            algnSwapBatch->pack ( info );
        }
    }

    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
    }

    /* finalize */
    algnmtCPUThreadDelegator.finalize();
    outputThreadDelegator.finalize();
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory();
    }

    delete[] algnThreadContext;
    delete outputBuf;
    // resultStream is used as an output
    // delete resultStream;
    int cnt=0;
    for ( int i=0;i<resultStream->dpSResult.size(); ++i)
        cnt += (*((resultStream->dpSResult)[i])).size();
    numDPAlignedRead = this->dpSAlignedRead;
    numDPAlignment = this->dpSAlignment;
}

void SingleEndAlignmentEngine::performAlignment (
    /* input */
    CandidateStream   *   canStream,
    DPParameters     *    dpPara,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    DPInfoForReads * dpInfoForReads,
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment,
    AlgnmtResultStream * &resultStream
)
{
    engine = new SingleEndAlignmentEngine();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, queryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy ( engine->, , queryComments);
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , dpInfoForReads );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy3 ( engine->, , accumReadNum, outputFormat, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    resultStream = engine->resultStream;
    delete engine;
}

SingleEndAlignmentEngine * SingleEndAlignmentEngine::engine;

void SingleDP_Space::algnmtCPUThread ( int threadId, void  *& empty )
{
//printf("singledp cpu thread %d start\n",threadId);
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    SingleEndAlignmentEngine::SingleEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++;
    sem_post ( & ( engine->algnThreadContext[threadId].ACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    //engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
    //sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    engine->algnThreadContext[threadId].semiGlobalAligner.performAlignment (
        batch->packedDNASeq, batch->DNALengths,
        batch->packedReadSeq, batch->lengths,
        batch->cutoffThresholds, batch->scores, batch->hitLocs,
        batch->maxScoreCounts,
        batch->pattern, batch->numOfThreads,
        batch->softClipLtSizes, batch->softClipRtSizes,
        batch->peLeftAnchorLocs, batch->peRightAnchorLocs );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    // Rearrange result and Output
    SingleDPResultBatch * resultBatch = new SingleDPResultBatch;
    for ( int i = 0; i < batch->numOfThreads; i++ )
    {
        if ( batch->scores[i] >= batch->cutoffThresholds[i] )
        {
            CigarStringEncoder<void> encoder;
            uchar lastType = 'N';

            for ( uchar * p = batch->pattern + i * engine->patternLength; *p != 0; p++ )
            {
                if ( *p == 'V' )
                {
                    encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                }
                else
                {
                    encoder.append ( *p, 1 );
                    lastType = *p;
                }
            }

            SingleAlgnmtResult result;
            result.readID = batch->canInfos[i].readID;
            result.strand = batch->canInfos[i].strand;
            result.seedAlignmentLength = batch->canInfos[i].seedAlignmentLength;
            result.algnmt = batch->canInfos[i].pos + batch->hitLocs[i];
            result.score = batch->scores[i];
            result.startPos = batch->canInfos[i].pos;
            result.refDpLength = batch->DNALengths[i];
            result.peLeftAnchor = batch->peLeftAnchorLocs[i];
            result.peRightAnchor = batch->peRightAnchorLocs[i];
            encoder.encodeCigarString ( openGapScore, extendGapScore );
            result.cigarString = encoder.cigarString;
            int L = batch->lengths[i] - encoder.charCount['I'] - encoder.charCount['S'];
            int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[i] ) /
                                ( matchScore - mismatchScore );
            result.editdist = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
            result.num_sameScore = batch->maxScoreCounts[i]; // TODO
            resultBatch->push_back ( result );
        }
    }

    // Output
    engine->algnThreadContext[threadId].resultBatch = resultBatch;
    int * pid = &threadId;
    engine->outputThreadDelegator.schedule ( pid );
    // fprintf(stderr, "ALgn CPU Thread done.\n");
    sem_wait ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
//printf("singledp cpu thread %d end\n",threadId);
}

void SingleDP_Space::DPSOutputThread ( int threadId, int *& pCallThreadId )
{
    int callThreadId = *pCallThreadId;
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    int batchID = engine->algnThreadContext[callThreadId].batchID;
    SingleDPResultBatch * resultBatch = engine->algnThreadContext[callThreadId].resultBatch;
    sem_post ( & ( engine->algnThreadContext[callThreadId].outputACKSem ) );
    vector<SingleDPResultBatch *> & dpResult = engine->resultStream->dpSResult;

    while ( dpResult.size() <= batchID )
    {
        dpResult.push_back ( NULL );
    }

    dpResult[batchID] = resultBatch;
/**/
#define MC_DPSOutputRead() { \
            engine->outputBuf->ready(); \
            if (engine->outputBuf->size > 0) { \
                engine->dpSAlignedRead += 1; \
                engine->dpSAlignment += engine->outputBuf->size; \
                engine->alignFlags->set(engine->lastReadID); \
            } \
        }
/**/
    uint numOut = engine->resultStream->numOut;

    while ( numOut < dpResult.size() && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        SingleDPResultBatch & batch = *dpResult[numOut];

        for ( int i = 0; i < batch.size(); i++ )
        {
            SingleAlgnmtResult & result = batch[i];
            int readID = result.readID;

            if ( readID != engine->lastReadID )
            {
                MC_DPSOutputRead();
                engine->outputBuf->clear();
                engine->lastReadID = readID;
            }
            
            engine->outputBuf->add ( result );
        }
        ++numOut;
    }

    engine->resultStream->numOut = numOut;
}

void SingleDP_Space::DPSOutputThreadFinalize()
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    MC_DPSOutputRead();
    engine->outputBuf->clear();
}

SingleAlgnmtResult * SingleDP_Space::fetchNextSingleDPResult ( SingleDP_Space::AlgnmtResultStream * singleDPResultStream, 
                                               uint &singleDPResultStreamIdx, uint &singleDPResultBatchIdx)
{
    while(true) {
        if ( singleDPResultStreamIdx < singleDPResultStream->dpSResult.size() )
        {
            SingleDPResultBatch & singleDPResultBatch = *(singleDPResultStream->dpSResult[singleDPResultStreamIdx]);
            if (singleDPResultBatch.size() == 0) {
                ++singleDPResultStreamIdx;
                singleDPResultBatchIdx = 0;
                continue;
            } else if ( singleDPResultBatchIdx < singleDPResultBatch.size() )
            {
                SingleAlgnmtResult * singleAlgnmtResult = &(singleDPResultBatch[singleDPResultBatchIdx]);
                ++singleDPResultBatchIdx;
                if ( singleDPResultBatchIdx == singleDPResultBatch.size() )
                { 
                    ++singleDPResultStreamIdx; 
                    singleDPResultBatchIdx = 0;
                }
                if (singleAlgnmtResult) 
                    return singleAlgnmtResult;
                else 
                    continue;
            } 
        } else return NULL;
    }
}

// both of singleDPResultStream and unalignedSingleReads are sorted by readID
void SingleDP_Space::DPSOutputUnpairedAlignment ( SingleDP_Space::AlgnmtResultStream * singleDPResultStream, UnalignedSingles * unalignedSingleReads,
                                 uint * queries, uint * upkdReadLengths, uint * origReadIDs, 
                                 char ** queryNames, char ** queryComments, char * upkdQualities, int inputMaxReadLength,
                                 uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                                 Soap3Index * index,
                                 int alignmentType,
                                 uint         &        dpSAlignedRead,
                                 uint         &        dpSAlignment   )
{
    OutputBuffer<SingleAlgnmtResult> * outputBuf = new OutputBuffer<SingleAlgnmtResult>();
    outputBuf->setAlignmentType ( alignmentType );
#define MC_DPSUnpairedOutputRead() { \
        outputBuf->ready(); \
        if (outputBuf->size > 0) { \
            outputDPSingleResult2( \
                                   outputBuf->elements, outputBuf->size, \
                                   queries, upkdReadLengths, origReadIDs, \
                                   queryNames, queryComments, upkdQualities, \
                                   inputMaxReadLength, accumReadNum, outputFormat, \
                                   samOutputDPFilePtr, index); \
            dpSAlignedRead += 1; \
            dpSAlignment += outputBuf->size; \
        } \
    }
    uint singleDPResultStreamIdx = 0, singleDPResultBatchIdx = 0, unalignedIter = 0;
    SingleAlgnmtResult * singleAlgnmtResult =  fetchNextSingleDPResult ( singleDPResultStream, singleDPResultStreamIdx, singleDPResultBatchIdx );
    QueryIDStream * unalignedIDStream = new QueryIDStream;
    while ( true )
    {
        while ( singleAlgnmtResult != NULL && singleAlgnmtResult->readID < unalignedSingleReads->readIDs[unalignedIter] )
            { singleAlgnmtResult = fetchNextSingleDPResult ( singleDPResultStream, singleDPResultStreamIdx, singleDPResultBatchIdx); }
        if ( singleAlgnmtResult == NULL )
            { break; }
        while ( unalignedIter < unalignedSingleReads->totalNum && unalignedSingleReads->readIDs[unalignedIter] < singleAlgnmtResult->readID )
        { 
            unalignedIDStream->data->push_back( unalignedSingleReads->readIDs[unalignedIter] );
            ++unalignedIter; 
        }
        if ( unalignedIter == unalignedSingleReads->totalNum )
            { break; }
        if ( singleAlgnmtResult->readID == unalignedSingleReads->readIDs[unalignedIter] )
        {   
            while (  singleAlgnmtResult != NULL && singleAlgnmtResult->readID == unalignedSingleReads->readIDs[unalignedIter] )
            {
                SingleAlgnmtResult & tmp = *singleAlgnmtResult ;
                outputBuf->add ( tmp );
                singleAlgnmtResult = fetchNextSingleDPResult ( singleDPResultStream, singleDPResultStreamIdx, singleDPResultBatchIdx);
            }           
            MC_DPSUnpairedOutputRead();
            outputBuf->clear();
            ++unalignedIter;
        }
        if ( unalignedIter == unalignedSingleReads->totalNum )
            { break; }
    }
    for ( ; unalignedIter < unalignedSingleReads->totalNum; ++unalignedIter )
    {
        unalignedIDStream->data->push_back( unalignedSingleReads->readIDs[unalignedIter] );
    }
    DPSOutputUnalignedReads ( unalignedIDStream, 
                              queries, upkdReadLengths, inputMaxReadLength,
                              index, 
                              queryNames, queryComments, origReadIDs, upkdQualities,
                              accumReadNum, outputFormat,
                              samOutputDPFilePtr );
    delete outputBuf;
    delete unalignedIDStream;
    return;
}


void SingleDP_Space::DPSOutputUnalignedReads (
    QueryIDStream * unalignedIDStream,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr
)
{
    // output unaligned result
#define MC_DPSOutputUnalgnRead() { \
        outputDPSingleResult2(buf, idx, \
                              queries, upkdReadLengths, origReadIDs, \
                              queryNames, queryComments, upkdQualities, \
                              inputMaxReadLength, accumReadNum, outputFormat, \
                              samOutputDPFilePtr, index); }
    SingleAlgnmtResult * buf;
    MC_CheckMalloc ( buf, SingleAlgnmtResult, 1024 );
    int idx = 0;

    for ( uint i = 0; i < unalignedIDStream->data->size(); i++ )
    {
        buf[idx].readID = ( * ( unalignedIDStream->data ) ) [i];
        buf[idx].algnmt = kMaxULL;
        buf[idx].cigarString = NULL;
        ++idx;

        if ( idx >= 1024 )
        {
            MC_DPSOutputUnalgnRead();
            idx = 0;
        }
    }

    if ( idx > 0 )
    { MC_DPSOutputUnalgnRead(); }

    free ( buf );
}



//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// default-dp space /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
using namespace DP_Space;

// ****
HalfEndOccStream::HalfEndOccStream ( ReadInputForDPArrays * input, BWT * bwt )
{
    this->data = input;
    this->bwt = bwt;
    arrayIndex = 0;
    iter_readInput = data->inputArrays[arrayIndex];
    iter_occ = iter_readInput->occ_list;
    end_occ = iter_occ + iter_readInput->occTotalNum;
    iter_sa = iter_readInput->sa_list;
    end_sa = iter_sa + iter_readInput->saRangeTotalNum;
    nextSAIndex = -1;
}

HalfEndOccStream::HalfEndOccStream ( SingleDP_Space::AlgnmtResultStream * input )
{
    this->data = NULL;

    this->singleDPData = input;
    this->singleDPDataArrayIndex=0;
    this->iter_singleDPData = 0;
    
}


int HalfEndOccStream::fetchNextOcc ( SRAOccurrence & occ )
{
    while ( true )
    {
#define SA2OCC() { \
        occ.readID = iter_sa->readID; \
        occ.mismatchCount = iter_sa->mismatchCount; \
        occ.strand = iter_sa->strand; \
        occ.ambPosition = BWTSaValue(bwt, nextSAIndex++); \
        if (nextSAIndex > iter_sa->saIndexRight) { \
            nextSAIndex = -1; \
            ++iter_sa; \
        } \
    }

        if ( nextSAIndex != -1 )
        {
            SA2OCC();
        }
        else if ( iter_occ == end_occ )
        {
            if ( iter_sa == end_sa )
            {
                ++arrayIndex;

                if ( arrayIndex >= data->numArrays )
                { return 0; }  // finished
                else
                {
                    iter_readInput = data->inputArrays[arrayIndex];
                    iter_occ = iter_readInput->occ_list;
                    end_occ = iter_occ + iter_readInput->occTotalNum;
                    iter_sa = iter_readInput->sa_list;
                    end_sa = iter_sa + iter_readInput->saRangeTotalNum;
                    continue;
                }
            }
            else
            {
                nextSAIndex = iter_sa->saIndexLeft;
                SA2OCC();
            }
        }
        else
        {
            if ( iter_sa == end_sa )
            {
                occ = * ( iter_occ++ );
            }
            else
            {
                if ( ( iter_occ->readID >> 1 ) < ( iter_sa->readID >> 1 ) )
                { occ = * ( iter_occ++ ); }
                else
                {
                    nextSAIndex = iter_sa->saIndexLeft;
                    SA2OCC();
                }
            }
        }

        return 1;
    }
}

int HalfEndOccStream::fetchNextSingleAlgnResult( SingleAlgnmtResult & singleAlgnmtResult)
{
    
    if ( singleDPDataArrayIndex >= singleDPData->dpSResult.size() )
        { return 0; }
    int limit = 4;
    while ( iter_singleDPData + limit < singleDPData->dpSResult[singleDPDataArrayIndex]->size() && ((*(singleDPData->dpSResult[singleDPDataArrayIndex]))[iter_singleDPData]).readID == ((*(singleDPData->dpSResult[singleDPDataArrayIndex]))[iter_singleDPData + limit]).readID )
        { ++iter_singleDPData; }
    if ( singleDPDataArrayIndex < singleDPData->dpSResult.size() )
    {
        memcpy( &singleAlgnmtResult, &((*(singleDPData->dpSResult[singleDPDataArrayIndex]))[iter_singleDPData]), sizeof( SingleAlgnmtResult ) );
        if ( lastFetchReadID == singleAlgnmtResult.readID )
            { ++lastFetchReadIDCount; }
        else
        {
            lastFetchReadID = singleAlgnmtResult.readID;
            lastFetchReadIDCount = 0;
        }
        iter_singleDPData++;
        ++counter;
        while ( singleDPDataArrayIndex < singleDPData->dpSResult.size() && iter_singleDPData >= singleDPData->dpSResult[singleDPDataArrayIndex]->size() )
        {
            ++singleDPDataArrayIndex;
            iter_singleDPData = 0;
        }
        return 1;
    }
    else // finished
    {
        return 0;
    }
}


// ****
HalfEndAlignmentEngine::HalfEndAlgnBatch::HalfEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, int inputMaxReadLength, uint * upkdReadLengths, int dPFromSingleEndDPResult
)
{
    MC_MemberCopy5 ( this->, , batchSize, peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low );
    MC_MemberCopy4 ( this->, , maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy4 ( this->, , index, queries, inputMaxReadLength, upkdReadLengths );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->isDoubleStrand        = ( peStrandLeftLeg == peStrandRightLeg );
    this->fullDNALength         = index->sraIndex->hsp->dnaLength;
    this->wordPerOldQuery       = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery          = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA            = MC_CeilDivide16 ( maxDNALength );
    if ( dPFromSingleEndDPResult )
    {
        MC_CheckMalloc ( singleEndDPCandidateInfo, SingleEndDPCandidateInfo, batchSize);
        canInfo = NULL;
    }
    else 
    {
        MC_CheckMalloc ( canInfo,             CandidateInfo,  batchSize );
        singleEndDPCandidateInfo = NULL;
    }
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( packedDNASequence,   uint,           batchSize * wordPerDNA );
    MC_CheckMalloc ( packedReadSequence,  uint,           batchSize * wordPerQuery );
    MC_CheckMalloc ( startLocs,           unsigned long long,           batchSize );
    MC_CheckMalloc ( hitLocs,             uint,           batchSize );
    MC_CheckMalloc ( scores,              int,            batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );
    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( peLeftAnchorLocs,    uint,           batchSize );
    MC_CheckMalloc ( peRightAnchorLocs,   uint,           batchSize );
    MC_CheckMalloc ( pattern,             uchar,          batchSize * patternLength );
    MC_CheckMalloc ( maxScoreCounts,      uint,           batchSize );
    clear();
}

HalfEndAlignmentEngine::HalfEndAlgnBatch::~HalfEndAlgnBatch()
{
    free ( canInfo );
    free ( singleEndDPCandidateInfo );
    free ( DNALengths );
    free ( packedDNASequence );
    free ( lengths );
    free ( packedReadSequence );
    free ( startLocs );
    free ( hitLocs );
    free ( scores );
    free ( cutoffThresholds );
    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( peLeftAnchorLocs );
    free ( peRightAnchorLocs );
    free ( pattern );
    free ( maxScoreCounts );
}

void HalfEndAlignmentEngine::HalfEndAlgnBatch::clear()
{
    numOfThreads = 0;
}

int HalfEndAlignmentEngine::HalfEndAlgnBatch::pack (
    SingleAlgnmtResult & singleAlgnmtResult
)
{
    uint alignedStrand = singleAlgnmtResult.strand;

    if ( alignedStrand != peStrandLeftLeg && alignedStrand != peStrandRightLeg )
    { return -2; }
    else
    {
        if ( numOfThreads + 1 + isDoubleStrand > batchSize )
        { return -1; }
    }

    uint alignedReadID = singleAlgnmtResult.readID;
    unsigned long long alignedPos = singleAlgnmtResult.algnmt;
    uint alignedReadLength = upkdReadLengths[alignedReadID];
    int  unalignedIsReadOrMate = 1 - ( alignedReadID & 1 );
    uint unalignedReadID = ( unalignedIsReadOrMate == 0 ?
                             alignedReadID - 1 : alignedReadID + 1 );
    uint unalignedReadLength = upkdReadLengths[unalignedReadID];
    
    #define MC_SetRead(strand) { \
            packRead(packedReadSequence, numOfThreads, \
                     unalignedReadID, unalignedReadLength, strand); \
            cutoffThresholds[numOfThreads] = std::max(DP_SCORE_THRESHOLD_RATIO * unalignedReadLength, DP_SCORE_THRESHOLD_LOWER_BOUND); \
            softClipLtSizes[numOfThreads] = (strand == 1) ? softClipLeft : softClipRight;  \
            softClipRtSizes[numOfThreads] = (strand == 1) ? softClipRight : softClipLeft; \
        }

    if ( peStrandLeftLeg == alignedStrand )
    {
        //aligned read: at left, unaligned read: at right
        unsigned long long rightEnd = alignedPos + insert_high;
        unsigned long long rightStart = alignedPos + insert_low - unalignedReadLength;
        // rightStart has to be >= alignedPos
        if (rightStart < alignedPos)
            rightStart = alignedPos;
        if ( rightStart < fullDNALength && rightEnd <= fullDNALength )
        {
            singleEndDPCandidateInfo[numOfThreads].refer = singleAlgnmtResult;
            singleEndDPCandidateInfo[numOfThreads].leftOrRight = 1;
            lengths[numOfThreads] = unalignedReadLength;
            startLocs[numOfThreads] = rightStart;
            DNALengths[numOfThreads] = rightEnd - rightStart;
            peLeftAnchorLocs[numOfThreads] = maxDNALength;
            peRightAnchorLocs[numOfThreads] = unalignedReadLength;
            repackDNA ( packedDNASequence, numOfThreads,
                        ( uint * ) index->sraIndex->hsp->packedDNA, rightStart, DNALengths[numOfThreads] );
            MC_SetRead ( peStrandRightLeg );
            ++numOfThreads;
        }
    }

    if ( peStrandRightLeg == alignedStrand )
    {
        //aligned read: at right, unaligned read: at left
        unsigned long long leftStart = alignedPos + alignedReadLength - insert_high;
        unsigned long long leftEnd = alignedPos + alignedReadLength - insert_low + unalignedReadLength;
        // leftEnd has to be < alignedPos + alignedReadLength
        if (leftEnd >= alignedPos + alignedReadLength)
            leftEnd = alignedPos + alignedReadLength - 1;

        if ( leftStart < fullDNALength && leftEnd <= fullDNALength )
        {
            singleEndDPCandidateInfo[numOfThreads].refer = singleAlgnmtResult;
            singleEndDPCandidateInfo[numOfThreads].leftOrRight = 0;
            lengths[numOfThreads] = unalignedReadLength;
            startLocs[numOfThreads] = leftStart;
            DNALengths[numOfThreads] = leftEnd - leftStart;
            peLeftAnchorLocs[numOfThreads] = insert_high - insert_low + 1;
            peRightAnchorLocs[numOfThreads] = 0;
            repackDNA ( packedDNASequence, numOfThreads,
                        ( uint * ) index->sraIndex->hsp->packedDNA, leftStart, DNALengths[numOfThreads] );
            MC_SetRead ( peStrandLeftLeg );
            ++numOfThreads;
        }
    }

    return numOfThreads;
}


inline void HalfEndAlignmentEngine::HalfEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}
inline void HalfEndAlignmentEngine::HalfEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, unsigned long long start, uint length
)
{
    size_t dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    uint *src = seq + start / 16;
    uint src_ofs = start % 16;
    uint *dst = packedSeq + dnaTPARA;
    uint dst_ofs = 1;
    *dst = 0;

    // we need a 0 at the beginning to reserve the first column of DP table
    while (length > 0) {
        uint len = min(min(16 - dst_ofs, 16 - src_ofs), length);
        *dst |= *src << (src_ofs * 2) >> (32 - len*2) << (32-(dst_ofs+len)*2);
        length -= len;
        dst_ofs += len;
        src_ofs += len;

        if (src_ofs == 16) { ++src; src_ofs = 0; }
        if (dst_ofs == 16) { dst += 32; dst_ofs = 0; *dst = 0; }
    }
}

// ****
void HalfEndAlignmentEngine::HalfEndAlgnThreadContext::init (
    HalfEndAlgnBatch * batch
)
{
    int batchSize = engine->DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                             engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
    sem_init ( &dispatchACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}

void HalfEndAlignmentEngine::HalfEndAlgnThreadContext::freeMemory()
{
    semiGlobalAligner.freeMemory ( );
    delete batch;
}

// ****
HalfEndAlignmentEngine::AlgnmtResultStream::AlgnmtResultStream()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}
HalfEndAlignmentEngine::AlgnmtResultStream::~AlgnmtResultStream()
{
}

// ****
void HalfEndAlignmentEngine::performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment )
{
    /* initialize */
    algnBatchCount = 0;
    dpAlignedRead = 0;
    dpAlignment = 0;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = insert_high - insert_low + inputMaxReadLength + 1;
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DP_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );
    algnSwapBatch =
        new HalfEndAlgnBatch ( DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                               peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                               maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                               index, queries, inputMaxReadLength, upkdReadLengths, engine->canStream->singleDPData != NULL);
    algnThreadContext = new HalfEndAlgnThreadContext[dpPara->numOfCPUThreads];

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        HalfEndAlgnBatch * batch =
            new HalfEndAlgnBatch ( DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                   peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                                   maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                   index, queries, inputMaxReadLength, upkdReadLengths, engine->canStream->singleDPData != NULL );
        // fprintf(stderr, "[%s] Batch size: %d\n", __func__, DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK);
        algnThreadContext[i].init ( batch );
    }

    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, algnmtCPUThread );

    numDPAlignedReads = ( uint * ) malloc ( dpPara->numOfCPUThreads * sizeof ( uint ) );
    memset ( numDPAlignedReads, 0, dpPara->numOfCPUThreads * sizeof ( uint ) );
    numDPAlignments = ( uint * ) malloc ( dpPara->numOfCPUThreads * sizeof ( uint ) );
    memset ( numDPAlignments, 0, dpPara->numOfCPUThreads * sizeof ( uint ) );

    /* perform alignment */
    int threadId;
    void * empty;

    if ( canStream->data != NULL )
    {
        assert(false);
    }
    else 
    {
        SingleAlgnmtResult singleAlgnmtResult;
        canStream->lastFetchReadID = -1;
        canStream->lastFetchReadIDCount = 0;
        canStream->counter = 0;

        SingleAlgnmtResult result;
        queue<SingleAlgnmtResult> results;
        int ret = canStream -> fetchNextSingleAlgnResult(result);
        while (ret) {
            results.push(result);
            uint readID = result.readID/2*2;
            while (ret = canStream->fetchNextSingleAlgnResult(result)) {
                if (result.readID/2*2==readID)
                    results.push(result);
                else 
                    break;
            }

            if (results.size() > algnSwapBatch->batchSize - algnSwapBatch->numOfThreads) {
                //launch one batch
                threadId = algnmtCPUThreadDelegator.schedule ( empty );
                sem_wait ( & ( algnThreadContext[threadId].dispatchACKSem ) );
                algnSwapBatch->clear();
            }
            while (!results.empty()) {
                SingleAlgnmtResult tmp = results.front();
                results.pop();
                inputFlags->set ((tmp.readID >> 1 ) << 1);
                algnSwapBatch->pack(tmp);
            }
        }
    }

    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].dispatchACKSem ) );
    }
    /* finalize */
    algnmtCPUThreadDelegator.finalize();
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;
    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory();
    }

    delete[] algnThreadContext;
    delete resultStream;
    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        numDPAlignedRead += numDPAlignedReads[i];
        numDPAlignment += numDPAlignments[i];
    }
    free ( numDPAlignedReads );
    free ( numDPAlignments );
}

void HalfEndAlignmentEngine::performAlignment (
    /* input */
    HalfEndOccStream   *  canStream,
    DPParameters     *    dpPara,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    DPInfoForReads * dpInfoForReads,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat,
    samfile_t * samOutputDPFilePtr, samfile_t ** currSamOutputFilePtr,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment
)
{
    engine = new HalfEndAlignmentEngine();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, queryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy ( engine->, , queryComments);
    MC_MemberCopy ( engine->, , dpInfoForReads );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy3 ( engine->, , accumReadNum, outputFormat, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    engine->currSamOutputFilePtr = currSamOutputFilePtr;
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    delete engine;
}

HalfEndAlignmentEngine * HalfEndAlignmentEngine::engine;
// ****
void DP_Space::algnmtCPUThread ( int threadId, void *& empty )
{
//printf("dp space cpu thread %d start\n",threadId);

    // Copy data, then ACK to dispatching thread
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++;
    HalfEndAlignmentEngine::HalfEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    sem_post ( & ( engine->algnThreadContext[threadId].dispatchACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    // launch kernel
    int * pThreadId = &threadId;
 //   engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
   // sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    //  timeRecorder.appendStart("GPUTime");
    engine->algnThreadContext[threadId].semiGlobalAligner.performAlignment ( batch->packedDNASequence, batch->DNALengths,
            batch->packedReadSequence, batch->lengths,
            batch->cutoffThresholds, batch->scores, batch->hitLocs,
            batch->maxScoreCounts,
            batch->pattern, batch->numOfThreads,
            batch->softClipLtSizes, batch->softClipRtSizes,
            batch->peLeftAnchorLocs, batch->peRightAnchorLocs );
// rearrange result and Output
    MC_MemberCopy2 ( int, engine->, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( uint *, batch->, hitLocs, DNALengths );
    MC_MemberCopy  ( unsigned long long *, batch->, startLocs );
    MC_MemberCopy2 ( uint *, batch->, peLeftAnchorLocs, peRightAnchorLocs );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    uchar * pattern = batch->pattern;
    CandidateInfo * canInfo =  batch->canInfo;
    SingleEndDPCandidateInfo * singleEndDPCandidateInfo = batch->singleEndDPCandidateInfo;
    DPResultBatch * resultBatch = &( engine->algnThreadContext[threadId].resultBatch );
    
    for ( int id = 0; id < batch->numOfThreads; id++ )
    {
        //Create record for AlgnmtDPResult;
        AlgnmtDPResult result;
        uint alignedID = (canInfo != NULL?canInfo[id].refer.readID:singleEndDPCandidateInfo[id].refer.readID);
        uint canInfoAmbPosition = (canInfo!=NULL?canInfo[id].refer.ambPosition : singleEndDPCandidateInfo[id].refer.algnmt );
        char canStrand = ( canInfo!=NULL?canInfo[id].refer.strand : singleEndDPCandidateInfo[id].refer.strand);
        uint canLeftOrRight = (canInfo!=NULL?canInfo[id].leftOrRight: singleEndDPCandidateInfo[id].leftOrRight );
        int canScore = ( canInfo!=NULL? canInfo[id].refer.mismatchCount: singleEndDPCandidateInfo[id].refer.score );
        
        uint alignedIsReadOrMate = alignedID & 1;
        result.readID = alignedID - alignedIsReadOrMate;
        unsigned long long dpAlgnmtPos;
        if ( batch->scores[id] >= batch->cutoffThresholds[id])
        {
            CigarStringEncoder<void> encoder;
            uchar lastType = 'N';

            for ( uchar * p = pattern + id * engine->patternLength; *p != 0; p++ )
            {
                if ( *p == 'V' )
                {
                    encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                }
                else
                {
                    encoder.append ( *p, 1 );
                    lastType = *p;
                }
            }

            encoder.encodeCigarString ( openGapScore, extendGapScore );
            result.cigarString = encoder.cigarString;
            // To get edit distance
            int L = batch->lengths[id] - encoder.charCount['I'] - encoder.charCount['S'];
            int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[id] ) /
                                ( matchScore - mismatchScore );
            result.editdist = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
            result.whichFromDP = 1 - alignedIsReadOrMate;
            dpAlgnmtPos = startLocs[id] + hitLocs[id];
            
            
            if ( dpAlgnmtPos < canInfoAmbPosition  )
            {
                // dp is on left
                result.insertSize = canInfoAmbPosition - dpAlgnmtPos +
                                    engine->upkdReadLengths[alignedID];
            }
            else
            {
                // dp is on right
                result.insertSize = dpAlgnmtPos - canInfoAmbPosition +
                                    batch->lengths[id] + encoder.charCount['D'] -
                                    encoder.charCount['I'] - encoder.charCount['S'];
            }
            result.startPos = startLocs[id];
            result.refDpLength = DNALengths[id];
            result.peLeftAnchor = peLeftAnchorLocs[id];
            result.peRightAnchor = peRightAnchorLocs[id];
            result.num_sameScore = batch->maxScoreCounts[id]; //TODO
        }
        else
        {
            result.cigarString = NULL;
            result.whichFromDP = 2;
            dpAlgnmtPos = kMaxULL;
        }
        if ( alignedIsReadOrMate == 0 )
        {
            // aligned is read, unaligned is mate
            result.algnmt_1 = canInfoAmbPosition;
            result.algnmt_2 = dpAlgnmtPos;
            result.score_1 = canScore;
            result.score_2 = batch->scores[id];
            result.strand_1 = canStrand;
            result.strand_2 = ( canLeftOrRight == 0 ? peStrandLeftLeg : peStrandRightLeg );
            if ( canInfo == NULL )
            {
               result.result_1 = ( SingleAlgnmtResult * ) malloc ( sizeof ( SingleAlgnmtResult ) );
               memcpy ( result.result_1, &(singleEndDPCandidateInfo[id].refer), sizeof ( SingleAlgnmtResult ) );
               // result.result_1 = &(singleEndDPCandidateInfo[id].refer);
            }
            else
                { result.result_1 = NULL; }
            result.result_2 = NULL;
        }
        else
        {
            // aligned is mate, unaligned is read
            result.algnmt_1 = dpAlgnmtPos;
            result.algnmt_2 = canInfoAmbPosition;
            result.score_1 = batch->scores[id];
            result.score_2 = canScore;
            result.strand_1 = ( canLeftOrRight == 0 ? peStrandLeftLeg : peStrandRightLeg );
            result.strand_2 = canStrand;
            result.result_1 = NULL;
            if ( canInfo == NULL )
            {
               result.result_2 = ( SingleAlgnmtResult * ) malloc ( sizeof ( SingleAlgnmtResult ) );
               memcpy ( result.result_2, &(singleEndDPCandidateInfo[id].refer), sizeof ( SingleAlgnmtResult ) );
               // result.result_2 = &(singleEndDPCandidateInfo[id].refer);
            }
            else
                { result.result_2 = NULL; }
        }
        resultBatch->push_back ( result );
    }
    // output thread
    int *pid = &threadId;
    DPOutputThread( threadId, pid );
//printf("dp space cpu thread %d end\n",threadId);

}

void DP_Space::DPOutputThread ( int outputThreadId, int *& pCallThreadId )
{
    int threadId = *pCallThreadId;
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    DPResultBatch * resultBatch = &(engine->algnThreadContext[threadId].resultBatch);
    DeepDPAlignResult * tmpBestResult = NULL;

#define MC_OutputRead() { \
        outputBuf->ready(1); \
        DeepDPAlignResult * alignResult = (DeepDPAlignResult *)malloc(outputBuf->size*sizeof(DeepDPAlignResult)); \
        AlgnmtDPResult * tmpDPPairResult = (AlgnmtDPResult *) outputBuf->elements; \
        if (outputBuf->size > 0 ) { \
            if ( tmpDPPairResult[0].result_1 == NULL && tmpDPPairResult[0].result_2 == NULL) { \
                assert(false); \
            } \
            else { \
                int validResult = 0; \
                for (int k=0,l=0;k<outputBuf->size;++k) { \
                    if ( outputBuf->elements[0].whichFromDP < 2 ) \
                    { \
                        alignResult[l].readID = tmpDPPairResult[k].readID; \
                        alignResult[l].insertSize = tmpDPPairResult[k].insertSize; \
                        if ( tmpDPPairResult[k].result_1 != NULL ) { \
                            alignResult[l].algnmt_1 = tmpDPPairResult[k].result_1->algnmt; \
                            alignResult[l].strand_1 = tmpDPPairResult[k].result_1->strand; \
                            alignResult[l].score_1 = tmpDPPairResult[k].result_1->score; \
                            alignResult[l].editdist_1 = tmpDPPairResult[k].result_1->editdist; \
                            alignResult[l].cigarString_1 = tmpDPPairResult[k].result_1->cigarString; \
                            alignResult[l].num_sameScore_1 = tmpDPPairResult[k].result_1->num_sameScore; \
                            alignResult[l].startPos_1 = tmpDPPairResult[k].result_1->startPos; \
                            alignResult[l].refDpLength_1 = tmpDPPairResult[k].result_1->refDpLength; \
                            alignResult[l].peLeftAnchor_1 = tmpDPPairResult[k].result_1->peLeftAnchor; \
                            alignResult[l].peRightAnchor_1 = tmpDPPairResult[k].result_1->peRightAnchor; \
                            alignResult[l].algnmt_2 = tmpDPPairResult[k].algnmt_2; \
                            alignResult[l].strand_2 = tmpDPPairResult[k].strand_2; \
                            alignResult[l].score_2 = tmpDPPairResult[k].score_2; \
                            alignResult[l].editdist_2 = tmpDPPairResult[k].editdist; \
                            alignResult[l].cigarString_2 = tmpDPPairResult[k].cigarString; \
                            alignResult[l].num_sameScore_2 = tmpDPPairResult[k].num_sameScore; \
                            alignResult[l].startPos_2 = tmpDPPairResult[k].startPos; \
                            alignResult[l].refDpLength_2 = tmpDPPairResult[k].refDpLength; \
                            alignResult[l].peLeftAnchor_2 = tmpDPPairResult[k].peLeftAnchor; \
                            alignResult[l].peRightAnchor_2 = tmpDPPairResult[k].peRightAnchor; \
                        } \
                        else if ( tmpDPPairResult[k].result_2 != NULL ) { \
                            alignResult[l].algnmt_1 = tmpDPPairResult[k].algnmt_1; \
                            alignResult[l].strand_1 = tmpDPPairResult[k].strand_1; \
                            alignResult[l].score_1 = tmpDPPairResult[k].score_1; \
                            alignResult[l].editdist_1 = tmpDPPairResult[k].editdist; \
                            alignResult[l].cigarString_1 = tmpDPPairResult[k].cigarString; \
                            alignResult[l].num_sameScore_1 = tmpDPPairResult[k].num_sameScore; \
                            alignResult[l].startPos_1 = tmpDPPairResult[k].startPos; \
                            alignResult[l].refDpLength_1 = tmpDPPairResult[k].refDpLength; \
                            alignResult[l].peLeftAnchor_1 = tmpDPPairResult[k].peLeftAnchor; \
                            alignResult[l].peRightAnchor_1 = tmpDPPairResult[k].peRightAnchor; \
                            alignResult[l].algnmt_2 = tmpDPPairResult[k].result_2->algnmt; \
                            alignResult[l].strand_2 = tmpDPPairResult[k].result_2->strand; \
                            alignResult[l].score_2 = tmpDPPairResult[k].result_2->score; \
                            alignResult[l].editdist_2 = tmpDPPairResult[k].result_2->editdist; \
                            alignResult[l].cigarString_2 = tmpDPPairResult[k].result_2->cigarString; \
                            alignResult[l].num_sameScore_2 = tmpDPPairResult[k].result_2->num_sameScore; \
                            alignResult[l].startPos_2 = tmpDPPairResult[k].result_2->startPos; \
                            alignResult[l].refDpLength_2 = tmpDPPairResult[k].result_2->refDpLength; \
                            alignResult[l].peLeftAnchor_2 = tmpDPPairResult[k].result_2->peLeftAnchor; \
                            alignResult[l].peRightAnchor_2 = tmpDPPairResult[k].result_2->peRightAnchor; \
                        } \
                        ++l; \
                    } \
                    validResult = l; \
                } \
                if ( validResult ) \
                { \
                    tmpBestResult = outputDeepDPResult2(alignResult, outputBuf->size, \
                                    engine->queries, engine->upkdReadLengths, \
                                    engine->origReadIDs, engine->queryNames, engine->queryComments, engine->upkdQualities, \
                                    engine->inputMaxReadLength, engine->accumReadNum, engine->outputFormat, \
                                    engine->currSamOutputFilePtr[threadId], engine->index, \
                                    engine->peStrandLeftLeg, engine->peStrandRightLeg, NULL); \
                    engine->dpInfoForReads->startPositions[tmpBestResult->readID] = tmpBestResult->startPos_1; \
                    engine->dpInfoForReads->strand_dpLengths[tmpBestResult->readID] = ( ( ( tmpBestResult->strand_1 - 1 ) << 15 ) | ( tmpBestResult->refDpLength_1 ) ); \
                    engine->dpInfoForReads->peLeftAnchors[tmpBestResult->readID] = tmpBestResult->peLeftAnchor_1; \
                    engine->dpInfoForReads->peRightAnchors[tmpBestResult->readID] = tmpBestResult->peRightAnchor_1; \
                    engine->dpInfoForReads->startPositions[tmpBestResult->readID + 1] = tmpBestResult->startPos_2; \
                    engine->dpInfoForReads->strand_dpLengths[tmpBestResult->readID + 1] = ( ( ( tmpBestResult->strand_2 - 1 ) << 15 ) | ( tmpBestResult->refDpLength_2 ) ); \
                    engine->dpInfoForReads->peLeftAnchors[tmpBestResult->readID+1] = tmpBestResult->peLeftAnchor_2; \
                    engine->dpInfoForReads->peRightAnchors[tmpBestResult->readID+1] = tmpBestResult->peRightAnchor_2; \
                } \
            } \
            engine->numDPAlignedReads[threadId] += 1; \
            engine->numDPAlignments[threadId] += outputBuf->size; \
            engine->alignFlags->set(lastReadID); \
        } \
        free (alignResult) ; \
    }

    OutputBuffer<AlgnmtDPResult> *outputBuf;
    outputBuf = new OutputBuffer<AlgnmtDPResult>();
    outputBuf->setAlignmentType ( engine->alignmentType );
    DPResultBatch & batch = *resultBatch;
    int lastReadID = -1;
    for ( int i = 0; i < batch.size(); i++ )
    {
        AlgnmtDPResult & result = batch[i];
        
        if ( result.readID != lastReadID )
        {
            MC_OutputRead();
            outputBuf->clear();
            lastReadID = result.readID;
        }

        outputBuf->add ( result );
    }
    MC_OutputRead();
    outputBuf->clear();

    delete outputBuf;
    for ( int i=0;i<batch.size();++i )
    {
        free ( batch[i].cigarString );
    }
    batch.clear();
}

void DP_Space::DPOutputThreadFinalize()
{
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// deep-dp space ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
using namespace DeepDP_Space;
#define DP2_SEEDING_BATCH_SIZE 4 * 1024 * 1024
#define DP2_MARGIN(l) ((l>100) ? 30 : 25)

// ****
DeepDP_Space::CandidateStream::CandidateStream()
{
    pthread_mutex_init ( &occupy_mutex, NULL );
}
void DeepDP_Space::CandidateStream::append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags )
{
    pthread_mutex_lock ( &occupy_mutex );

    for ( vector<CandidateInfo>::iterator it = canInfo->begin();
            it < canInfo->end(); ++it )
    {
        data.push_back ( *it );
        alignFlags->set ( ( it->readIDLeft >> 1 ) << 1 );
    }

    pthread_mutex_unlock ( &occupy_mutex );
}

// ****
PairEndSeedingEngine::PairEndSeedingBatch::PairEndSeedingBatch (
    uint batchSize, DPParameters * dpPara,
    uint * queries, uint * readLengths, uint inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg, BWT * bwt, PairEndSeedingEngine * home
)
{
    numberOfShortcut = 0;
    lazyCnt = 0;
    this->home = home;
    MC_MemberCopy4 ( this->, , batchSize, queries, readLengths, bwt );
    MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    this->numOfCPUForSeeding = dpPara->numOfCPUForSeeding;

    for ( int i = 0; i < 2; i++ )
    {
        this->maxHitNum[i]  = dpPara->paramRead[i].maxHitNum;
    }

    this->maxSeedLength = inputMaxReadLength;
    this->wordPerSeed   = MC_CeilDivide16 ( maxSeedLength );
    this->wordPerQuery  = getWordPerQuery ( inputMaxReadLength );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        MC_CheckMalloc ( readIDs[lOr],    uint,   batchSize );
        MC_CheckMalloc ( lengths[lOr],    uint,   batchSize );
        MC_CheckMalloc ( minLengths[lOr], uint,   batchSize );
        MC_CheckMalloc ( offsets[lOr],    uint,   batchSize );
        MC_CheckMalloc ( seeds[lOr],      uint,   batchSize * wordPerSeed );
    }

    this->numNBMismatch = dpPara->seedProperties.numNBMismatch;
    this->restrictedNBM = dpPara->seedProperties.restrictedNBM ;

    MC_CheckMalloc ( seedPositions,       int,    inputMaxReadLength );
    clear();
}

PairEndSeedingEngine::PairEndSeedingBatch::~PairEndSeedingBatch()
{
    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        free ( readIDs[lOr] );
        free ( lengths[lOr] );
        free ( minLengths[lOr] );
        free ( offsets[lOr] );
        free ( seeds[lOr] );
    }
    free ( seedPositions );
}

void PairEndSeedingEngine::PairEndSeedingBatch::clear()
{
    for ( int i = 0; i < 2; i++ )
    {
        numQueries[i] = 0;
        inPosArr[i].clear();
    }

    lastPairID = -1;
}
uint PairEndSeedingEngine::PairEndSeedingBatch::findRevStart (
    SeedPos * arr, uint len
)
{
    if ( len == 0 || ! ( arr[len - 1].strand_readID >> 31 ) )
    { return len; }

    uint start = 0;
    uint end = len - 1;

    while ( start < end )
    {
        uint mid = ( start + end ) / 2;

        if ( ( arr[mid].strand_readID >> 31 ) )
        {
            // reverse
            end = mid;
        }
        else
        {
            // forward
            start = mid + 1;
        }
    }

    return start;
}

inline void PairEndSeedingEngine::PairEndSeedingBatch::pack (
    uint evenReadID, uint readID, int off, int seedLength, int varMinLength, int readOrMate
)
{   
    int seedID = numQueries[readOrMate];
    readIDs[readOrMate][seedID] = evenReadID;
    lengths[readOrMate][seedID] = seedLength;
    minLengths[readOrMate][seedID] = varMinLength;
    offsets[readOrMate][seedID] = off;
#define MC_OldReadUnpackIn(X,i) ((X[oldReadTPARA + ((i>>4)<<5)] >> ((i & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerQuery + ( readID % 32 );
    uint seedTPARA = ( seedID / 32 ) * 32 * wordPerSeed + ( seedID % 32 );

    for ( int i = 0; i < wordPerSeed; i++ )
    {
        seeds[readOrMate][seedTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 0; i < seedLength; i++ )
    {
        int pos = off + i;
        seeds[readOrMate][seedTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) MC_OldReadUnpackIn ( queries, pos ) << ( ( i & 0xF ) << 1 ) ;
    }

    ++numQueries[readOrMate];
}

inline void PairEndSeedingEngine::PairEndSeedingBatch::packForMmpSeeding ( QueryIDStream * queryIDStream, int readOrMate )
{
    vector<int> & evenReadIDs = *(queryIDStream->data);
    for ( int i=0;i<evenReadIDs.size(); ++i )
    {
        readIDs[readOrMate][i] = evenReadIDs[i];
        home->inputFlags->set( evenReadIDs[i] );
    }
    numQueries[readOrMate] = evenReadIDs.size();
}

int PairEndSeedingEngine::PairEndSeedingBatch::packSeeds (
    uint evenReadID, int stage
)
{
    for ( int i = 0; i < 2; i++ )
    {
        uint readID = evenReadID + i;
        uint readLength = readLengths[readID];
        int seedNum, seedLength;
        // getSeedPositions ( STAGE_DEEP_DP_ROUND1, readLength, &seedLength, seedPositions, &seedNum );
        SeedingProperties * seedingProperties = &(home->dpPara->seedProperties);
        seedNum = seedingProperties->numberOfSeed;
        if ( numQueries[i] + seedNum > batchSize )
        {
            return 0;
        }
  
        for ( int j = 0; j < seedNum; j++ )
        {
            if ( seedingProperties->seedPos[j] + seedingProperties->maxSeedLength[j] <= readLength )
            {   
                // normal seeding
                pack (evenReadID, readID, seedingProperties->seedPos[j], seedingProperties->maxSeedLength[j], seedingProperties->minSeedLength[j], i); 
                lastPairID = evenReadID >> 1;
            }
            // pack ( evenReadID, readID, seedPositions[j], seedLength, HALF_SEED_SIZE_FOR_V_LONG_READ_LB<<1, i );
        }

    }
    
    return 1;
}

int PairEndSeedingEngine::PairEndSeedingBatch::packSeeds (
    SRAOccurrence & occ, int stage
)
{
    uint readID = occ.readID;
    uint pairID = readID >> 1;
    int readOrMate = readID & 1;
    uint evenReadID = readID - readOrMate;

    if ( pairID != lastPairID )
    {
        if ( !packSeeds ( evenReadID, stage ) )
        {
            return 0;
        }
    }

    SeedPos pos;
    pos.pos = occ.ambPosition;
    pos.strand_readID = ( ( occ.strand - 1 ) << 31 ) | evenReadID;
    inPosArr[readOrMate].push_back ( pos );
    return 1;
}

void PairEndSeedingEngine::PairEndSeedingBatch::pairEndMerge (
    vector<CandidateInfo> * &pairEndPos,
    SeedPos * readPos, SeedPos * matePos,
    int isMatePositive
)
{
#define MC_DecodePos(x) ((x)->pos)
#define MC_DecodeID(x) ((x)->strand_readID & 0x7FFFFFFF)
#define MC_ReadID() MC_DecodeID(readIter)
#define MC_MateID() MC_DecodeID(mateIter)
    SeedPos * readIter = readPos;
    SeedPos * mateIter = matePos;
    uint readID;
    uint mateID;
    
    uint roundNum = engine->dpPara->roundNum;
    while ( true )
    {
        mateID = MC_MateID();

        while ( MC_ReadID() < mateID )
        { ++readIter; }

        readID = MC_ReadID();

        while ( MC_MateID() < readID )
        { ++mateIter; }

        mateID = MC_MateID();

        if ( mateID == 0x7FFFFFFF )
        { break; }
        else if ( readID < mateID )
        { continue; }

        SeedPos * readStart = readIter;
        SeedPos * mateStart = mateIter;

        // assert : readID == mateID
        while ( MC_ReadID() == readID )
        { ++readIter; }

        while ( MC_MateID() == mateID )
        { ++mateIter; }

        SeedPos * readEnd = readIter;
        SeedPos * mateEnd = mateIter;
#define MC_Compress(start, end, divideGap) { \
            SeedPos *cmprReadIter = start; \
            register unsigned long long prevLoc = MC_DecodePos(cmprReadIter); \
            for (SeedPos* p = start+1; p < end; p++) { \
                register unsigned long long curLoc = MC_DecodePos(p); \
                if (prevLoc + divideGap < curLoc) { \
                    *(++cmprReadIter) = *p; \
                    prevLoc = curLoc; \
                } \
            } \
            end = cmprReadIter + 1; \
        }
        MC_Compress ( readStart, readEnd, DP2_DIVIDE_GAP );
        // get the read length of negative strand !
        int readLength = readLengths[readID /2*2 + 1 - isMatePositive];
        //fprintf(stderr,"my readlen %d wrong readlen !! %d\n", readLength, readLengths[readID]);
        int margin = DP2_MARGIN ( readLength );
        int length_low = insert_low - readLength - margin;
        if (length_low < 0 )
            length_low = 0;
        int length_high = insert_high - readLength + margin;
        SeedPos * readP = readStart;
        SeedPos * mateP = mateStart;
        register unsigned long long readLoc = MC_DecodePos ( readP );
        register unsigned long long mateLoc = MC_DecodePos ( mateP );

        SeedPos * matePreStart = mateP;
        while ( readP < readEnd )
        {
            bool updatePreStart = false;
            readLoc = MC_DecodePos ( readP );
            for ( SeedPos * matei = matePreStart; matei < mateEnd; ++matei )
            {
                mateLoc = MC_DecodePos ( matei );
                if ( readLoc + length_high < mateLoc )
                {
                    break;
                }
                else if ( readLoc + length_low <= mateLoc )
                {
                    CandidateInfo ci;
                    ci.pos[0] = readLoc;
                    ci.pos[1] = mateLoc;
                    ci.readIDLeft = MC_DecodeID ( readP ) + isMatePositive;
                    pairEndPos->push_back ( ci );
                    if ( !updatePreStart )
                    {
                        matePreStart = matei;
                    }
                }    
            }
            ++readP;
            readLoc = MC_DecodePos( readP );
        }
    }
}

struct RadixTraitsCandInfo {
    static const int nBytes = 4;
    int kth_byte(const DeepDP_Space::CandidateInfo &x, int k) {
        return x.readIDLeft >> (k << 3) & 0xFF;
    }
    bool compare(const DeepDP_Space::CandidateInfo &x, const DeepDP_Space::CandidateInfo &y) {
        return x.readIDLeft < y.readIDLeft;
    }
};

struct CompareCandInfo {
    bool operator()(const DeepDP_Space::CandidateInfo &x, const DeepDP_Space::CandidateInfo &y) {
        return x.readIDLeft < y.readIDLeft;
    }
};

vector<DeepDP_Space::CandidateInfo> * PairEndSeedingEngine::PairEndSeedingBatch::mergeAndPairPairedEnd ( SeedPos * readPos, SeedPos * matePos, int readPosLen, int matePosLen )
{
    SeedPos * readArr[2], *mateArr[2];
    // 0 -- forward, 1 -- reverse
    readArr[0] = readPos;
    mateArr[0] = matePos;
    readArr[1] = readPos + findRevStart ( readPos, readPosLen );
    mateArr[1] = matePos + findRevStart ( matePos, matePosLen );
    vector<CandidateInfo> * canInfo = new vector<CandidateInfo>;
    // pair positvie read and negative mate
    pairEndMerge ( canInfo, readArr[peStrandLeftLeg - 1], mateArr[peStrandRightLeg - 1], 0 );
    // pair positive mate and negative read
    pairEndMerge ( canInfo, mateArr[peStrandLeftLeg - 1], readArr[peStrandRightLeg - 1], 1 );
    free ( readPos );
    free ( matePos );
    // for ( int i=0;i<canInfo->size();++i )
    //     { fprintf( stderr, "%u %u\n", (*canInfo)[i].pos[0], (*canInfo)[i].pos[1] ); }
    // Sort the candidates so that readID will be in order
    // To be revised
    vector<CandidateInfo> & candArr = *canInfo;
    uint arrLength = candArr.size();
    double start_time = setStartTime();
    // kx::radix_sort(candArr.begin(), candArr.end(), RadixTraitsCandInfo());
    __gnu_parallel::stable_sort(candArr.begin(), candArr.end(), CompareCandInfo());
    // CandidateInfo * auxCandArr;
    // stable sort looks more faster after sorting. so do not use kxsort here
    // MC_CheckMalloc ( auxCandArr, CandidateInfo, arrLength );
    // MC_RadixSort_32_16 ( candArr, readIDLeft, auxCandArr, arrLength );
    // free ( auxCandArr );
    fprintf(stderr, "[%s] Sorting n & time: %u, %f\n", __func__, candArr.size(), getElapsedTime(start_time));
    return canInfo;    
}

// ****
void PairEndSeedingEngine::PairEndSeedingThreadContext::init (
    PairEndSeedingBatch * batch
)
{
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    this->batch = batch;
    this->batch->clear();
}
void PairEndSeedingEngine::PairEndSeedingThreadContext::freeMemory()
{
    delete batch;
}

// ****
PairEndSeedingEngine::PairEndSeedingEngine()
{
    queryIDStream = NULL;
    halfEndOccStream = NULL;
    tooManyHitIDStream = NULL;
}

struct SeedAlign {
  unsigned long long offset;
  uint multiplicity : 16;
  int length : 16;
  uint query_offset : 16;
  SeedAlign() {}
  SeedAlign(unsigned long long o, uint m, int l, uint qo) {
    offset = o;
    multiplicity = m;
    length = l;
    query_offset = qo;
  }

  bool operator< (const SeedAlign &rhs) const {
    return offset < rhs.offset;
  }
};

struct SeedSAalign {
    uint query_offset : 10;
    unsigned long long sa_l;
    uint sa_diff : 10;
    uint seed_len : 12;

    SeedSAalign (uint qo = 0, unsigned long long sa_l = 0, uint sa_diff = 0, uint seed_len = 0)
        : query_offset(qo), sa_l(sa_l), sa_diff(sa_diff), seed_len(seed_len) {}

    bool valid() { return seed_len > 0; }
};

#define QUERY_BLOCK 32

template <int T>
void mmp(uint*, vector<SeedSAalign>&, BWT*, LT *lkt, int, int, unsigned long long, MmpProperties);

#define DEBUG_RESEED do { \
    string s; \
    for (int i = 0; i < len; ++i) {\
        char c = (query[i/CHAR_PER_WORD] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK; \
        s.push_back("ACGT"[c]); \
    } \
    fprintf(stderr, "[%s] %s Reseeding: Len %d->%d, range %d->%d\n", __func__, s.c_str(), seed_len, last_seed_len, r-l+1, last_r-last_l+1); \
} while (0)

#define CHECK_AND_SET_LAST \
    do { \
        if (seed_len >= mmpPara.seedMinLength && nextr - nextl < r - l) { \
            last_r = r; \
            last_l = l; \
            last_seed_len = seed_len; \
        } \
    } while (0)

#define CHECK_AND_ADD_RANGE(_x_)                                          \
      do {                                                              \
        int diff = 0;                                                    \
        int x = _x_;                                                   \
        if (seed_len >= mmpPara.seedMinLength) {                             \
            if (seed_len >= mmpPara.reseedLen && last_r - last_l + 1 <= mmpPara.seedSAsizeThreshold && \
                (seed_len - last_seed_len <= mmpPara.reseedAbsDiff || seed_len * mmpPara.reseedRLTratio < last_seed_len)) { \
                diff = seed_len - last_seed_len; \
                l = last_l, r = last_r, seed_len = last_seed_len;\
            } \
            result.emplace_back(SeedSAalign(x, l, min(1ULL * mmpPara.seedSAsizeThreshold, r - l), seed_len)); \
            if (false) fprintf(stderr, "[%s] i: %d, x: %d, sl: %d, l: %llu, r: %llu\n", __func__, i, x, seed_len, l, r); \
        } /* else discard this range */                                 \
        i -= diff; \
        i -= min(seed_len, mmpPara.seedMinLength); \
        if (false) fprintf(stderr, "[%s] i: %d, sl: %d, r-l: %d, diff: %d\n", __func__, i, seed_len, r-l, diff); \
        l = 0; \
        r = text_length; \
        revl = 0; \
        revr = text_length; \
        seed_len = 0; \
        last_l = l, last_r = r, last_seed_len = 0;  \
    } while (0)

#define BIT_PER_CHAR LOOKUP_BIT_PER_CHAR
#define CHAR_MASK LOOKUP_CHAR_MASK

template<>
//positive strand backward
void mmp<0>(uint* query, vector<SeedSAalign> &result, BWT* bwt, LT *lkt, int start, int len, unsigned long long text_length, MmpProperties mmpPara) {
    unsigned char c;
    int i = 0, seed_len = 0, n_ranges = 0;
    unsigned long long l = 0, r = text_length, nextl, nextr, revl, revr;
    unsigned long long last_l = l, last_r = r, last_seed_len = 0;

    vector<int> lkp(start + len);
    for (int i = start, p = 0; i < start + len; ++i) {
        int j = start+len-(i-start)-1;
        int c = (query[j/CHAR_PER_WORD] >> ( j % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
        p = (p >> LOOKUP_BIT_PER_CHAR) | (c << LOOKUP_BIT_PER_CHAR * (LOOKUP_SIZE-1));
        if (i >= LOOKUP_SIZE - 1) {
            lkp[i-LOOKUP_SIZE+1] = p;
        }
    }

    for (i=start; i<start+len ; ++i) {
        if (seed_len == 0) {
            if (start + len - i < mmpPara.seedMinLength) { break; }
            nextl = lkp[i] == 0 ? 1 : (lkt->table[lkp[i]-1]+1);
            nextr = lkt->table[lkp[i]];

            i += LOOKUP_SIZE - 1;
            seed_len = LOOKUP_SIZE - 1;
        } else {
            int j = start+len-(i-start)-1;
            c = (query[j/CHAR_PER_WORD] >> ( j % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
            nextl = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
            nextr = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r+1, c);
        }

        if (nextl <= nextr) {
            CHECK_AND_SET_LAST;
            l = nextl;
            r = nextr;
            ++seed_len;
        } else {
            CHECK_AND_ADD_RANGE(start+len-(i-start));
        }
    }
    CHECK_AND_ADD_RANGE(start);
}

template<> 
//positive strand forward
void mmp<1>(uint* query, vector<SeedSAalign> &result, BWT* bwt, LT *lkt, int start, int len, unsigned long long text_length, MmpProperties mmpPara) {
    unsigned char c;
    int i = 0, seed_len = 0, n_ranges = 0;
    unsigned long long l = 0, r = text_length, revl = 0, revr = text_length,  nextl, nextr, oL[4], oR[4], oCount[4];
    unsigned long long last_l = l, last_r = r, last_seed_len = 0;

    vector<int> lkp(start + len);
    for (int i = start, p = 0; i < start + len; ++i) {
        int c = (query[i/CHAR_PER_WORD] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
        p = (p >> LOOKUP_BIT_PER_CHAR) | (c << LOOKUP_BIT_PER_CHAR * (LOOKUP_SIZE-1));
        if (i >= LOOKUP_SIZE - 1) {
            lkp[i-LOOKUP_SIZE+1] = p;
        }
    }

    for (i=start; i<start+len; ++i) {
        if (seed_len == 0) {
            if (start + len - i < mmpPara.seedMinLength) { break; }
            nextl = lkp[i] == 0 ? 1 : (lkt->table[lkp[i]-1]+1);
            nextr = lkt->table[lkp[i]];
            i += LOOKUP_SIZE - 1;
            seed_len = LOOKUP_SIZE - 1;
        } else {
            c = (query[i/CHAR_PER_WORD] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
            BWTAllOccValue(bwt,revl,oL);
            BWTAllOccValue(bwt,revr+1,oR);
        	oCount[3] = 0;
        	for (int k=2;k>=0;k--) oCount[k] = oCount[k+1]+oR[k+1] - oL[k+1];
            nextl = bwt->cumulativeFreq[c] + oL[c] + 1;
            nextr = bwt->cumulativeFreq[c] + oR[c];
        }

        if (nextl <= nextr) {
            CHECK_AND_SET_LAST;
            revl = nextl;
            revr = nextr;
    	    r = r - oCount[c];
    	    l = r - (revr-revl);
            ++seed_len;
        } else {
            CHECK_AND_ADD_RANGE(i-seed_len);
        }
    }
    CHECK_AND_ADD_RANGE(start+len-seed_len);
}

template<> 
//negative strand backward
void mmp<2>(uint* query, vector<SeedSAalign> &result, BWT* bwt, LT *lkt, int start, int len, unsigned long long text_length, MmpProperties mmpPara) {
    unsigned char c;
    int i = 0, seed_len = 0, n_ranges = 0;
    unsigned long long l = 0, r = text_length, nextl, nextr, revl, revr;
    unsigned long long last_l = l, last_r = r, last_seed_len = 0;

    vector<int> lkp(start + len);
    for (int i = start, p = 0; i < start + len; ++i) {
        int c = (query[i/CHAR_PER_WORD] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
        p = (p >> LOOKUP_BIT_PER_CHAR) | ((3-c) << LOOKUP_BIT_PER_CHAR * (LOOKUP_SIZE-1));
        if (i >= LOOKUP_SIZE - 1) {
            lkp[i-LOOKUP_SIZE+1] = p;
        }
    }

    for (i=start; i<start+len ; ++i) {
        if (seed_len == 0) {
            if (start + len - i < mmpPara.seedMinLength) { break; }
            nextl = lkp[i] == 0 ? 1 : (lkt->table[lkp[i]-1]+1);
            nextr = lkt->table[lkp[i]];

            // validate correctness of LKT
            // unsigned long long l = 0, r = text_length;
            // for (int k = 0; k < LOOKUP_SIZE; ++k) {
            //     char a = (lkp[i] >> LOOKUP_BIT_PER_CHAR * k) & CHAR_MASK;
            //     char b = (query[(i+k)/CHAR_PER_WORD] >> ( (i+k) % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
            //     b = 3 - b;
            //     assert(a == b);
            //     l = bwt->cumulativeFreq[a] + BWTOccValue(bwt, l, a) + 1;
            //     r = bwt->cumulativeFreq[a] + BWTOccValue(bwt, r+1, a);
            // }
            // fprintf(stderr, "lktL: %llu, lktR: %llu, bwtL: %llu, bwtR: %llu\n", nextl, nextr, l, r);

            // for (unsigned k = 0; k < 1 << LOOKUP_BIT_PER_CHAR * LOOKUP_SIZE; ++k) {
            //     if (lkt->table[k] == r) {
            //         fprintf(stderr, "%x  <-->  %x\n", k, lkp[i]);
            //     }
            // }

            i += LOOKUP_SIZE - 1;
            seed_len = LOOKUP_SIZE - 1;
        } else {
            c = (query[i/CHAR_PER_WORD] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR)) & CHAR_MASK;
    	    c = 3-c;
            nextl = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
            nextr = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r+1, c);
        }

        if (nextl <= nextr) {
            CHECK_AND_SET_LAST;
            l = nextl;
            r = nextr;
            ++seed_len;
        } else {
            CHECK_AND_ADD_RANGE(i-seed_len);
        }
    }
    CHECK_AND_ADD_RANGE(start+len-seed_len);
}

typedef struct MMP_ARG {
    BWT* lbwt;
    BWT* rbwt;
    LT* lkt;
    LT* rlkt;
    unsigned long long text_len;
    uint* read_lens;
    uint wpq;
    uint* queries;
    vector<vector<vector<SeedSAalign> > >* results;
    uint st;
    uint ed;
    MmpProperties mmpPara;
} MMP_ARG;

void* mmp_worker(void* mmp_args) {
	MMP_ARG* args = (MMP_ARG*) mmp_args;
  for (int i=args->st;i<args->ed;++i) {
    int read_len = args->read_lens[i];

    mmp<0>(args->queries+i * args->wpq, (*(args->results))[1][i], args->lbwt, args->lkt, 0, read_len, args->text_len, args->mmpPara);
    mmp<2>(args->queries+i * args->wpq, (*(args->results))[2][i], args->lbwt, args->lkt, 0, read_len, args->text_len, args->mmpPara);
  }
}

uint64_t PairEndSeedingEngine::PairEndSeedingBatch::mmpSeeding( int words_per_query, SeedPos** readPos, SeedPos** matePos, int n_threads, SeedPool *seedPool, MmpProperties mmpPara) {
  BWT *lbwt = bwt;
  BWT *rbwt = engine->index->sraIndex->rev_bwt;

  unsigned long long text_len = lbwt->textLength;
  // build queries
  uint nQ0 = numQueries[0];
  uint nQ1 = numQueries[1];
  uint nQ = nQ0 + nQ1;
  uint nQ_round = nQ;
  uint *packedQueries = (uint*) malloc( nQ_round * words_per_query * sizeof(uint) ); // is actually unpacked queries!
  uint *read_lens = (uint*) malloc( nQ_round * sizeof(uint) );
  // pack even-numbered reads
  omp_set_num_threads(n_threads);

#pragma omp parallel for
  for (uint i = 0; i < nQ0; ++i) {
    int id = readIDs[0][i];
    uint *src_p = queries + (id / QUERY_BLOCK * QUERY_BLOCK * words_per_query + id % QUERY_BLOCK);
    uint *dst_p = packedQueries + words_per_query * i;
    for (uint j = 0; j < words_per_query; ++j) {
        *(dst_p + j) = *(src_p + QUERY_BLOCK * j);
    }
    read_lens[i] = readLengths[id];
  }
  // pack odd-numbered reads. note that read IDs give are even, so need +1
#pragma omp parallel for
  for (uint i = 0; i < nQ1; ++i) {
    int id = readIDs[1][i] + 1;
    uint *src_p = queries + (id / QUERY_BLOCK * QUERY_BLOCK * words_per_query + id % QUERY_BLOCK);
    uint *dst_p = packedQueries + words_per_query * (i+nQ0);
    for (uint j = 0; j < words_per_query; ++j) {
      *(dst_p + j) = *(src_p + QUERY_BLOCK * j);
    }
    read_lens[i+nQ0] = readLengths[id];
  }

  vector<vector<vector<SeedSAalign> > > results(4, vector<vector<SeedSAalign> >(nQ));
  
  pthread_t* threads = new pthread_t[n_threads];

  MMP_ARG* args = new MMP_ARG[n_threads];
  int pp = nQ / n_threads;
  for (int i=0;i<n_threads;++i) {
	args[i].queries = packedQueries;
	args[i].results = &results;
	args[i].lbwt = lbwt;
	args[i].rbwt = rbwt;
    args[i].lkt = engine->index->sraIndex->lookupTable;
    args[i].rlkt = engine->index->sraIndex->rev_lookupTable;
	args[i].read_lens = read_lens;
	args[i].text_len = text_len;
	args[i].wpq = words_per_query;
    args[i].st = pp * i ;
	args[i].ed = pp * (i+1);
    args[i].mmpPara = mmpPara;
  }
  args[n_threads-1].ed = nQ;

  for (int i=0;i<n_threads;i++) 
	pthread_create(&threads[i], NULL, mmp_worker, (void*) &args[i]);	

  for (int i=0;i<n_threads;++i)
	pthread_join(threads[i],NULL);

  vector<vector<SeedPos> > vReadPos(n_threads);
  vector<vector<SeedPos> > vMatePos(n_threads);
  vector<vector<SeedPos> > vReadNeg(n_threads);
  vector<vector<SeedPos> > vMateNeg(n_threads);

#pragma omp parallel for schedule(static) // must be static to reserve the order of read id
  for (int i = 0 ;i < nQ; ++i) {
    vector<SeedAlign> seed_aligns[4];
    for (int a = 1; a <= 2; ++a) { // a=1 pos; a=2 neg
      for (size_t j=0;j<results[a][i].size();++j) {
        auto &res = results[a][i][j];

        unsigned long long l = res.sa_l;
        unsigned long long r = res.sa_l + res.sa_diff;
        if (r > l + mmpPara.seedSAsizeThreshold) { r = mmpPara.seedSAsizeThreshold + l - 1; }

        uint off = res.query_offset;
        uint seedlen = res.seed_len;

        for (unsigned long long k=l; k<=r; ++k) {
          unsigned long long target_ofs = a == 1 ? (BWTSaValue (lbwt, k) - off) : (BWTSaValue (lbwt, k) - (read_lens[i] - seedlen - off));
          if (seedlen >= mmpPara.goodSeedLen || seedlen >= read_lens[i] / 2) { // treat it as unique if seed is long enough
            seed_aligns[a].push_back(SeedAlign(target_ofs, 1, seedlen, off));
          } else {
            seed_aligns[a].push_back(SeedAlign(target_ofs, r-l+1, seedlen, off));
          }
        }
      }
    } // end seed collection

    sort(seed_aligns[1].begin(), seed_aligns[1].end());
    sort(seed_aligns[2].begin(), seed_aligns[2].end());

    vector<SeedPos> s_pos;
    uint max_seed_len = 0;
    
    // merge the seed alignments
    for (int a = 1; a <= 2; ++a) {
        for (size_t m = 0 ; m < seed_aligns[a].size(); ++m) {
          bool has_unique = seed_aligns[a][m].multiplicity <= mmpPara.uniqThreshold && seed_aligns[a][m].length >= mmpPara.seedMinLength;
          SeedPos sp;
          sp.pos = seed_aligns[a][m].offset;

          vector<pair<uint, uint> > itv(1, std::make_pair((uint)seed_aligns[a][m].query_offset, seed_aligns[a][m].length + seed_aligns[a][m].query_offset));

          while (m+1 < seed_aligns[a].size() && seed_aligns[a][m+1].offset <= sp.pos + mmpPara.indelFuzz) {
            ++m;
            has_unique |= seed_aligns[a][m].multiplicity <= mmpPara.uniqThreshold && seed_aligns[a][m].length >= mmpPara.seedMinLength;
            itv.push_back(std::make_pair((uint)seed_aligns[a][m].query_offset, seed_aligns[a][m].length + seed_aligns[a][m].query_offset));
          }

          sort(itv.begin(), itv.end());
          uint total_len = 0;
          uint cur_start = 0;
          uint cur_end = 0;
          for (uint i = 0; i < itv.size(); ++i) {
            if (itv[i].first >= cur_end) {
                total_len += cur_end - cur_start;
                cur_start = itv[i].first;
            }
            cur_end = max(cur_end, itv[i].second);
          }
          total_len += cur_end - cur_start;

          if (max_seed_len < total_len) { max_seed_len = total_len; }

          if (has_unique || total_len >= mmpPara.goodSeedLen) { // filtering criteria
            sp.paired_seedLength = total_len;
            sp.strand_readID = i < nQ0 ? readIDs[0][i] : readIDs[1][i-nQ0];
            sp.strand_readID |= (a != 1) * 1ULL << 31;
            s_pos.push_back(sp);
          }
        }
    }

    uint strand_readID = i < nQ0 ? readIDs[0][i] : readIDs[1][i-nQ0];
    int tid = omp_get_thread_num();
    for (size_t j = 0; j < s_pos.size(); ++j) {
        if (s_pos[j].paired_seedLength >= mmpPara.shortSeedRatio * max_seed_len) {
            bool isNeg = (s_pos[j].strand_readID >> 31);
            auto &v = i < nQ0 ? (isNeg ? vReadNeg[tid] : vReadPos[tid]) : (isNeg ? vMateNeg[tid] : vMatePos[tid]);
            v.push_back(s_pos[j]);
        }
    }
  }

  // sentinels
  SeedPos t;
  t.pos = kMaxULL;
  t.paired_seedLength = 0xffffffff;

  // merge all results
  uint n_read_pos = 0, n_mate_pos = 0;
  for (int i = 0; i < n_threads; ++i) {
    n_read_pos += vReadPos[i].size() + vReadNeg[i].size();
    n_mate_pos += vMatePos[i].size() + vMateNeg[i].size();
  }

  *readPos = (SeedPos*) malloc( (n_read_pos + 2) * sizeof(SeedPos) );
  *matePos = (SeedPos*) malloc( (n_mate_pos + 2) * sizeof(SeedPos) );

  uint nr_cur = 0, mr_cur = 0;

  for (int i = 0; i < n_threads; ++i) {
    memcpy(*readPos + nr_cur, &vReadPos[i][0], vReadPos[i].size() * sizeof(SeedPos));
    nr_cur += vReadPos[i].size();

    memcpy(*matePos + mr_cur, &vMatePos[i][0], vMatePos[i].size() * sizeof(SeedPos));
    mr_cur += vMatePos[i].size();
  }

  t.strand_readID = 0x7fffffff;
  (*readPos)[nr_cur++] = t;
  (*matePos)[mr_cur++] = t;

  for (int i = 0; i < n_threads; ++i) {
    memcpy(*readPos + nr_cur, &vReadNeg[i][0], vReadNeg[i].size() * sizeof(SeedPos));
    nr_cur += vReadNeg[i].size();

    memcpy(*matePos + mr_cur, &vMateNeg[i][0], vMateNeg[i].size() * sizeof(SeedPos));
    mr_cur += vMateNeg[i].size();
  }

  t.strand_readID = 0xffffffff;
  (*readPos)[nr_cur++] = t;
  (*matePos)[mr_cur++] = t;

  fprintf(stderr, "[%s] read: %u, mate %u\n", __func__, n_read_pos, n_mate_pos);
  free(read_lens);
  free(packedQueries);
  delete[] threads;
  delete[] args;

  double start_time = setStartTime();
  // add to seed pool
  for (int i = 0; i < n_threads; ++i) {
    addSeedPosToSeedPool(seedPool, &vReadPos[i][0], 0, 0, vReadPos[i].size());
    addSeedPosToSeedPool(seedPool, &vReadNeg[i][0], 0, 0, vReadNeg[i].size());

    addSeedPosToSeedPool(seedPool, &vMatePos[i][0], 1, 0, vMatePos[i].size());
    addSeedPosToSeedPool(seedPool, &vMateNeg[i][0], 1, 0, vMateNeg[i].size());
  }
  fprintf(stderr, "[%s] %u %u %u %u\n", __func__, vReadPos[0].size(), vReadNeg[0].size(), vMatePos[0].size(), vMateNeg[0].size());
  fprintf(stderr, "[%s] Add to seed pool: %f\n", __func__, getElapsedTime(start_time));

  return (uint64_t(n_mate_pos + 1) << 32 | (n_read_pos + 1));
}

void PairEndSeedingEngine::performMmpSeeding()
{
    double start_time = setStartTime();
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    int word_per_query = getWordPerQuery( inputMaxReadLength );
    seedingSwapBatch =
        new PairEndSeedingBatch ( queryIDStream->data->size(), dpPara,
                                  queries, queryLengths, inputMaxReadLength,
                                  insert_high, insert_low,
                                  peStrandLeftLeg, peStrandRightLeg, index->sraIndex->bwt, this );
    seedingSwapBatch->packForMmpSeeding( queryIDStream, 0 );
    seedingSwapBatch->packForMmpSeeding( queryIDStream, 1 );

    fprintf(stderr, "[%s] Pre-seeding time: %f\n", __func__, getElapsedTime(start_time));
    start_time = setStartTime();
    // Necessary Parameters && Input
    // ======================================================================================
    //   [Variable Name]: [Description]
    // ======================================================================================
    //   seedingSwapBatch: it is PairEndSeedingBatch which has set the readIDs and numQueries
    //   this->soap3Wrapper: contains _bwt, _revBwt, _occ, _revOcc which are already malloc in GPU
    //   word_per_query: as its name
    // ======================================================================================

    // Expected Output
    // ======================================================================================
    //  struct SeedPos
    //  {
    //      uint pos;                  // position
    //      uint strand_readID;        // ( (strand<<31) | evenReadID ), strand: 0 = '+' and 1 = '-'
    //      uint paired_seedLength;    // ( (0<<31) | seedLength), seedLength = seed alignment Length, paired part is determined in pairing stage 
    //  };
    // ======================================================================================
    // SeedPos * readPos;
    // SeedPos * matePos;
    // 
    // Note:
    //   Please push two SeedPos, one with strand_readID=0x7fffffff and one with strand_readID=0xffffffff into the end of readPos
    //   and also for matePos
    //   They will be then used to seperate SeedPos by strand and the length of readPos/matePos can be found
    //   
    
    // Chun will write the code for sorting them

    // Maybe we can call mmpSeeding as follows
    SeedPos * readPos;
    SeedPos * matePos;
    uint64_t arrayLens = seedingSwapBatch->mmpSeeding( word_per_query, &readPos, &matePos, dpPara->numOfCPUThreads, engine->seedPool, dpPara->mmpProperties);


    fprintf(stderr, "[%s] Seeding time: %f\n", __func__, getElapsedTime(start_time));
    start_time = setStartTime();

    vector<CandidateInfo> * candidates = seedingSwapBatch->mergeAndPairPairedEnd ( readPos, matePos, arrayLens & 0xffffffff, arrayLens >> 32);

    fprintf(stderr, "[%s] MergePE time: %f\n", __func__, getElapsedTime(start_time));
    start_time = setStartTime();

    engine->canStream->append ( candidates, engine->alignFlags );
    alignFlags->getXOR ( inputFlags, unseededIDStream->data );


    fprintf(stderr, "[%s] Post-seeding time: %f\n", __func__, getElapsedTime(start_time));
    start_time = setStartTime();

    delete candidates;
    delete inputFlags;
    delete alignFlags;
    delete seedingSwapBatch;
}

void PairEndSeedingEngine::performMmpSeeding (
    /* input */
    QueryIDStream    *    queryIDStream,
    DPParameters     *    dpPara,
    uint * queries, uint * queryLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    /* soap3 seeding related */
    SOAP3Wrapper<void>  * soap3Wrapper,
    Soap3Index      *     index,
    SeedPool * seedPool,
    /* output */
    CandidateStream   *   canStream,
    QueryIDStream    *    unseededIDStream
)
{
    engine = new PairEndSeedingEngine();
    MC_MemberCopy5 ( engine->, , queryIDStream, dpPara, queries, queryLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy4 ( engine->, , soap3Wrapper, index, canStream, unseededIDStream );
    engine->seedPool = seedPool;
    engine->performMmpSeeding();
    delete engine;
}

PairEndSeedingEngine * PairEndSeedingEngine::engine;

void DeepDP_Space::SeedingGPUThreadInit()
{
    //  showGPUMemInfo("seeding enter");
}
void DeepDP_Space::SeedingGPUThread ( int threadId, int *& pCallThreadId )
{
    assert(false);
}
void DeepDP_Space::SeedingGPUThreadFinalize()
{
    //  showGPUMemInfo("seeding exit");
}
void DeepDP_Space::SeedingCPUThread ( int threadId, void *& empty )
{
    assert(false);
}

void DeepDP_Space::DPExtensionPassUnalignedReads (
        QueryIDStream * unalignedIDStream,
        uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
        Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
        char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
        uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
        unsigned int maxHitNumForDP, unsigned int maxHitNumForDP2,
        ReadInputForDP * dpInput, ReadInputForDP * dpInputForNewDefault,
        ReadInputForDP * otherSoap3Result, BothUnalignedPairs * bothUnalignedPairs)
{
    assert(false);
}


void DeepDP_Space::DP2OutputUnalignedReads (
    QueryIDStream * unalignedIDStream, SeedPool * seedPool,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr
)
{
    // output unaligned result
#define MC_DP2OutputUnalgnRead() { \
        outputDeepDPResult2(buf, idx, \
                            queries, upkdReadLengths, \
                            origReadIDs, queryNames, queryComments, upkdQualities, \
                            inputMaxReadLength, accumReadNum, outputFormat, \
                            samOutputDPFilePtr, index, \
                            peStrandLeftLeg, peStrandRightLeg, NULL); }
    DeepDPAlignResult * buf;
    MC_CheckMalloc ( buf, DeepDPAlignResult, 1024 );
    int idx = 0;

    for ( uint i = 0; i < unalignedIDStream->data->size(); i++ )
    {
        buf[idx].readID = ( * ( unalignedIDStream->data ) ) [i];
        buf[idx].algnmt_1 = kMaxULL;
        buf[idx].algnmt_2 = kMaxULL;
        buf[idx].cigarString_1 = NULL;
        buf[idx].cigarString_2 = NULL;
        ++idx;

        if ( idx >= 1024 )
        {
            MC_DP2OutputUnalgnRead();
            idx = 0;
        }
    }

    if ( idx > 0 )
    { MC_DP2OutputUnalgnRead(); }

    free ( buf );
}

// ****
PairEndAlignmentEngine::PairEndAlgnBatch::PairEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths,
    DPInfoForReads * dpInfoForReads
)
{
    MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy4 ( this->, , peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low );
    MC_MemberCopy3 ( this->, , queries, inputMaxReadLength, upkdLengths );
    MC_MemberCopy ( this->, , dpInfoForReads );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->wordPerOldQuery   = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery      = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA        = MC_CeilDivide16 ( maxDNALength );
    this->packedDNA         = index->sraIndex->hsp->packedDNA;
    this->fullDNALength     = index->sraIndex->hsp->dnaLength;
    this->index             = index;
    MC_CheckMalloc ( packedDNASeq,        uint,           batchSize * wordPerDNA );
    MC_CheckMalloc ( packedReadSeq,       uint,           batchSize * wordPerQuery );
    MC_CheckMalloc ( canInfos,            CandidateInfo,  batchSize );
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        MC_CheckMalloc ( scores[lOr],     int,            batchSize );
        MC_CheckMalloc ( hitLocs[lOr],    uint,           batchSize );
        MC_CheckMalloc ( pattern[lOr],    uchar,          batchSize * patternLength );
        MC_CheckMalloc ( maxScoreCounts[lOr],     uint,           batchSize );
    }

    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( peLeftAnchorLocs,    uint,           batchSize );
    MC_CheckMalloc ( peRightAnchorLocs,   uint,           batchSize );
    clear();
}
PairEndAlignmentEngine::PairEndAlgnBatch::~PairEndAlgnBatch()
{
    free ( packedDNASeq );
    free ( packedReadSeq );
    free ( canInfos );
    free ( DNALengths );
    free ( lengths );
    free ( cutoffThresholds );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        free ( scores[lOr] );
        free ( hitLocs[lOr] );
        free ( pattern[lOr] );
        free ( maxScoreCounts[lOr] );
    }

    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( peLeftAnchorLocs );
    free ( peRightAnchorLocs );
}
void PairEndAlignmentEngine::PairEndAlgnBatch::clear()
{
    numOfThreads = 0;
}

int PairEndAlignmentEngine::PairEndAlgnBatch::packLeft (
    CandidateInfo & canInfo
)
{
    if ( numOfThreads >= batchSize )
    {
        return 0;
    }

    uint readIDLeft = canInfo.readIDLeft;
    uint readLength = upkdLengths[readIDLeft];
    int margin = DP2_MARGIN ( readLength );
    unsigned long long DNAStartLeft = ( engine->dpPara->isExtensionDP ? canInfo.pos[0] : ( canInfo.pos[0] - margin ) );

    if ( DNAStartLeft >= fullDNALength )
    {
        DNAStartLeft = 0;
    }

    uint DNALength = readLength + margin * 2;

    if ( DNAStartLeft + DNALength > fullDNALength )
    {
        DNALength = fullDNALength - DNAStartLeft;
    }

    // no anchor requirement
    peLeftAnchorLocs[numOfThreads] = maxDNALength;
    peRightAnchorLocs[numOfThreads] = 0;

    if ( engine->dpPara->isExtensionDP )
    {
        packRead ( packedReadSeq, numOfThreads,
                   readIDLeft, readLength, 1 );
        if ( peStrandLeftLeg == 1 )
        {
            repackDNA ( packedDNASeq, numOfThreads,
                        packedDNA, DNAStartLeft, DNALength );
        }
        else
        {
            repackDNA_Reverse_Complement ( packedDNASeq, numOfThreads, 
                                           packedDNA, DNAStartLeft + readLength, DNALength);
            DNAStartLeft += readLength;
        }
    }
    else
    {
        packRead ( packedReadSeq, numOfThreads,
                   readIDLeft, readLength, peStrandLeftLeg );
        repackDNA ( packedDNASeq, numOfThreads,
                    packedDNA, DNAStartLeft, DNALength );
    }
    softClipLtSizes[numOfThreads] = ( peStrandLeftLeg == 1 ) ?
                                    softClipLeft : softClipRight;
    softClipRtSizes[numOfThreads] = ( peStrandLeftLeg == 1 ) ?
                                    softClipRight : softClipLeft;
    DNALengths[numOfThreads] = DNALength;
    lengths[numOfThreads] = readLength;
    cutoffThresholds[numOfThreads] = std::max(DP_SCORE_THRESHOLD_RATIO * readLength, DP_SCORE_THRESHOLD_LOWER_BOUND);
    canInfo.pos[0] = DNAStartLeft;
    canInfo.dnaLength[0] = DNALength;
    canInfo.peLeftAnchor[0] = peLeftAnchorLocs[numOfThreads];
    canInfo.peRightAnchor[0] = peRightAnchorLocs[numOfThreads];
    canInfos[numOfThreads] = canInfo;
    ++numOfThreads;
    return 1;
}

void PairEndAlignmentEngine::PairEndAlgnBatch::packRight()
{
    for ( int i = 0; i < numOfThreads; i++ )
    {
        uint readIDLeft = canInfos[i].readIDLeft;
        uint leftIsOdd = readIDLeft & 1;

        if ( scores[0][i] >= cutoffThresholds[i] )
        {
            uint readIDLeft = canInfos[i].readIDLeft;
            uint readIDRight = ( leftIsOdd ) ? ( readIDLeft - 1 ) : ( readIDLeft + 1 );
            uint readLength = upkdLengths[readIDRight];
            // fprintf(stderr, "%d: %u\n", readIDRight, readLength);
            uint margin = DP2_MARGIN ( readLength );
            unsigned long long DNAStartRight = ( engine->dpPara->isExtensionDP ? canInfos[i].pos[1] : ( canInfos[i].pos[1] - margin ) );

            if ( DNAStartRight >= fullDNALength )
            {
                DNAStartRight = 0;
            }

            uint DNALength = readLength + margin * 2;

            if ( DNAStartRight + DNALength > fullDNALength )
            {
                DNALength = fullDNALength - DNAStartRight;
            }

            unsigned long long hitPosLeft = canInfos[i].pos[0] + hitLocs[0][i];
            // restrict maximum insert size
            unsigned long long boundedLength = hitPosLeft + insert_high - DNAStartRight;

            if ( boundedLength < DNALength ) \
                DNALength = boundedLength;

            if ( engine->dpPara->isExtensionDP )
            {
                // set pair-end anchor boundary, restrict minimum insert size
                
                packRead ( packedReadSeq, i,
                           readIDRight, readLength, 1 );
                if ( peStrandRightLeg == 1 )
                {
                    repackDNA ( packedDNASeq, i,
                                packedDNA, DNAStartRight, DNALength );
                }
                else
                {
                    repackDNA_Reverse_Complement ( packedDNASeq, i, 
                                                   packedDNA, DNAStartRight + readLength, DNALength);
                    DNAStartRight += readLength;
                }
                
                // assume no anchor for DP Extension
                peLeftAnchorLocs[i] = maxDNALength;
                peRightAnchorLocs[i] = 0;
            }
            else
            {
                // set pair-end anchor boundary, restrict minimum insert size
                peLeftAnchorLocs[i] = maxDNALength;
                long long rightAnchor = hitPosLeft + insert_low - DNAStartRight;
                peRightAnchorLocs[i] = rightAnchor > 0 ? rightAnchor : 0;
                packRead ( packedReadSeq, i,
                           readIDRight, readLength, peStrandRightLeg );
                repackDNA ( packedDNASeq, i,
                            packedDNA, DNAStartRight, DNALength );
            }
            softClipLtSizes[i] = ( peStrandRightLeg == 1 ) ?
                                     softClipLeft : softClipRight;
            softClipRtSizes[i] = ( peStrandRightLeg == 1 ) ?
                                     softClipRight : softClipLeft;
            DNALengths[i] = DNALength;
            lengths[i] = readLength;
            cutoffThresholds[i] = std::max(DP_SCORE_THRESHOLD_RATIO * readLength, DP_SCORE_THRESHOLD_LOWER_BOUND);
            canInfos[i].pos[1] = DNAStartRight;
            canInfos[i].dnaLength[1] = DNALength;
            canInfos[i].peLeftAnchor[1] = peLeftAnchorLocs[i];
            canInfos[i].peRightAnchor[1] = peRightAnchorLocs[i]; 
        }
    }
}

inline void PairEndAlignmentEngine::PairEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}

inline void PairEndAlignmentEngine::PairEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, unsigned long long start, uint length
)
{
    size_t dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    uint *src = seq + start / 16;
    uint src_ofs = start % 16;
    uint *dst = packedSeq + dnaTPARA;
    uint dst_ofs = 1;
    *dst = 0;

    // we need a 0 at the beginning to reserve the first column of DP table
    while (length > 0) {
        uint len = min(min(16 - dst_ofs, 16 - src_ofs), length);
        *dst |= *src << (src_ofs * 2) >> (32 - len*2) << (32-(dst_ofs+len)*2);
        length -= len;
        dst_ofs += len;
        src_ofs += len;

        if (src_ofs == 16) { ++src; src_ofs = 0; }
        if (dst_ofs == 16) { dst += 32; dst_ofs = 0; *dst = 0; }
    }
}

// Variable:
//      start : is the last position of positive strand DNA
inline void PairEndAlignmentEngine::PairEndAlgnBatch::repackDNA_Reverse_Complement (
    uint * packedSeq, uint threadId,
    uint * seq, unsigned long long start, uint length
)
{
#define MC_OldDnaUnpack(X,i) ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
    size_t dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[dnaTPARA + ( i << 5 )] = 0;
    }
    const char dnaMap[] = {'A','C','G','T'};
    for ( int i = 1; i <= length; i++ )
    { packedSeq[dnaTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) ( ( 3 - MC_OldDnaUnpack ( seq, start - i ) ) ) << ( ( 15 - ( i & 0xF ) ) << 1 ); } 
}

// ****
void PairEndAlignmentEngine::PairEndAlgnThreadContext::init (
    PairEndAlgnBatch * batch
)
{
    int batchSize = engine->DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                             engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
    batchID = -1;
    this->dp2AlignedRead = 0;
    this->dp2Alignment = 0;
    this->lastReadID = -1;
    this->alignFlag = new AlgnmtFlags;
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}
void PairEndAlignmentEngine::PairEndAlgnThreadContext::freeMemory()
{
    semiGlobalAligner.freeMemory();
    delete this->alignFlag;
    delete batch;
}

// ****
PairEndAlignmentEngine::AlgnmtResultStream::AlgnmtResultStream()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}

PairEndAlignmentEngine::AlgnmtResultStream::~AlgnmtResultStream()
{
    for ( int i = 0; i < dp2Result.size(); i++ )
    {
        DP2ResultBatch & resultBatch = * ( dp2Result[i] );

        for ( int j = 0; j < resultBatch.size(); j++ )
        {
            free ( resultBatch[j].cigarString_1 );
            free ( resultBatch[j].cigarString_2 );
        }

        delete dp2Result[i];
    }

    dp2Result.clear();
}

// ****
PairEndAlignmentEngine::PairEndAlignmentEngine() {}

void PairEndAlignmentEngine::performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment )
{
    /* initialize */
    algnBatchCount = 0;
    dp2AlignedRead = 0;
    dp2Alignment = 0;
    lastReadID = -1;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    outputBuf = new OutputBuffer<DeepDPAlignResult>();
    /**/
    resultStreams = ( AlgnmtResultStream ** ) malloc ( dpPara->numOfCPUThreads * sizeof ( AlgnmtResultStream * ) );
    outputBufs = ( OutputBuffer<DeepDPAlignResult> ** ) malloc ( dpPara->numOfCPUThreads * sizeof ( OutputBuffer<DeepDPAlignResult> * ) );
    resultOutputStrings = ( char *** ) malloc ( dpPara->numOfCPUThreads * sizeof ( char ** ) );
    numberOfResultOutputString = ( int * ) malloc ( dpPara->numOfCPUThreads * sizeof ( int ) );
    /**/
    outputBuf->setAlignmentType ( alignmentType );
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = maxReadLength + 2 * DP2_MARGIN ( inputMaxReadLength ) + 8;
    
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DP2_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );

    resultOutputStringBufferIter = 0;
    resultOutputStringBufferSize = DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK * 2 * 7;
    resultOutputStringBuffer = ( char ** ) malloc ( resultOutputStringBufferSize * sizeof ( char * ) );
    for ( int i=0; i<dpPara->numOfCPUThreads; ++i )
    { 
        resultStreams[i] = new AlgnmtResultStream;
        outputBufs[i] = new OutputBuffer<DeepDPAlignResult>();
        outputBufs[i]->setAlignmentType( alignmentType );
        resultOutputStrings[i] = ( char **  ) malloc ( DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK * 2 * sizeof ( char * ) );
        numberOfResultOutputString[i] = 0;
    }

    algnSwapBatch =
        new PairEndAlgnBatch ( DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                               peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                               maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                               index, queries, inputMaxReadLength, upkdReadLengths,
                               dpInfoForReads );
    algnThreadContext = new PairEndAlgnThreadContext[dpPara->numOfCPUThreads];
    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        PairEndAlgnBatch * batch =
            new PairEndAlgnBatch ( DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                   peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                                   maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                   index, queries, inputMaxReadLength, upkdReadLengths,
                                   dpInfoForReads
                                  );
        algnThreadContext[i].init ( batch );
    }
    
    if ( !( samOutputDPFilePtr != NULL && samOutputDPFilePtr->type & TYPE_BAM ) )
    {
        outputThreadDelegator.init ( 1, DP2OutputThread2,
                                     NULL, DP2OutputThread2Finalize );
    }
    else 
    {
        outputThreadDelegator.init ( 1, DP2OutputThread,
                                     NULL, DP2OutputThreadFinalize );
    }
    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, DP2CPUAlgnThread );
    /* perform alignment */
    int threadId;
    void * empty;
    int cntRound=0,cntPairs=0,cntP=0;
    double startTime=setStartTime(), lastEventTime=0,elapsedTime, maxWait=0, minWait=100, tmp, totalWait=0;
    int numberOfDPInstancesRemainedOfCurrentReadPair = 0;
    int lastReadId=-1;
    fprintf(stderr, "[PairEndAlignmentEngine::performAlignment] Candidate: %d\n", canStream->data.size());
    for ( uint i = 0; i < canStream->data.size(); i++ )
    {
        if ( lastReadId != canStream->data[i].readIDLeft )
            { lastReadId = canStream->data[i].readIDLeft; ++cntP; }
        if ( numberOfDPInstancesRemainedOfCurrentReadPair == 0 )
        { 
            uint j=i;
            for ( j=i; j < canStream->data.size() && ( canStream->data[j].readIDLeft & 0xfffffffe ) == (canStream->data[i].readIDLeft & 0xfffffffe); ++j );
            numberOfDPInstancesRemainedOfCurrentReadPair = j-i;
            ++cntPairs;
        }
        CandidateInfo & info = canStream->data[i];
        inputFlags->set ( ( info.readIDLeft >> 1 ) << 1 );
#ifdef ALIGNMENT_DEBUG
        assert ( numberOfDPInstancesRemainedOfCurrentReadPair <= algnSwapBatch->batchSize );
#else
        if ( numberOfDPInstancesRemainedOfCurrentReadPair > algnSwapBatch->batchSize )
        {
            fprintf(stderr, "[PairEndAlignmentEngine::performAlignment] Wrong configuration between GPU batch size and Max. # of seed alignment pairs for each read pair\n");
        }
#endif
        uint availableNumberOfCaninfoForThisBatch = algnSwapBatch->batchSize - algnSwapBatch->numOfThreads;
        if ( numberOfDPInstancesRemainedOfCurrentReadPair > availableNumberOfCaninfoForThisBatch )
        {
            lastEventTime = getElapsedTime(startTime) ;
            // launch one batch
            threadId = algnmtCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
            tmp = getElapsedTime (startTime);
            if ( maxWait < tmp - lastEventTime )
                maxWait = tmp - lastEventTime;
            if ( minWait > tmp - lastEventTime )
                minWait = tmp -lastEventTime;
            totalWait += tmp - lastEventTime;
            lastEventTime = tmp;
            ++cntRound;
            algnSwapBatch->clear();
        }
        
        algnSwapBatch->packLeft ( info ) ;
        --numberOfDPInstancesRemainedOfCurrentReadPair;  

    }
    /*
    fprintf(stderr, "Min wait: %f seconds\n", minWait);
    fprintf(stderr, "Max wait: %f seconds\n", maxWait);
    fprintf(stderr, "Total wait: %f seconds\n", totalWait);
    fprintf(stderr, "Avg wait: %f seconds\n", totalWait*1.0/cntRound);
    fprintf(stderr, "%d rounds for %d pairs among %d pairs\n", cntRound, cntPairs, cntP);
    */
    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
    }
    /* finalize */
    algnmtCPUThreadDelegator.finalize();
    outputThreadDelegator.finalize();
    /**/
    for ( int i=0;i<dpPara->numOfCPUThreads;++i )
        { alignFlags->XOR ( algnThreadContext[i].alignFlag ); }
    /**/
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );

    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;

    /**/
    for ( int i=0;i<dpPara->numOfCPUThreads;++i )
    { 
        dp2AlignedRead += algnThreadContext[i].dp2AlignedRead;
        dp2Alignment += algnThreadContext[i].dp2Alignment;
        delete outputBufs[i]; 
        delete resultStreams[i];
        free ( resultOutputStrings[i] );
    }
    free ( outputBufs );
    free ( resultStreams );
    free ( resultOutputStrings );
    free ( resultOutputStringBuffer );
    free ( numberOfResultOutputString );
    /**/
    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory();
    }
    delete[] algnThreadContext;
    delete outputBuf;
    delete resultStream;
    // fprintf(stderr, "GPU %f seconds\n",gpuAlignmentTime);
    numDPAlignedRead = this->dp2AlignedRead;
    numDPAlignment = this->dp2Alignment;
}

void PairEndAlignmentEngine::performAlignment (
    /* input */
    CandidateStream   *   canStream,
    DPParameters     *    dpPara,
    bool shortDnaLength,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    DPInfoForReads * dpInfoForReads, 
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, SeedPool * seedPool, samfile_t ** samOutputDPFilePtrs,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment
)
{
    engine = new PairEndAlignmentEngine();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, queryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy ( engine->, , queryComments);
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , dpInfoForReads );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy4 ( engine->, , accumReadNum, outputFormat, samOutputDPFilePtrs, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    MC_MemberCopy2 ( engine->, , shortDnaLength, seedPool);
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    delete engine;
}
PairEndAlignmentEngine * PairEndAlignmentEngine::engine;

// ****

inline void DeepDP_Space::reverseCigarString ( char * cigar )
{
#define SWAP_TWO_CHAR(a,b) (a^=b^=a^=b)
    // reverse the cigar String
    int length = strlen (cigar);
    int i, j;
    for ( i=0;i<(length>>1);++i ) 
        { SWAP_TWO_CHAR ( cigar[i], cigar[length-i-1] ); }

    for ( i=length-1, j = length; i >= 0; --i )
    { 
        if ( !( '0' <= cigar[i] && cigar[i] <= '9' ) )
        {
            for ( int ll = 0; ll < ( j-i ) >> 1 ; ++ll )
                { SWAP_TWO_CHAR ( cigar[i+ll], cigar[j-ll-1] ); }
            j = i;
        }
    }
}

inline int DeepDP_Space::cigarAlignmentLength_Extension ( char * cigar )
{
    int l=0, range=0;
    for ( int i=0;cigar[i];++i )
    {
        if ( '0' <= cigar[i] && cigar[i] <= '9' )
            { range = range * 10 + cigar[i] - '0'; }
        else if ( cigar[i] == 'M' || cigar[i] == 'D' || ( cigar[i+1] == 0 && cigar[i] == 'S' ) )
            { l += range; range = 0; }
        else 
            { range = 0; }
    }
    return l;
}


void DeepDP_Space::DP2CPUAlgnThread ( int threadId, void *& empty )
{
//fprintf(stderr,"deepdp cpu thread %d start\n", threadId);
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    PairEndAlignmentEngine::PairEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    if ( !( engine->samOutputDPFilePtr != NULL && engine->samOutputDPFilePtr->type & TYPE_BAM ) ) 
        { engine->algnThreadContext[threadId].batchID++; }
    else
        { engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++; }
    sem_post ( & ( engine->algnThreadContext[threadId].ACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    // align left side
    batch->leftOrRight = 0;
        engine->algnThreadContext[threadId].semiGlobalAligner.performAlignment (
            batch->packedDNASeq, batch->DNALengths,
            batch->packedReadSeq, batch->lengths,
            batch->cutoffThresholds, batch->scores[batch->leftOrRight], batch->hitLocs[batch->leftOrRight],
            batch->maxScoreCounts[batch->leftOrRight],
            batch->pattern[batch->leftOrRight], batch->numOfThreads, 
            batch->softClipLtSizes, batch->softClipRtSizes,
            batch->peLeftAnchorLocs, batch->peRightAnchorLocs
        );    

    // align right side
    batch->packRight();
    batch->leftOrRight = 1;
    engine->algnThreadContext[threadId].semiGlobalAligner.performAlignment (
            batch->packedDNASeq, batch->DNALengths,
            batch->packedReadSeq, batch->lengths,
            batch->cutoffThresholds, batch->scores[batch->leftOrRight], batch->hitLocs[batch->leftOrRight],
            batch->maxScoreCounts[batch->leftOrRight],
            batch->pattern[batch->leftOrRight], batch->numOfThreads, 
            batch->softClipLtSizes, batch->softClipRtSizes,
            batch->peLeftAnchorLocs, batch->peRightAnchorLocs
        );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    // rearrange result and Output
    vector<DeepDPAlignResult> * resultBatch = new vector<DeepDPAlignResult>();
    for ( int i = 0; i < batch->numOfThreads; i++ )
    {
        int readSide = batch->canInfos[i].readIDLeft & 1;
        int mateSide = 1 - readSide;

        // DO NOT USE values in batch->packedDNASeq, packedReadSeq, DNALengths, length, cutoffThresholds
        // they are overwritten in packRight()
        if ( batch->scores[0][i] >= (int)std::max(DP_SCORE_THRESHOLD_RATIO * engine->upkdReadLengths[batch->canInfos[i].readIDLeft], DP_SCORE_THRESHOLD_LOWER_BOUND) &&
                batch->scores[1][i] >= (int)std::max(DP_SCORE_THRESHOLD_RATIO * engine->upkdReadLengths[batch->canInfos[i].readIDLeft ^1], DP_SCORE_THRESHOLD_LOWER_BOUND))
        {
            char * cigarString[2];
            int editdist[2], DIS[2];
            
            for ( int lOr = 0; lOr < 2; lOr++ )
            {
                CigarStringEncoder<void> encoder;
                uchar lastType = 'N';

                for ( uchar * p = batch->pattern[lOr] + i * engine->patternLength; *p != 0; p++ )
                {
                    if ( *p == 'V' )
                    {
                        encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                    }
                    else
                    {
                        encoder.append ( *p, 1 );
                        lastType = *p;
                    }
                }

                encoder.encodeCigarString ( openGapScore, extendGapScore );
                cigarString[lOr] = encoder.cigarString;
                // To get edit distance
                int L = batch->lengths[i] - encoder.charCount['I'] - encoder.charCount['S'];
                int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[lOr][i] ) /
                                    ( matchScore - mismatchScore );
                editdist[lOr] = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
                DIS[lOr] = encoder.charCount['D'] - encoder.charCount['I'] - encoder.charCount['S'];
            }
            
            //#define MC_GetMateID(x) (((x)&1)?((x)-1):((x)+1))
            DeepDPAlignResult result;
            result.readID = batch->canInfos[i].readIDLeft - readSide;
            result.strand_1 = ( ( readSide == 0 ) ? engine->peStrandLeftLeg : engine->peStrandRightLeg );
            result.strand_2 = ( ( mateSide == 0 ) ? engine->peStrandLeftLeg : engine->peStrandRightLeg );

            if ( engine->dpPara->isExtensionDP )
            {
                if ( result.strand_1 != 1 )
                {
                    reverseCigarString ( cigarString[readSide] );
                    result.algnmt_1 = batch->canInfos[i].pos[readSide] - cigarAlignmentLength_Extension ( cigarString[readSide] );
                }
                else
                {
                    result.algnmt_1 = batch->canInfos[i].pos[readSide] + batch->hitLocs[readSide][i];
                }
                
                if ( result.strand_2 != 1 )
                {
                    reverseCigarString ( cigarString[mateSide] );        
                    result.algnmt_2 = batch->canInfos[i].pos[mateSide] - cigarAlignmentLength_Extension ( cigarString[mateSide] );
                }
                else 
                {
                    result.algnmt_2 = batch->canInfos[i].pos[mateSide] + batch->hitLocs[mateSide][i];
                }
            }
            else
            {
                result.algnmt_1 = batch->canInfos[i].pos[readSide] + batch->hitLocs[readSide][i];
                result.algnmt_2 = batch->canInfos[i].pos[mateSide] + batch->hitLocs[mateSide][i];
            }           
            result.cigarString_1 = cigarString[readSide];
            result.cigarString_2 = cigarString[mateSide];
            result.score_1 = batch->scores[readSide][i];
            result.score_2 = batch->scores[mateSide][i];
            result.editdist_1 = editdist[readSide];
            result.editdist_2 = editdist[mateSide];
            result.startPos_1 = batch->canInfos[i].pos[readSide];
            result.startPos_2 = batch->canInfos[i].pos[mateSide];
            result.refDpLength_1 = batch->canInfos[i].dnaLength[readSide];
            result.refDpLength_2 = batch->canInfos[i].dnaLength[mateSide];
            result.peLeftAnchor_1 = batch->canInfos[i].peLeftAnchor[readSide];
            result.peLeftAnchor_2 = batch->canInfos[i].peLeftAnchor[mateSide];
            result.peRightAnchor_1 = batch->canInfos[i].peRightAnchor[readSide];
            result.peRightAnchor_2 = batch->canInfos[i].peRightAnchor[mateSide];
            if ( result.algnmt_1 < result.algnmt_2 )
                result.insertSize = result.algnmt_2 - result.algnmt_1 +
                                    batch->lengths[i] + DIS[1];
            else
                result.insertSize = result.algnmt_1 - result.algnmt_2 +
                                    batch->lengths[i] + DIS[1];
            result.num_sameScore_1 = batch->maxScoreCounts[readSide][i]; //TODO
            result.num_sameScore_2 = batch->maxScoreCounts[mateSide][i];
            resultBatch->push_back ( result );
        }
    }
    
    // Output
    engine->algnThreadContext[threadId].resultBatch = resultBatch;
    int * pid = &threadId;
    if ( !( engine->samOutputDPFilePtr != NULL && engine->samOutputDPFilePtr->type & TYPE_BAM ) )
        { engine->DP2Output( pid ); }
    engine->outputThreadDelegator.schedule ( pid );
    sem_wait ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
//fprintf(stderr, "deepdp cpu thread %d end\n", threadId);
}

void DeepDP_Space::DP2OutputThread ( int threadId, int *& pCallThreadId )
{
    double startTime = setStartTime(), elapsedTime;
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    int callThreadId = *pCallThreadId;
    int batchID = engine->algnThreadContext[callThreadId].batchID;
    DP2ResultBatch * resultBatch = engine->algnThreadContext[callThreadId].resultBatch;
    sem_post ( & ( engine->algnThreadContext[callThreadId].outputACKSem ) );
    vector<DP2ResultBatch *> & dpResult = engine->resultStream->dp2Result;
    OCC * occ = OCCConstruct();

    while ( dpResult.size() <= batchID )
    {
        dpResult.push_back ( NULL );
    }
    DeepDPAlignResult * tmpBestAlgn = NULL;
    dpResult[batchID] = resultBatch;
#define MC_DP2OutputRead() { \
        engine->outputBuf->ready(); \
        if (engine->outputBuf->size > 0) { \
            occ->occPositionCacheCount = 0; \
            tmpBestAlgn = outputDeepDPResult2(engine->outputBuf->elements, engine->outputBuf->size, \
                                engine->queries, engine->upkdReadLengths, \
                                engine->origReadIDs, engine->queryNames, engine->queryComments, engine->upkdQualities, \
                                engine->inputMaxReadLength, engine->accumReadNum, engine->outputFormat, \
                                engine->samOutputDPFilePtr, engine->index, \
                                engine->peStrandLeftLeg, engine->peStrandRightLeg, engine->seedPool, occ ); \
            engine->dp2AlignedRead += 1; \
            engine->dp2Alignment += engine->outputBuf->size; \
            engine->alignFlags->set(engine->lastReadID << 1); \
            engine->dpInfoForReads->startPositions[tmpBestAlgn->readID] = tmpBestAlgn->startPos_1; \
            engine->dpInfoForReads->strand_dpLengths[tmpBestAlgn->readID] =  ( ( ( tmpBestAlgn->strand_1 - 1 ) << 15 ) | ( tmpBestAlgn->refDpLength_1 ) ); \
            engine->dpInfoForReads->peLeftAnchors[tmpBestAlgn->readID] = tmpBestAlgn->peLeftAnchor_1; \
            engine->dpInfoForReads->peRightAnchors[tmpBestAlgn->readID] = tmpBestAlgn->peRightAnchor_1; \
            engine->dpInfoForReads->startPositions[tmpBestAlgn->readID+1] = tmpBestAlgn->startPos_2; \
            engine->dpInfoForReads->strand_dpLengths[tmpBestAlgn->readID+1] =  ( ( ( tmpBestAlgn->strand_2 - 1 ) << 15 ) | ( tmpBestAlgn->refDpLength_2 ) ); \
            engine->dpInfoForReads->peLeftAnchors[tmpBestAlgn->readID+1] = tmpBestAlgn->peLeftAnchor_2; \
            engine->dpInfoForReads->peRightAnchors[tmpBestAlgn->readID+1] = tmpBestAlgn->peRightAnchor_2; \
        } \
    }
    uint numOut = engine->resultStream->numOut;
    while ( numOut < dpResult.size() && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        DP2ResultBatch & batch = *dpResult[numOut];
        for ( int i = 0; i < batch.size(); i++ )
        {
            DeepDPAlignResult & result = batch[i];
            int pairID = result.readID >> 1;
            if ( pairID != engine->lastReadID )
            {
                MC_DP2OutputRead();
                engine->outputBuf->clear();
                engine->lastReadID = pairID;
            }

            engine->outputBuf->add ( result );
        }

        ++numOut;
    }
    OCCFree( occ );
    engine->resultStream->numOut = numOut;
    elapsedTime = getElapsedTime ( startTime );
    // fprintf(stderr, "Output thread: %f %d\n",elapsedTime, callThreadId);
}

void DeepDP_Space::DP2OutputThreadFinalize()
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    // last read
    OCC * occ = OCCConstruct();
    DeepDPAlignResult * tmpBestAlgn = NULL;
    MC_DP2OutputRead();
    OCCFree( occ );
    engine->outputBuf->clear();
}

void PairEndAlignmentEngine::DP2Output ( int *& threadId )
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    int callThreadId = *threadId;  
    int batchID = engine->algnThreadContext[callThreadId].batchID;
    DP2ResultBatch * resultBatch = engine->algnThreadContext[callThreadId].resultBatch;
    vector<DP2ResultBatch *> & dpResult = engine->resultStreams[callThreadId]->dp2Result;
    engine->outputBufs[callThreadId]->clear();
    while ( dpResult.size() <= batchID )
    {
        dpResult.push_back ( NULL );
    }
    OCC * occ = OCCConstruct();
    DynamicUint8Array * charArray = DynamicUint8ArrayConstruct();
    DeepDPAlignResult * tmpBestAlgn = NULL;
    engine->numberOfResultOutputString[callThreadId] = 0;   
    dpResult[batchID] = resultBatch;
#define MC_DP2OutputReads(tid) { \
        engine->outputBufs[tid]->ready(); \
        if (engine->outputBufs[tid]->size > 0) { \
            OCCReset( occ );\
            DynamicUint8ArrayReset( charArray ); \
            tmpBestAlgn = outputDeepDPResult2(engine->outputBufs[tid]->elements, engine->outputBufs[tid]->size, \
                                engine->queries, engine->upkdReadLengths, \
                                engine->origReadIDs, engine->queryNames, engine->queryComments, engine->upkdQualities, \
                                engine->inputMaxReadLength, engine->accumReadNum, engine->outputFormat, \
                                engine->samOutputDPFilePtrs[tid], engine->index, \
                                engine->peStrandLeftLeg, engine->peStrandRightLeg, engine->seedPool, \
                                occ, charArray, engine->resultOutputStrings[tid] + engine->numberOfResultOutputString[tid], tid ); \
            engine->algnThreadContext[tid].dp2AlignedRead += 1; \
            engine->algnThreadContext[tid].dp2Alignment += engine->outputBufs[tid]->size; \
            engine->algnThreadContext[tid].alignFlag->set(engine->algnThreadContext[tid].lastReadID << 1); \
            engine->dpInfoForReads->startPositions[tmpBestAlgn->readID] = tmpBestAlgn->startPos_1; \
            engine->dpInfoForReads->strand_dpLengths[tmpBestAlgn->readID] =  ( ( ( tmpBestAlgn->strand_1 - 1 ) << 15 ) | ( tmpBestAlgn->refDpLength_1 ) ); \
            engine->dpInfoForReads->peLeftAnchors[tmpBestAlgn->readID] = tmpBestAlgn->peLeftAnchor_1; \
            engine->dpInfoForReads->peRightAnchors[tmpBestAlgn->readID] = tmpBestAlgn->peRightAnchor_1; \
            engine->dpInfoForReads->startPositions[tmpBestAlgn->readID+1] = tmpBestAlgn->startPos_2; \
            engine->dpInfoForReads->strand_dpLengths[tmpBestAlgn->readID+1] =  ( ( ( tmpBestAlgn->strand_2 - 1 ) << 15 ) | ( tmpBestAlgn->refDpLength_2 ) ); \
            engine->dpInfoForReads->peLeftAnchors[tmpBestAlgn->readID+1] = tmpBestAlgn->peLeftAnchor_2; \
            engine->dpInfoForReads->peRightAnchors[tmpBestAlgn->readID+1] = tmpBestAlgn->peRightAnchor_2; \
            engine->numberOfResultOutputString[tid] += 2; \
        } \
    }
    uint numOut = engine->resultStreams[callThreadId]->numOut;
    while ( numOut < dpResult.size() && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        DP2ResultBatch & batch = *dpResult[numOut];
        for ( int i = 0; i < batch.size(); i++ )
        {
            DeepDPAlignResult & result = batch[i];
            int pairID = result.readID >> 1;
            if ( pairID != engine->algnThreadContext[callThreadId].lastReadID )
            {
                MC_DP2OutputReads( callThreadId );
                engine->outputBufs[callThreadId]->clear();
                engine->algnThreadContext[callThreadId].lastReadID = pairID;
            }
            engine->outputBufs[callThreadId]->add ( result );
        }
        MC_DP2OutputReads( callThreadId );
        engine->outputBufs[callThreadId]->clear();
        ++numOut;
    }
    OCCFree( occ );
    DynamicUint8ArrayFree( charArray );
    // fprintf(stderr, "callThreadId: %d, NumOut: %d, batchID: %d #: %d\n",callThreadId,numOut,batchID, engine->numberOfResultOutputString[callThreadId]);
    engine->resultStreams[callThreadId]->numOut = numOut;
}

void DeepDP_Space::DP2OutputThread2 ( int threadId, int *& pCallThreadId )
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    int callThreadId = *pCallThreadId;

#define MC_DP2OutputAllSamStrings() do { \
        if (!engine->index->sraIndex->hspaux->megapathMode) { \
            for ( int i=0;i<engine->resultOutputStringBufferIter; ++i) \
            { \
                fputs ( engine->resultOutputStringBuffer[i], engine->samOutputDPFilePtr->x.tamw); fputc('\n', engine->samOutputDPFilePtr->x.tamw); \
                free ( engine->resultOutputStringBuffer[i] ); \
            } \
        } \
        engine->resultOutputStringBufferIter = 0; \
    } while (0);

    for ( int i=0; i < engine->numberOfResultOutputString[callThreadId]; ++i )
    { 
        if ( engine->resultOutputStringBufferIter == engine->resultOutputStringBufferSize ) 
            { MC_DP2OutputAllSamStrings(); }
        engine->resultOutputStringBuffer[engine->resultOutputStringBufferIter++] = engine->resultOutputStrings[callThreadId][i];
    }
   
    sem_post ( & ( engine->algnThreadContext[callThreadId].outputACKSem ) );
}

void DeepDP_Space::DP2OutputThread2Finalize()
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    MC_DP2OutputAllSamStrings ();
}


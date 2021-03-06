/*
 *
 *    DV-DPfunctions.h
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

#ifndef _DP_FUNCTIONS_H
#define _DP_FUNCTIONS_H

#include "AlgnResult.h"
#include "alignment.h"
#include "CPUfunctions.h"
#include "CPU_DPfunctions.h"

#include <sys/time.h>
#include <pthread.h>
#include <semaphore.h>

#include <string>
#include <stack>
#include <map>
#include <vector>

#define TYPE_BAM 1 // should be same as TYPE_BAM in "samtools-0.1.18/sam.c"
using namespace std;
  
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// For seeding ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
class SOAP3Wrapper
{
        Soap3Index * index;
    public:
        SOAP3Wrapper ( Soap3Index * index,
                       int indexInside = 1 )
        {
	           this->index = index;
        }
};



//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// For output ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

template <class ResultType>
int ScoreCompare ( ResultType & a, ResultType & b );
template <class ResultType>
bool ResultCompare ( const ResultType & a, const ResultType & b );

template <class ResultType>
int isValid ( ResultType & a )
{
    return 1;
}

template <class ResultType>
struct OutputBuffer
{
    vector<ResultType> rawBuffer;
    ResultType * elements;
    int capacity;
    int size;
    int alignmentType;
    public:
    OutputBuffer()
    {
        // srand(time(NULL));
        capacity = 128;
        elements = ( ResultType * ) malloc ( capacity * sizeof ( ResultType ) );
        clear();
    }
    ~OutputBuffer()
    {
        free ( elements );
    }
    void clear()
    {
        rawBuffer.clear();
        size = 0;
    }
    void setAlignmentType ( int alignmentType )
    {
        this->alignmentType = alignmentType;
    }
    void add ( ResultType & result )
    {
        rawBuffer.push_back ( result );
    }
    inline void filterBest()
    {
        if ( size == 0 ) { return; }

        ResultType maxElement = elements[0];
        int arrSize = 1;

        for ( int i = 1; i < size; i++ )
        {
            int retval = ScoreCompare ( elements[i], maxElement );

            if ( retval > 0 )
            {
                maxElement = elements[i];
                arrSize = 0;
            }
            else if ( retval < 0 )
            {
                continue;
            }

            elements[arrSize++] = elements[i];
        }

        size = arrSize;
    }
    inline void filterInvalid()
    {
        int arrSize = 0;

        for ( int i = 0; i < size; i++ )
        {
            if ( isValid ( elements[i] ) )
            {
                elements[arrSize++] = elements[i];
            }
        }

        size = arrSize;
    }
    inline void arrayCopy()
    {
        size = rawBuffer.size();

        if ( size > 0 )
        {
            if ( size >= capacity )
            {
                free ( elements );
                capacity = size * 2;
                elements = ( ResultType * ) malloc ( capacity * sizeof ( ResultType ) );
            }

            copy ( rawBuffer.begin(), rawBuffer.end(), elements );
        }
    }
    inline void arrayCopyNRemoveDuplicate()
    {
        size = rawBuffer.size();

        if ( size > 0 )
        {
            if ( size >= capacity )
            {
                free ( elements );
                capacity = size * 2;
                elements = ( ResultType * ) malloc ( capacity * sizeof ( ResultType ) );
            }

            int arrPointer = 0;
            sort ( rawBuffer.begin(), rawBuffer.end(), ResultCompare<ResultType> );
            elements[arrPointer] = rawBuffer[0];

            for ( uint i = 1; i < rawBuffer.size(); i++ )
            {
                if ( ResultCompare ( elements[arrPointer], rawBuffer[i] ) )
                {
                    elements[++arrPointer] = rawBuffer[i];
                }
            }

            size = arrPointer + 1;
        }
    }

    void ready ( int optionFlag = 0 )
    {
        // alignmentType :  1 -- all valid
        //                  2 -- all best
        //                  3 -- unique best
        //                  4 -- random best
        //-----------------------------------
        // alignmentType == 1 --> do nothing
#define DO_NOT_OUTPUT_HALF_ALIGNED  1
#define OUTPUT_AS_INPUT_ORDER       2

        if ( optionFlag & OUTPUT_AS_INPUT_ORDER )
        { arrayCopy(); }
        else
        { arrayCopyNRemoveDuplicate(); }

        if ( size == 0 ) { return; }

        if ( alignmentType == 1 )
        {
            if ( optionFlag & DO_NOT_OUTPUT_HALF_ALIGNED )
            {
                filterInvalid();
            }
        }
        else
        {
            /* alignmentType > 1 */
            filterBest();

            if ( alignmentType == 2 )
            {
                // do nothing
            }
            else if ( alignmentType == 3 )
            {
                if ( size != 1 )
                { size = 0; }
            }
            else if ( alignmentType == 4 )
            {
                // elements[0] = elements[rand() % size];
                size = 1;
            }
        }
    }
};


struct AlgnmtFlags
{
    uint MASK[32];
    uint size;
    uint * flags;
    pthread_mutex_t occupy_mutex;

    AlgnmtFlags ( uint range = 16777216 );
    void clear();
    inline void increaseSize ( uint newSize );
    inline void reserveSize ( AlgnmtFlags * algnFlags );

    // set(int readID) uses mutex and is thread-safe
    void set ( int readID );

    void get ( vector<int> * diff );
    void getXOR ( AlgnmtFlags * algnFlags, vector<int> * diff );

    void XOR ( AlgnmtFlags * algnFlags );
    void AND ( AlgnmtFlags * algnFlags );

    ~AlgnmtFlags();
};

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Utility class declaration//////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////// class TimeRecorder ///////////////////////////////////
template <class T>
class TimeRecorder
{
    private:
        double referenceTime;
        struct TimeLine
        {
            int state;
            int length;
            TimeLine();
            vector<double> startTime;
            vector<double> endTime;
        };
        map<string, TimeLine> recorder;
    private:
        double setStartTime();
        double getElapsedTime ( double startTime );
    public:
        TimeRecorder();
        void reset();
        void appendStart ( string label = "default", string type = "CPU" );
        void appendEnd ( string label = "default", string type = "CPU" );
        double getTotalTime ( string label = "default" );
        void printTimeLine ( string label = "default", FILE * output = stdout );
};

////////////////////////////////// class MultiThreadDelegator ////////////////////////////////
template <class T>
void * MultiThreadDelegator_call_func ( void * classRef );

template <class ArgType>
class MultiThreadDelegator
{
    private:
        typedef stack<int> ThreadPool;
        int numOfThreads;
        ArgType * argList;
        pthread_t * threadHandles;
        sem_t * threadSems;
        ThreadPool * threadPool;
        int * threadFlags;
        sem_t availableSem, finishSem;
        pthread_mutex_t occupy_mutex, finish_mutex;
        int createThreadCnt;
        int finishFlag;
        void ( *funcWrapper ) ( int, ArgType & );
        void ( *threadInit ) ( void );
        void ( *threadFinalize ) ( void );
    private:
        void clear();
        int allocThread();
        void releaseThread ( int threadId );
        int checkFinish ( int threadId );
    protected:
        void thread_run_func();
        template <class T>
        friend void * MultiThreadDelegator_call_func ( void * classRef );
    public:
        MultiThreadDelegator();
        ~MultiThreadDelegator();

        void init ( int numOfThreads, void ( *threadRunFunc ) ( int, ArgType & ),
                    void ( *threadInitFunc ) ( void ) = NULL, void ( *threadFinalizeFunc ) ( void ) = NULL );
        int schedule ( ArgType & arg );
        void finalize();
};

////////////////////////////////// struct MyCigarStringEncoder ////////////////////////////////
template <class T>
struct CigarStringEncoder
{
    char cigarType[4096];
    int cigarCnt[4096];
    char lastType;
    int lastCnt;
    int cIndex;
    //Encoded results
    char * cigarString;
    int charCount[128], gapPenalty;

    CigarStringEncoder();
    void append ( char type, int cnt );
    void encodeCigarString ( int GapOpenScore, int GapExtendScore );
};

////////////////////////////////// struct CigarStringDecoder ////////////////////////////////
template <class T>
struct CigarStringDecoder
{
    char * cigarString;
    int index;

    CigarStringDecoder ( char * cigarString );
    bool isEmpty();
    void decodeNext ( char & cigarType, int & cigarCnt );
    int totalCount ( char c1, char c2 = 'N' );
};



//////////////////////////////////// struct CigarStringEncoder //////////////////////////////
template <class T>
CigarStringEncoder<T>::CigarStringEncoder()
{
    lastType = 'N';
    lastCnt = 0;
    cIndex = 0;
}
template <class T>
void CigarStringEncoder<T>::append ( char type, int cnt )
{
    if ( lastType == type )
    { lastCnt += cnt; }
    else
    {
        cigarType[cIndex] = lastType;
        cigarCnt[cIndex] = lastCnt;
        ++cIndex;
        lastType = type;
        lastCnt = cnt;
    }
}
template <class T>
void CigarStringEncoder<T>::encodeCigarString ( int GapOpenScore, int GapExtendScore )
{
    memset ( charCount, 0, 128 * sizeof ( int ) );
    gapPenalty = 0;
    char buf[1024];
    int len = 0;
    append ( 'N', 0 );

    for ( int i = cIndex - 1; i > 0; i-- )
    {
        char type = cigarType[i];
        int cnt = cigarCnt[i];

        if ( cnt > 0 )
        {
            len += sprintf ( buf + len, "%u%c", cnt, type );
            charCount[type] += cnt;

            if ( type == 'I' || type == 'D' )
            {
                gapPenalty += GapOpenScore + ( cnt - 1 ) * GapExtendScore;
            }
        }
    }

    cigarString = ( char * ) malloc ( len + 1 );
    memcpy ( cigarString, buf, len );
    cigarString[len] = 0;
}

//////////////////////////////////// struct CigarStringDecoder //////////////////////////////
template <class T>
CigarStringDecoder<T>::CigarStringDecoder ( char * cigarString )
{
    this->cigarString = cigarString;
    index = 0;
}
template <class T>
bool CigarStringDecoder<T>::isEmpty()
{
    return ( cigarString[index] == 0 );
}
template <class T>
void CigarStringDecoder<T>::decodeNext ( char & cigarType, int & cigarCnt )
{
    cigarCnt = 0;
    char c = cigarString[index];

    while ( c >= '0' && c <= '9' )
    {
        cigarCnt = cigarCnt * 10 + c - '0';
        c = cigarString[++index];
    }

    cigarType = cigarString[index++];
}
template <class T>
int CigarStringDecoder<T>::totalCount ( char c1, char c2 )
{
    char type;
    int cnt;
    uint total = 0;

    while ( !isEmpty() )
    {
        decodeNext ( type, cnt );

        if ( type == c1 || type == c2 )
        {
            total += cnt;
        }
    }

    index = 0;
    return total;
}

//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Utility class definition ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////// class TimeRecorder /////////////////////////////////
template <class T>
TimeRecorder<T>::TimeLine::TimeLine()
{
    state = 0;
    length = 0;
}
template <class T>
TimeRecorder<T>::TimeRecorder()
{
    reset();
}
template <class T>
double TimeRecorder<T>::setStartTime()
{
    struct timeval tp;
    gettimeofday ( &tp, NULL );
    return ( double ) tp.tv_sec + ( double ) tp.tv_usec / ( double ) 1000000;
}
template <class T>
double TimeRecorder<T>::getElapsedTime ( double startTime )
{
    struct timeval tp;
    gettimeofday ( &tp, NULL );
    return ( double ) tp.tv_sec + ( double ) tp.tv_usec / ( double ) 1000000 - startTime;
}
template <class T>
void TimeRecorder<T>::reset()
{
    referenceTime = setStartTime();
}
template <class T>
void TimeRecorder<T>::appendStart ( string label, string type )
{
    TimeLine & timeLine = recorder[label];

    if ( timeLine.state == 1 )
    { return; }

    if ( type == "GPU" )
    {
        cudaEventCreate ( &timeLine.start );
        cudaEventCreate ( &timeLine.stop );
        cudaEventRecord ( timeLine.start, 0 );
    }

    timeLine.startTime.push_back ( getElapsedTime ( referenceTime ) );
    timeLine.state = 1;
}
template <class T>
void TimeRecorder<T>::appendEnd ( string label, string type )
{
    TimeLine & timeLine = recorder[label];

    if ( timeLine.state == 0 )
    { return; }

    if ( type == "GPU" )
    {
        float elapsedTime;
        cudaEventRecord ( timeLine.stop, 0 );
        cudaEventSynchronize ( timeLine.stop );
        cudaEventElapsedTime ( &elapsedTime, timeLine.start, timeLine.stop );
        cudaEventDestroy ( timeLine.start );
        cudaEventDestroy ( timeLine.stop );
        timeLine.endTime.push_back ( timeLine.startTime[timeLine.length] + elapsedTime / 1000.0 );
    }
    else
    { timeLine.endTime.push_back ( getElapsedTime ( referenceTime ) ); }

    timeLine.length += 1;
    timeLine.state = 0;
}
template <class T>
double TimeRecorder<T>::getTotalTime ( string label )
{
    TimeLine & timeLine = recorder[label];
    double totalTime = 0;

    for ( int i = 0; i < timeLine.length; i++ )
    {
        totalTime += timeLine.endTime[i] - timeLine.startTime[i];
    }

    return totalTime;
}
template <class T>
void TimeRecorder<T>::printTimeLine ( string label, FILE * output )
{
    TimeLine & timeLine = recorder[label];
    fprintf ( output, "Label: %s  --  Total %d records.\n", label.c_str(), timeLine.length );

    for ( uint i = 0; i < timeLine.length; i++ )
    {
        fprintf ( output, "[%u] %lf -> %lf\n", i, timeLine.startTime[i], timeLine.endTime[i] );
    }
}

//////////////////////////////////// class MultiThreadDelegator //////////////////////////////
template <class T>
void * MultiThreadDelegator_call_func ( void * classRef )
{
    MultiThreadDelegator<T> * delegator = ( MultiThreadDelegator<T> * ) classRef;
    delegator -> thread_run_func();
    return NULL;
}
template <class ArgType>
MultiThreadDelegator<ArgType>::MultiThreadDelegator()
{
    threadHandles = NULL;
}
template <class ArgType>
MultiThreadDelegator<ArgType>::~MultiThreadDelegator()
{
    clear();
}
template <class ArgType>
void MultiThreadDelegator<ArgType>::clear()
{
    if ( threadHandles != NULL )
    {
        delete[] threadHandles;
        delete[] threadSems;
        delete threadPool;
        delete[] threadFlags;
        delete[] argList;
    }

    threadHandles = NULL;
}
template <class ArgType>
int MultiThreadDelegator<ArgType>::allocThread()
{
    pthread_mutex_lock ( &occupy_mutex );
    int threadId = threadPool->top();
    threadPool->pop();
    threadFlags[threadId] = 1;
    pthread_mutex_unlock ( &occupy_mutex );
    return threadId;
}
template <class ArgType>
void MultiThreadDelegator<ArgType>::releaseThread ( int threadId )
{
    pthread_mutex_lock ( &occupy_mutex );
    threadPool->push ( threadId );
    threadFlags[threadId] = 0;
    pthread_mutex_unlock ( &occupy_mutex );
}
template <class ArgType>
void MultiThreadDelegator<ArgType>::init ( int num, void ( *inFuncWrapper ) ( int, ArgType & ),
        void ( *inThreadInit ) ( void ),
        void ( *inThreadFinalize ) ( void ) )
{
    clear();
    numOfThreads = num;
    threadHandles = new pthread_t[numOfThreads];
    threadSems = new sem_t[numOfThreads];
    threadPool = new ThreadPool;
    threadFlags = new int[numOfThreads];
    argList = new ArgType[numOfThreads];
    funcWrapper = inFuncWrapper;
    threadInit = inThreadInit;
    threadFinalize = inThreadFinalize;
    finishFlag = 0;

    for ( int i = 0; i < numOfThreads; i++ )
    {
        sem_init ( threadSems + i, 0, 0 );
        threadFlags[i] = 0; // empty
    }

    sem_init ( &availableSem, 0, 0 );
    sem_init ( &finishSem, 0, 0 );
    pthread_mutex_init ( &occupy_mutex, NULL );
    pthread_mutex_init ( &finish_mutex, NULL );
    createThreadCnt = 0;

    for ( int i = 0; i < numOfThreads; i++ )
    {
        pthread_create ( threadHandles + i, NULL, MultiThreadDelegator_call_func<ArgType>, this );
    }
}
template <class ArgType>
int MultiThreadDelegator<ArgType>::checkFinish ( int threadId )
{
    int flag;
    pthread_mutex_lock ( &occupy_mutex );
    flag = finishFlag && ( threadFlags[threadId] == 0 );
    pthread_mutex_unlock ( &occupy_mutex );
    return flag;
}
template <class ArgType>
void MultiThreadDelegator<ArgType>::thread_run_func()
{
    pthread_mutex_lock ( &occupy_mutex );
    int threadId = createThreadCnt++;
    pthread_mutex_unlock ( &occupy_mutex );

    if ( threadInit != NULL )
    { threadInit(); }

    releaseThread ( threadId );
    sem_post ( &availableSem );
    sem_t * pThreadSem = threadSems + threadId;
    sem_wait ( pThreadSem );

    while ( checkFinish ( threadId ) == 0 )
    {
        funcWrapper ( threadId, argList[threadId] );
        releaseThread ( threadId );
        sem_post ( &availableSem );
        sem_wait ( pThreadSem );
    }

    if ( threadFinalize != NULL )
    { threadFinalize(); }
}
template <class ArgType>
int MultiThreadDelegator<ArgType>::schedule ( ArgType & arg )
{
    sem_wait ( &availableSem );
    int threadId = allocThread();
    argList[threadId] = arg;
    sem_post ( threadSems + threadId );
    return threadId;
}
template <class ArgType>
void MultiThreadDelegator<ArgType>::finalize()
{
    pthread_mutex_lock ( &occupy_mutex );
    finishFlag = 1;
    pthread_mutex_unlock ( &occupy_mutex );

    for ( int i = 0; i < numOfThreads; i++ )
    { sem_post ( threadSems + i ); }

    for ( int i = 0; i < numOfThreads; i++ )
    { pthread_join ( threadHandles[i], NULL ); }

    clear();
}


/////////////////////////////////////// single-dp space //////////////////////////////////////

namespace SingleDP_Space
{

struct SeedPos
{
    uint strand : 2;  // value 1 = '+'
    unsigned long long pos;
    uint readID;
    uint seedAlignmentLength ;  // for statistical usage
    uint peLeftAnchor : 31;
    uint peRightAnchor : 31;
};

typedef struct SeedPos CandidateInfo;

struct CandidateStream
{
    vector<CandidateInfo> data;
    pthread_mutex_t occupy_mutex;
    CandidateStream();
    void append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags );
};

class SingleEndSeedingEngine
{
    protected:
#define DPS_DIVIDE_GAP 5

        struct SingleEndSeedingBatch
        {
            SingleEndSeedingEngine * home;
            
            /* SOAP3 seeding parameters
             * *********************/
            uint * queries;
            uint * queryLengths;
            int wordPerQuery;
            int * seedPositions;

            /* SOAP3 seeding, Input
             * *********************/
            int numOfCPUForSeeding, maxHitNum;
            uint batchSize, maxSeedLength, wordPerSeed, numQueries;
            uint * readIDs, *seeds, *lengths, *offsets;

            /* SOAP3 seeding, Output
             * *********************/
            uint64 numOfAnswer;
            uint numOfAlignedRead;

            SingleEndSeedingBatch (
                uint batchSize, DPParameters * dpPara,
                uint * queries, uint * queryLengths, uint inputMaxReadLength, SingleEndSeedingEngine* home
            );
            ~SingleEndSeedingBatch();

            void clear();
            inline void pack ( uint readID, int off, int seedLength );
            int packSeeds ( uint readID, int stage );

            SeedPos * decodePositions ( BWT * bwt );
            vector<CandidateInfo> * decodeMergePositions ( BWT * bwt );
        };

        struct SingleEndSeedingThreadContext
        {
            SingleEndSeedingBatch * batch;
            sem_t ACKSem;
            sem_t GPUFinishSem;

            void init ( SingleEndSeedingBatch * batch );
            void freeMemory();
        };

        static SingleEndSeedingEngine * engine;


        QueryIDStream  *  queryIDStream;
        DPParameters   *  dpPara;
        uint       *      queries;
        uint       *      queryLengths;
        int inputMaxReadLength;
        Soap3Index    *   index;
        CandidateStream * canStream;
        QueryIDStream  *  unseededIDStream;
        SOAP3Wrapper<void>  * soap3Wrapper;

        
        AlgnmtFlags * inputFlags, *alignFlags;
        SingleEndSeedingBatch      *      seedingSwapBatch;
        SingleEndSeedingThreadContext  *  seedingThreadContext;
        MultiThreadDelegator<void *>    seedingCPUThreadDelegator;
        MultiThreadDelegator<int *>     seedingGPUThreadDelegator;

    public:
        SingleEndSeedingEngine();
        vector<CandidateInfo> * singleMerge ( SeedPos * readPos, QueryIDStream * eliminatedReadIDs = NULL );
};



void algnmtCPUThread ( int threadId, void *& empty );
void algnmtGPUThreadInit();
void algnmtGPUThread ( int threadId, int *& pCallThreadId );
void algnmtGPUThreadFinalize();

void DPSOutputThread ( int threadId, int *& pCallThreadId );
void DPSOutputThreadFinalize();
void DPSOutputIndelRealignmentThread ( int threadId, int *& pCallThreadId );
void DPSOutputIndelRealignmentThreadFinalize();


// fetching SingleDPResult
SingleAlgnmtResult * fetchNextSingleDPResult ( SingleDP_Space::AlgnmtResultStream * singleDPResultStream, 
                                               uint &singleDPResultStreamIdx, uint &singleDPResultBatchIdx );

// Output Alignment reads cannot be be paired
// including output unaligned reads
// both of singleDPResultStream and unalignedSingleReads are sorted by readID
void DPSOutputUnpairedAlignment ( AlgnmtResultStream * singleDPResultStream, UnalignedSingles * unalignedSingleReads,
                                 uint * queries, uint * upkdReadLengths, uint * origReadIDs, 
                                 char ** queryNames, char ** queryComments, char * upkdQualities, int inputMaxReadLength,
                                 uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
                                 Soap3Index * index,
                                 int alignmentType,
                                 uint         &        dpSAlignedRead,
                                 uint         &        dpSAlignment   );


class SingleEndAlignmentEngine
{
    protected:

        struct SingleEndAlgnBatch
        {
            Soap3Index * index;
            uint * packedDNA;
            uint * queries;
            uint inputMaxReadLength;
            uint * upkdLengths;
            DPInfoForReads * dpInfoForReads;
            int wordPerOldQuery, wordPerQuery, wordPerDNA;
            unsigned long long fullDNALength;

            CandidateInfo * canInfos;
            int softClipLeft, softClipRight;
            int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
            uint * packedDNASeq, *packedReadSeq;
            uint * DNALengths, *lengths, *hitLocs;
            uchar * pattern;
            int * scores, *cutoffThresholds;
            uint * softClipLtSizes, *softClipRtSizes;
            uint * peLeftAnchorLocs, * peRightAnchorLocs;
            uint * maxScoreCounts;

            int batchSize;
            uint numOfThreads;

            // Array index of Modified Reference for Indel Realignment
            unsigned int * modRefIDs;

            SingleEndAlgnBatch (
                int batchSize, DPParameters * dpPara,
                int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
                Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths,
                DPInfoForReads * dpInfoForReads );
            ~SingleEndAlgnBatch();

            void clear();
            int pack ( CandidateInfo & canInfo );
            inline void packRead ( uint * packedSeq, uint threadId, uint readID, uint length, int strand );
            inline void repackDNA ( uint * packedSeq, uint threadId, uint * seq, unsigned long long start, uint length );
            inline void repackModifiedDNA ( uint * packedSeq, uint threadId, uint * modSeq, uint length );
        };

        struct SingleEndAlgnThreadContext
        {
	    SemiGlobalAligner semiGlobalAligner;

            SingleEndAlgnBatch * batch;
            SingleDPResultBatch * resultBatch;
            int batchID;
            sem_t ACKSem;
            sem_t GPUFinishSem;
            sem_t outputACKSem;
            void init ( SingleEndAlgnBatch * batch );
            void freeMemory();
        };

        

        friend void algnmtCPUThread ( int threadId, void *& empty );
        friend void algnmtGPUThreadInit();
        friend void algnmtGPUThread ( int threadId, int *& pCallThreadId );
        friend void algnmtGPUThreadFinalize();

        friend void DPSOutputThread ( int threadId, int *& pCallThreadId );
        friend void DPSOutputThreadFinalize();
        friend void DPSOutputIndelRealignmentThread ( int threadId, int *& pCallThreadId );
        friend void DPSOutputIndelRealignmentThreadFinalize();

        static SingleEndAlignmentEngine * engine;

        int algnBatchCount;

        //  for output
        int dpSAlignedRead, dpSAlignment;
        int lastReadID;
        AlgnmtFlags * inputFlags, *alignFlags;
        AlgnmtResultStream * resultStream;
        OutputBuffer<SingleAlgnmtResult> * outputBuf;

        CandidateStream * canStream;
        QueryIDStream * unalignedIDStream;

        SingleEndAlgnBatch * algnSwapBatch;
        SingleEndAlgnThreadContext * algnThreadContext;

        int inputMaxReadLength;
        int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
        int DPS_ALGN_NUM_OF_BLOCKS;

        //BWT structure
        Soap3Index * index;
        uint * queries;
        char ** queryNames;
        char ** queryComments;
        char *upkdQualities;
        uint * origReadIDs, *upkdReadLengths;
        DPInfoForReads * dpInfoForReads;

        //GPU data structure
        SemiGlobalAligner semiGlobalAligner;

        //parameters
        int alignmentType;
        DPParameters * dpPara;
        //for output
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;

        MultiThreadDelegator<void *> algnmtCPUThreadDelegator;
        MultiThreadDelegator<int *>  algnmtGPUThreadDelegator;
        MultiThreadDelegator<int *>  outputThreadDelegator;

        void performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment );
        void performAlignmentForPreparedRegion ( uint & numDPAlignedRead, uint & numDPAlignment );

    public:
        static void performAlignment (
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
        );

        static void performAlignmentForPreparedRegion (
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
        );
};

void DPSOutputUnalignedReads (
    QueryIDStream * unalignedIDStream,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr
);
}

// TODO
/////////////////////////////////////// default-dp space //////////////////////////////////////

namespace DP_Space
{
typedef vector<AlgnmtDPResult> DPResultBatch;

struct HalfEndOccStream
{
    ReadInputForDPArrays * data;
    SingleDP_Space::AlgnmtResultStream * singleDPData;
    uint lastFetchReadID;   // selection count for Default DP 
    int lastFetchReadIDCount; // selection count for Default DP 
    int counter; // selection count for Default DP
    int singleDPDataArrayIndex;
    int iter_singleDPData;
    
    BWT * bwt;

    int arrayIndex;      
    ReadInputForDP      *     iter_readInput;
    SRAOccurrence      *      iter_occ;
    SRAOccurrence      *      end_occ;
    PESRAAlignmentResult   *  iter_sa;
    PESRAAlignmentResult   *  end_sa;

    int nextSAIndex;     

    HalfEndOccStream ( ReadInputForDPArrays * input, BWT * bwt );
    HalfEndOccStream ( SingleDP_Space::AlgnmtResultStream * input );
    int fetchNextOcc ( SRAOccurrence & occ );
    int fetchNextSingleAlgnResult( SingleAlgnmtResult & singleAlgnResult );
};

struct CandidateInfo
{
    SRAOccurrence refer;
    int leftOrRight;
};

struct SingleEndDPCandidateInfo
{
    SingleAlgnmtResult refer;
    int leftOrRight;
};

void algnmtCPUThread ( int threadId, void *& empty );
void algnmtGPUThreadInit();
void algnmtGPUThread ( int threadId, int *& pCallThreadId );
void algnmtGPUThreadFinalize();

void DPOutputThread ( int threadId, int *& pCallThreadId );
void DPOutputThreadFinalize();

class HalfEndAlignmentEngine
{
    protected:

        struct HalfEndAlgnBatch
        {        
            Soap3Index * index;
            uint inputMaxReadLength;
            unsigned long long fullDNALength;
            uint * queries;
            uint * upkdReadLengths;
            int wordPerOldQuery, wordPerQuery, wordPerDNA;

            int peStrandLeftLeg, peStrandRightLeg;
            int insert_high, insert_low;
            int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
            int softClipLeft, softClipRight;
            int isDoubleStrand;

            CandidateInfo * canInfo;
            SingleEndDPCandidateInfo * singleEndDPCandidateInfo;

            uint * packedDNASequence, *packedReadSequence;
            uint * DNALengths, *lengths, *hitLocs;
            unsigned long long *startLocs;
            uint * peLeftAnchorLocs, *peRightAnchorLocs;
            uint * softClipLtSizes, *softClipRtSizes;
            int * scores, *cutoffThresholds;
            uchar * pattern;
            uint * maxScoreCounts;

            int batchSize;
            uint numOfThreads;

            HalfEndAlgnBatch (
                int batchSize, DPParameters * dpPara,
                int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
                int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
                Soap3Index * index, uint * queries, int inputMaxReadLength, uint * upkdReadLengths, int dPFromSingleEndDPResult
            );
            ~HalfEndAlgnBatch();
            void clear();

            int pack ( SingleAlgnmtResult & curSingleAlgnmtResult );
            inline void packRead ( uint * packedSeq, uint threadId, uint readID, uint length, int strand );
            inline void repackDNA ( uint * packedSeq, uint threadId, uint * seq, unsigned long long start, uint length );
        };

        struct HalfEndAlgnThreadContext
        {
            HalfEndAlgnBatch * batch;
            DPResultBatch resultBatch;
            int batchID;
            sem_t dispatchACKSem;
            sem_t GPUFinishSem;
            sem_t outputACKSem;

        SemiGlobalAligner semiGlobalAligner;
            void init ( HalfEndAlgnBatch * batch );
            void freeMemory();
        };

        struct AlgnmtResultStream
        {
            uint numOut;
            pthread_mutex_t occupy_mutex;
            AlgnmtResultStream();
            ~AlgnmtResultStream();
        };

        friend void algnmtCPUThread ( int threadId, void *& empty );
        friend void algnmtGPUThreadInit();
        friend void algnmtGPUThread ( int threadId, int *& pCallThreadId );
        friend void algnmtGPUThreadFinalize();

        friend void DPOutputThread ( int threadId, int *& pCallThreadId );
        friend void DPOutputThreadFinalize();

        static HalfEndAlignmentEngine * engine;

        int algnBatchCount;

        //  for output
        int dpAlignedRead, dpAlignment;
        AlgnmtFlags * inputFlags, *alignFlags;
        AlgnmtResultStream * resultStream;

        HalfEndOccStream * canStream;
        QueryIDStream * unalignedIDStream;

        HalfEndAlgnBatch * algnSwapBatch;
        HalfEndAlgnThreadContext * algnThreadContext;

        int inputMaxReadLength;
        int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
        int insert_high, insert_low;
        int peStrandLeftLeg, peStrandRightLeg;
        int DP_ALGN_NUM_OF_BLOCKS;

        //BWT structure
        Soap3Index * index;
        uint * queries;
        char ** queryNames;
        char ** queryComments;
        char * upkdQualities;
        uint * origReadIDs, *upkdReadLengths;
        uint * occValue, *revOccValue;
        DPInfoForReads * dpInfoForReads;

        //GPU data structure
        SemiGlobalAligner semiGlobalAligner;

        //parameters
        int alignmentType;
        DPParameters * dpPara;
        //for output
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;
        samfile_t ** currSamOutputFilePtr;

        uint * numDPAlignedReads;
        uint * numDPAlignments;

        MultiThreadDelegator<void *> algnmtCPUThreadDelegator;
        MultiThreadDelegator<int *>  algnmtGPUThreadDelegator;
        MultiThreadDelegator<int *>  outputThreadDelegator;

        void performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment );

    public:
        static void performAlignment (
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
            uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, samfile_t ** currSamOutputFilePtr,
            /* output */
            QueryIDStream    *    unalignedIDStream,
            uint         &        numDPAlignedRead,
            uint         &        numDPAlignment
        );
};
}


// TODO
//////////////////////////////////////// deep-dp space ///////////////////////////////////////

namespace DeepDP_Space
{
typedef vector<DeepDPAlignResult> DP2ResultBatch;

struct CandidateInfo
{
    //  readID of left read
    uint readIDLeft;
    unsigned long long pos[2];
    uint dnaLength[2];
    uint peLeftAnchor[2];
    uint peRightAnchor[2];
};

struct CandidateStream
{
    vector<CandidateInfo> data;
    pthread_mutex_t occupy_mutex;
    CandidateStream();
    void append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags );
};

void SeedingCPUThread ( int threadId, void *& empty );
void SeedingGPUThreadInit();
void SeedingGPUThread ( int threadId, int *& pCallThreadId );
void SeedingGPUThreadFinalize();

class PairEndSeedingEngine
{
    protected:
        struct PairEndSeedingBatch
        {
            int numberOfShortcut;
            int lazyCnt;
            PairEndSeedingEngine * home;
#define DP2_DIVIDE_GAP 5

            /* SOAP3 seeding parameters
            * *********************/
            BWT * bwt;
            uint * queries, *readLengths;
            int insert_high, insert_low;
            int peStrandLeftLeg, peStrandRightLeg;
            int wordPerQuery;
            int numOfCPUForSeeding;
            int * seedPositions;
            int maxHitNum[2];

            /* SOAP3 seeding, Input
            * *********************/
            uint batchSize, maxSeedLength, wordPerSeed;
            uint numQueries[2];
            uint * readIDs[2], *lengths[2], *minLengths[2], *offsets[2];
            uint * seeds[2];
            uint numNBMismatch;
            int restrictedNBM;

            /* SOAP3 seeding, Output
            * *********************/
            uint64 numOfAnswer[2];
            uint numOfAlignedRead[2];

            /* internal */
            vector<SeedPos> inPosArr[2];
            int lastPairID;

            PairEndSeedingBatch (
                uint batchSize, DPParameters * dpPara,
                uint * queries, uint * readLengths, uint inputMaxReadLength,
                int insert_high, int insert_low,
                int peStrandLeftLeg, int peStrandRightLeg,
                BWT * bwt, PairEndSeedingEngine * home
            );
            ~PairEndSeedingBatch();
            void clear();

            inline void pack ( uint evenReadID, uint readID, int off, int seedLength, int varMinLength, int readOrMate );
            inline void packForMmpSeeding ( QueryIDStream * queryIDStream, int readOrMate );
            int packSeeds ( uint evenReadID, int stage );
            int packSeeds ( SRAOccurrence & occ, int stage );
            inline int packSeedsOneSide ( uint evenReadID, int readOrMate, int stage );
            int packSeedsOneSide ( SRAOccurrence & occ, int stage );

            uint findRevStart ( SeedPos * arr, uint len );
            void pairEndMerge (
                vector<CandidateInfo> * &pairEndPos,
                SeedPos * readPos, SeedPos * matePos,
                int leftReadOrMate
            );
            uint64_t mmpSeeding( int words_per_query, SeedPos** readPos, SeedPos** matePos, int n_threads, SeedPool *seedPool, MmpProperties mmpPara);
            vector<CandidateInfo> * mergeAndPairPairedEnd ( SeedPos * readPos, SeedPos * matePos, int readPosLen, int matePosLen );
        };

        struct PairEndSeedingThreadContext
        {
            PairEndSeedingBatch * batch;
            sem_t ACKSem;
            sem_t GPUFinishSem;
            void init ( PairEndSeedingBatch * batch );
            void freeMemory();
        };

        static PairEndSeedingEngine * engine;
        QueryIDStream  *  queryIDStream;
        DPParameters   *  dpPara;
        uint       *      queries;
        uint       *      queryLengths;
        int inputMaxReadLength;
        Soap3Index    *   index;
        CandidateStream * canStream;
        QueryIDStream  *  tooManyHitIDStream;
        QueryIDStream  *  unseededIDStream;
        SOAP3Wrapper<void>  * soap3Wrapper;

        DP_Space::HalfEndOccStream * halfEndOccStream;

        int insert_high, insert_low;
        int peStrandLeftLeg, peStrandRightLeg;

        int seedingStage;
        SeedPool* seedPool;

        AlgnmtFlags * inputFlags, *alignFlags, *tooManyHitFlags;
        PairEndSeedingBatch     *     seedingSwapBatch;
        PairEndSeedingThreadContext * seedingThreadContext;
        MultiThreadDelegator<void *>    seedingCPUThreadDelegator;
        MultiThreadDelegator<int *>     seedingGPUThreadDelegator;

        PairEndSeedingEngine();
        friend void SeedingCPUThread ( int threadId, void *& empty );
        friend void SeedingGPUThreadInit();
        friend void SeedingGPUThread ( int threadId, int *& pCallThreadId );
        friend void SeedingGPUThreadFinalize();

        void performMmpSeeding();

    public:
        static void performMmpSeeding (
            /* input */
            QueryIDStream    *    queryIDStream,
            DPParameters     *    dpPara,
            uint * queries, uint * queryLengths, int inputMaxReadLength,
            int insert_high, int insert_low,
            int peStrandLeftLeg, int peStrandRightLeg,
            /* soap3 seeding related */
            SOAP3Wrapper<void>  * soap3Wrapper,
            Soap3Index      *      index, SeedPool * seedPool,
            /* output */
            CandidateStream   *   canStream,
            QueryIDStream    *    unseededIDStream
        ) ;
};

void DP2CPUAlgnThread ( int threadId, void *& empty );
void DP2GPUAlgnThreadInit();
void DP2GPUAlgnThread ( int threadId, int *& pCallThreadId );
void DP2GPUAlgnThreadFinalize();

void DP2OutputThread ( int threadId, int *& pCallThreadId );
void DP2OutputThreadFinalize();

void DP2OutputThread2 ( int threadId, int *& pCallThreadId );
void DP2OutputThread2Finalize();

inline void reverseCigarString ( char * cigar );
inline int cigarAlignmentLength_Extension ( char * cigar );

class PairEndAlignmentEngine
{
    protected:
        struct PairEndAlgnBatch
        {
            int batchSize;

            Soap3Index * index;
            uint * packedDNA;
            uint * queries;
            uint inputMaxReadLength;
            uint * upkdLengths;
            DPInfoForReads * dpInfoForReads;

            int peStrandLeftLeg, peStrandRightLeg;
            int softClipLeft, softClipRight;
            int insert_high, insert_low;
            int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
            int wordPerOldQuery, wordPerQuery, wordPerDNA;
            unsigned long long fullDNALength;

            CandidateInfo * canInfos;
            uint * packedDNASeq;
            uint * packedReadSeq;
            uint * DNALengths;
            uint * lengths;
            int * cutoffThresholds;

            uint * hitLocs[2];
            uchar * pattern[2];
            int * scores[2];
            uint * maxScoreCounts[2];

            uint * peLeftAnchorLocs, *peRightAnchorLocs;
            uint * softClipLtSizes, *softClipRtSizes;

            int leftOrRight;
            uint numOfThreads;

            PairEndAlgnBatch (
                int batchSize, DPParameters * dpPara,
                int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
                int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
                Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths,
                DPInfoForReads * dpInfoForReads
            );
            ~PairEndAlgnBatch();
            void clear();

            int packLeft ( CandidateInfo & canInfo );
            void packRight();

            inline void packRead ( uint * packedSeq, uint threadId, uint readID, uint length, int strand );
            inline void repackDNA ( uint * packedSeq, uint threadId, uint * seq, unsigned long long start, uint length );
            inline void repackDNA_Reverse_Complement ( uint * packedSeq, uint threadId, uint * seq, unsigned long long start, uint length );
        };

        struct PairEndAlgnThreadContext
        {
            PairEndAlgnBatch * batch;
            DP2ResultBatch * resultBatch;
            AlgnmtFlags * alignFlag;
            int lastReadID;
            int dp2AlignedRead, dp2Alignment;
            int batchID;
	        SemiGlobalAligner semiGlobalAligner;
            sem_t ACKSem;
            sem_t GPUFinishSem;
            sem_t outputACKSem;

            void init ( PairEndAlgnBatch * batch );
            void freeMemory();
        };

        struct AlgnmtResultStream
        {
            vector<DP2ResultBatch *> dp2Result;
            uint numOut;
            pthread_mutex_t occupy_mutex;
            AlgnmtResultStream();
            ~AlgnmtResultStream();
        };

        PairEndAlignmentEngine();

        friend void DP2CPUAlgnThread ( int threadId, void *& empty );
        friend void DP2GPUAlgnThreadInit();
        friend void DP2GPUAlgnThread ( int threadId, int *& pCallThreadId );
        friend void DP2GPUAlgnThreadFinalize();

        friend void DP2OutputThread ( int threadId, int *& pCallThreadId );
        friend void DP2OutputThreadFinalize();
        friend void DP2OutputThread2 ( int threadId, int *& pCallThreadId );
        friend void DP2OutputThread2Finalize();

        static PairEndAlignmentEngine * engine;

        int algnBatchCount;

        SeedPool * seedPool;

        //  for output
        int dp2AlignedRead, dp2Alignment;
        int lastReadID;
        AlgnmtFlags * inputFlags, *alignFlags;
        AlgnmtResultStream * resultStream;
        AlgnmtResultStream ** resultStreams;
        OutputBuffer<DeepDPAlignResult> * outputBuf;
        OutputBuffer<DeepDPAlignResult> ** outputBufs;
        char *** resultOutputStrings;
        int * numberOfResultOutputString;
        char ** resultOutputStringBuffer;
        int resultOutputStringBufferIter;
        int resultOutputStringBufferSize;

        CandidateStream * canStream;
        QueryIDStream * unalignedIDStream;

        PairEndAlgnBatch * algnSwapBatch;
        PairEndAlgnThreadContext * algnThreadContext;

        int inputMaxReadLength;
        int insert_high, insert_low;
        int peStrandLeftLeg, peStrandRightLeg;
        int maxReadLength, maxDNALength, maxDPTableLength, patternLength;
        int DP2_ALGN_NUM_OF_BLOCKS;

        //BWT structure
        Soap3Index * index;
        uint * queries;
        char ** queryNames;
        char ** queryComments;
        char * upkdQualities;
        uint * origReadIDs, *upkdReadLengths;   
        DPInfoForReads * dpInfoForReads;

        //GPU data structure
        SemiGlobalAligner semiGlobalAligner;

        //parameters
        bool shortDnaLength;
        int alignmentType;
        DPParameters * dpPara;
        //for output
        uint accumReadNum;
        int outputFormat;
        samfile_t * samOutputDPFilePtr;
        samfile_t ** samOutputDPFilePtrs;
	uint * numDPAlignedReads;
        uint * numDPAlignments;


        MultiThreadDelegator<void *> algnmtCPUThreadDelegator;
        MultiThreadDelegator<int *>  algnmtGPUThreadDelegator;
        MultiThreadDelegator<int *>  outputThreadDelegator;

        void DP2Output ( int *& pCallThreadId );
        void DP2OutputFinalize( int threadId );

        void performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment );

    public:
        static void performAlignment (
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
            uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr, SeedPool * seedPool,
            samfile_t ** samOutputDPFilePtrs,
            /* output */
            QueryIDStream    *    unalignedIDStream,
            uint         &        numDPAlignedRead,
            uint         &        numDPAlignment
        );
};

void DPExtensionPassUnalignedReads (
        QueryIDStream * unalignedIDStream,
        uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
        Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
        char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
        uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr,
        unsigned int maxHitNumForDP, unsigned int maxHitNumForDP2,
        ReadInputForDP * dpInput, ReadInputForDP * dpInputForNewDefault,
        ReadInputForDP * otherSoap3Result, BothUnalignedPairs * bothUnalignedPairs
);

void DP2OutputUnalignedReads (
    QueryIDStream * unalignedIDStream, SeedPool* seedPool,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
    char ** queryNames, char ** queryComments, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat, samfile_t * samOutputDPFilePtr
);
}

#endif

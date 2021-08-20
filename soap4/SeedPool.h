// ------------------------------------------------
//  SeedPool.h
//  -----------------------------------------------
//
//  Use for store overall seed alignments generated
//  in different Seeding Scheme
//
//
#ifndef __SEED_POOL_H__
#define __SEED_POOL_H__

#include "AlgnResult.h"
#include <vector>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Alignment modules //////////////////////////////////////
/////////////////// The following code better be placed in separate files ///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// standard space ///////////////////////////////////////

struct QueryIDStream
{
    vector<int> * data;

    QueryIDStream();
    QueryIDStream ( BothUnalignedPairsArrays * input );
    ~QueryIDStream();
    void append ( QueryIDStream * stream );
    void setBuffer ( vector<int> * input );
};

void mergeTwoSortedQueryIDStream( QueryIDStream * streamA, QueryIDStream * streamB, QueryIDStream * resultStream );

//////////////////////////////////////////////////////////////////////////////////////////////

struct SingleAlgnmtResult;
namespace SingleDP_Space{
typedef vector<SingleAlgnmtResult> SingleDPResultBatch;
struct AlgnmtResultStream
{
    vector<SingleDPResultBatch *> dpSResult;
    uint numOut;
    pthread_mutex_t occupy_mutex;
    AlgnmtResultStream();
    ~AlgnmtResultStream();
};
}

namespace DeepDP_Space {
struct SeedPos
{
    unsigned long long pos;
    uint strand_readID;
    uint paired_seedLength;
};
}

#define SEED_POOL_DYNAMIC_ARRAY_INIT_SIZE 32

template <class Element>
class SeedPoolDynamicArray{
    public:
        Element* element;
        unsigned int size;
        unsigned int availableSize;
        SeedPoolDynamicArray();
        ~SeedPoolDynamicArray();
        void add(Element* e);
        void append(SeedPoolDynamicArray<Element> * array);
        Element get(int index){return element[index];}
        void clear() { size = 0; }
};

typedef struct SeedAlignmentsOfRead{
    // dynamic array of OccRecord Array
    SeedPoolDynamicArray<OccRecord> * occList;
    int numberOfAlignments;
}SeedAlignmentsOfRead;

typedef struct SeedPool{
    // array of seed alignments pointer, readId as index 
    //      seedAlignmentsOfRead[readId] MAY equal to NULL 
    //      if the corresponding read pair passed DP
    SeedAlignmentsOfRead *** seedAlignmentsOfReadArrays;
    unsigned int * flag; // 1 = lazy seeding
    unsigned int numberOfThread;
    unsigned int numberOfReads;
    unsigned int * totalNumberOfAlignmentsPerThread;
}SeedPool;

// construct SeedAlignmentOfRead with corresponding alignment result
SeedAlignmentsOfRead * constructSeedAlignmentsOfRead();

// construct a SeedPool with the totalNumberOfReads
SeedPool * constructSeedPool(unsigned int numberOfReads, unsigned int numberOfThread);

// after adding the alignment to array, return the number of alignments
int addSeedAlignmentsOfRead( SeedAlignmentsOfRead* seedAlignments, OccRecord * occList, unsigned int occListSize);

void addSeedPosToSeedPool( SeedPool * seedPool, DeepDP_Space::SeedPos * seedPositions, int readOrMate, int threadId, size_t size);

void addSingleSeedPosToSeedPool(SeedPool *seedPool, DeepDP_Space::SeedPos &sp, int readOrMate, unsigned int tid);

// add seedAlignmentOfRead into seedPoolThreadArray
void addSeedAlignmentsOfReadToSeedPoolThreadArray ( SeedAlignmentsOfRead ** seedAlignmentsOfReadArray, unsigned int readId, SeedAlignmentsOfRead * seedAlignmentsOfRead );

// add seedAlignmentOfRead into seedPool
void addSeedAlignmentsOfReadToSeedPool( SeedPool * seedPool , unsigned int readId, SeedAlignmentsOfRead * seedAlignmentsOfRead, int threadId );

// free SeedAlignmentsOfRead without freeing up the data
void freeSeedAlignmentsOfRead(SeedAlignmentsOfRead * seedAlignmentsOfRead);

// free up SeedAlignmentOfRead memory with corresponding readId
// the occlist data will be free up
int freeSeedAlignmentsOfReadWithReadId( SeedPool * seedPool, unsigned int readId );

// free up SeedPool
int freeSeedPool( SeedPool * seedPool );
void resetSeedPool(SeedPool *seedPool);

void filterOutUnpairedSingleReads ( UnalignedSingles * unalignedSingleReads, BothUnalignedPairs * bothUnalignedPairs, vector<int> * unalignedReadsFromSingleEndDP );


#endif

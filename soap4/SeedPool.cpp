// ------------------------------------------------
//  SeedPool.cpp
//  -----------------------------------------------
//
//  Use for store overall seed alignments generated
//  in different Seeding Scheme
//

#include "SeedPool.h"

QueryIDStream::QueryIDStream()
{
    data = new vector<int>;
}
QueryIDStream::QueryIDStream ( BothUnalignedPairsArrays * input )
{
    data = new vector<int>;

    for ( int arrIndex = 0; arrIndex < input->arrayNum; arrIndex++ )
    {
        BothUnalignedPairs * array = input->array[arrIndex];

        for ( int pairIndex = 0; pairIndex < array->totalNum; pairIndex++ )
        {
            data->push_back ( array->readIDs[pairIndex] );
        }
    }
}
QueryIDStream::~QueryIDStream()
{
    delete data;
}
void QueryIDStream::append ( QueryIDStream * stream )
{
    data->insert ( data->end(), stream->data->begin(), stream->data->end() );
}
void QueryIDStream::setBuffer ( vector<int> * input )
{
    delete data;
    data = input;
}

void mergeTwoSortedQueryIDStream( QueryIDStream * streamA, QueryIDStream * streamB, QueryIDStream * resultStream )
{
    int i,j;     
    for (i=0,j=0; i<streamA->data->size() && j<streamB->data->size();)
    {
        if ( (*(streamA->data))[i] < (*(streamB->data))[j] )
        {
            resultStream->data->push_back( (*(streamA->data))[i] );
            ++i;
        }
        else
        {
            resultStream->data->push_back( (*(streamB->data))[j] );
            ++j;
        }
    }
    for ( ; i<streamA->data->size() ; ++i )
        { resultStream->data->push_back( (*(streamA->data))[i] ); }
    for ( ; j<streamB->data->size() ; ++j )
        { resultStream->data->push_back( (*(streamB->data))[j] ); }
    
    for ( i=1;i<resultStream->data->size(); ++i )
        if ( (*(resultStream->data))[i-1] >= (*(resultStream->data))[i] ){
            fprintf(stderr, "Problem:%u %u %u\n", i, (*(resultStream->data))[i-1], (*(resultStream->data))[i]);
        }
}


template <class Element>
SeedPoolDynamicArray<Element>::SeedPoolDynamicArray()
{
    size=0;
    availableSize=SEED_POOL_DYNAMIC_ARRAY_INIT_SIZE;
    element = (Element*)malloc(availableSize*sizeof(Element));
}

template <class Element>
SeedPoolDynamicArray<Element>::~SeedPoolDynamicArray()
{
    free(element);
}

template <class Element>
void SeedPoolDynamicArray<Element>::add(Element* e)
{
    if (size+1>=availableSize)
    {
        while (size+1 >= availableSize)
        {
            availableSize *= 2;
        }
        element = (Element*)realloc(element, availableSize*sizeof(Element));
    }
    element[size++] = (*e);
}

template <class Element>
void SeedPoolDynamicArray<Element>::append(SeedPoolDynamicArray<Element> * array)
{
    if (size+array->size>=availableSize)
    {
        while (size+array->size >= availableSize)
        {
            availableSize *= 2;
        }

        element = (Element*)realloc(element, availableSize*sizeof(Element));
    }
    memcpy(element+size,array->element,array->size * sizeof(Element));
    size+=array->size;
}

// construct SeedAlignmentOfRead with corresponding alignment result
SeedAlignmentsOfRead * constructSeedAlignmentsOfRead()
{
    SeedAlignmentsOfRead * seedAlignmentsOfRead = (SeedAlignmentsOfRead *)malloc(sizeof (SeedAlignmentsOfRead ));
    seedAlignmentsOfRead->occList = new SeedPoolDynamicArray<OccRecord>();
    seedAlignmentsOfRead->numberOfAlignments = 0;
    return seedAlignmentsOfRead;
}

int addSeedAlignmentsOfRead( SeedAlignmentsOfRead* seedAlignments, OccRecord * occList, unsigned int occListSize)
{
    for ( int i=0; i<occListSize; ++i )
        { seedAlignments->occList->add(&(occList[i])); }
    seedAlignments->numberOfAlignments += occListSize;
    return seedAlignments->numberOfAlignments;
}

void addSeedPosToSeedPool( SeedPool * seedPool, DeepDP_Space::SeedPos * seeds, int readOrMate, int threadId, size_t size )
{
    for (size_t i = 0; i < size; ++i)
    {
        addSingleSeedPosToSeedPool(seedPool, seeds[i], readOrMate, threadId);
    }
}

SeedPool * constructSeedPool(unsigned int numberOfReads, unsigned int numberOfThread )
{
    SeedPool * seedPool = (SeedPool*) malloc (sizeof(SeedPool));
    seedPool->numberOfReads = numberOfReads;
    seedPool->numberOfThread = numberOfThread;
    seedPool->totalNumberOfAlignmentsPerThread = (uint *)malloc(numberOfThread * sizeof(uint));
    seedPool->flag = (uint* ) malloc ( numberOfReads * sizeof ( uint ) );
    seedPool->seedAlignmentsOfReadArrays = (SeedAlignmentsOfRead***) malloc(numberOfThread * sizeof(SeedAlignmentsOfRead**));
    for (int threadId=0; threadId<numberOfThread; ++threadId)
    {
        seedPool->totalNumberOfAlignmentsPerThread[threadId] = 0;
        seedPool->seedAlignmentsOfReadArrays[threadId] = (SeedAlignmentsOfRead**)malloc(numberOfReads*sizeof(SeedAlignmentsOfRead*));
        memset(seedPool->seedAlignmentsOfReadArrays[threadId], NULL, numberOfReads*sizeof(SeedAlignmentsOfRead*));
    }
    memset( seedPool->flag, 0, numberOfReads * sizeof( uint ) );

    // for (int i = 0; i < numberOfReads; ++i) {
    //     if (seedPool->seedAlignmentsOfReadArrays[0][i] == NULL)
    //     {
    //         SeedAlignmentsOfRead * tmp = constructSeedAlignmentsOfRead();
    //         seedPool->seedAlignmentsOfReadArrays[0][i] = tmp;
    //     }
    // }

    return seedPool;
}

void addSeedAlignmentsOfReadToSeedPoolThreadArray ( SeedAlignmentsOfRead ** seedAlignmentsOfReadArray, unsigned int readId, SeedAlignmentsOfRead * seedAlignmentsOfRead )
{
    if ( seedAlignmentsOfReadArray[readId] == NULL )
    {
        SeedAlignmentsOfRead * tmp = constructSeedAlignmentsOfRead();
        seedAlignmentsOfReadArray[readId] = tmp;
    }
    seedAlignmentsOfReadArray[readId]->occList->append ( seedAlignmentsOfRead->occList );
    seedAlignmentsOfReadArray[readId]->numberOfAlignments += seedAlignmentsOfRead->numberOfAlignments;
}

// add or append seedAlignmentOfRead into seedPool
void addSeedAlignmentsOfReadToSeedPool( SeedPool * seedPool , unsigned int readId, SeedAlignmentsOfRead * seedAlignmentsOfRead, int threadId )
{
    if (seedPool->seedAlignmentsOfReadArrays[threadId][readId] == NULL)
    {
        SeedAlignmentsOfRead * tmp = constructSeedAlignmentsOfRead();
        seedPool->seedAlignmentsOfReadArrays[threadId][readId] = tmp;
    }
    seedPool->seedAlignmentsOfReadArrays[threadId][readId]->occList->append(seedAlignmentsOfRead->occList);
    seedPool->seedAlignmentsOfReadArrays[threadId][readId]->numberOfAlignments += seedAlignmentsOfRead->numberOfAlignments;
    seedPool->totalNumberOfAlignmentsPerThread[threadId] += seedAlignmentsOfRead->numberOfAlignments;
}

void addSingleSeedPosToSeedPool(SeedPool *seedPool, DeepDP_Space::SeedPos &sp, int readOrMate, unsigned int tid) {
    OccRecord occ;
    int readId = (sp.strand_readID & 0x7FFFFFFF) + readOrMate;
    occ.strand = (sp.strand_readID >> 31) + 1;
    occ.pos = sp.pos;
    occ.seedLength = sp.paired_seedLength & 0x7FFFFFFF;

    if (seedPool->seedAlignmentsOfReadArrays[tid][readId] == NULL) {
        SeedAlignmentsOfRead * tmp = constructSeedAlignmentsOfRead();
        seedPool->seedAlignmentsOfReadArrays[tid][readId] = tmp;
    }

    seedPool->seedAlignmentsOfReadArrays[tid][readId]->occList->add(&occ);
    seedPool->seedAlignmentsOfReadArrays[tid][readId]->numberOfAlignments++;
    seedPool->totalNumberOfAlignmentsPerThread[tid]++;
}

// free SeedAlignmentOfRead without freeing up the data
void freeSeedAlignmentsOfRead(SeedAlignmentsOfRead * seedAlignmentsOfRead){
    delete seedAlignmentsOfRead->occList;
    free ( seedAlignmentsOfRead );
}

// free up SeedAlignmentOfRead memory with corresponding readId
// the salist and occlist data will be free up
int freeSeedAlignmentsOfReadWithReadId( SeedPool * seedPool, unsigned int readId )
{
    int cnt=0;
    for (int threadId=0; threadId<seedPool->numberOfThread; ++threadId)
    {   
        if (seedPool->seedAlignmentsOfReadArrays[threadId][readId] != NULL){
            SeedAlignmentsOfRead * seedAlignmentsOfRead = seedPool->seedAlignmentsOfReadArrays[threadId][readId];
            // free occList
            delete seedAlignmentsOfRead->occList;
            seedPool->totalNumberOfAlignmentsPerThread[threadId] -= seedAlignmentsOfRead->numberOfAlignments;
            free ( seedPool->seedAlignmentsOfReadArrays[threadId][readId] );
            seedPool->seedAlignmentsOfReadArrays[threadId][readId] = NULL;
            cnt = 1;
        }
    }
    return cnt;
}

void resetSeedAlignmentsOfReadWithReadId( SeedPool * seedPool, unsigned int readId )
{
    for (int threadId=0; threadId<seedPool->numberOfThread; ++threadId)
    {   
        if (seedPool->seedAlignmentsOfReadArrays[threadId][readId] != NULL) {
            SeedAlignmentsOfRead * seedAlignmentsOfRead = seedPool->seedAlignmentsOfReadArrays[threadId][readId];
            seedAlignmentsOfRead->occList->clear();
            seedPool->totalNumberOfAlignmentsPerThread[threadId] -= seedAlignmentsOfRead->numberOfAlignments;
            seedAlignmentsOfRead->numberOfAlignments = 0;
        }
    }
}

// free up SeedPool
int freeSeedPool( SeedPool * seedPool ){
    int cnt=0;
    for ( unsigned int i=0;i<seedPool->numberOfReads;++i )
        { cnt += freeSeedAlignmentsOfReadWithReadId( seedPool, i ); }
    for (int threadId=0; threadId<seedPool->numberOfThread; ++threadId)
        { free ( seedPool->seedAlignmentsOfReadArrays[threadId] ); }
    free ( seedPool->totalNumberOfAlignmentsPerThread );
    free ( seedPool->seedAlignmentsOfReadArrays );
    free ( seedPool->flag );
    free( seedPool );
    return cnt;
}

void resetSeedPool(SeedPool *seedPool) {
    for ( unsigned int i=0;i<seedPool->numberOfReads;++i )
        { resetSeedAlignmentsOfReadWithReadId( seedPool, i ); }
    memset( seedPool->flag, 0, seedPool->numberOfReads * sizeof( uint ) );
}

void filterOutUnpairedSingleReads ( UnalignedSingles * unalignedSingleReads, BothUnalignedPairs * bothUnalignedPairs, vector<int> * unalignedReadsFromSingleEndDP )
{
    int new_size = bothUnalignedPairs->totalNum * 2 + unalignedReadsFromSingleEndDP->size();
    unalignedSingleReads->readIDs = (uint *) malloc ( new_size * sizeof (uint) );
    unalignedSingleReads->size = new_size;
    unalignedSingleReads->totalNum = 0;
    int cnt=0;
    std::vector<int>::iterator iter_unalignedFromSingleEndDP = unalignedReadsFromSingleEndDP->begin();
    int iter_bothUnalignedPairs=0;
    if ( unalignedReadsFromSingleEndDP->size() >= 2 )
    {
        while ( iter_unalignedFromSingleEndDP < unalignedReadsFromSingleEndDP->end() - 1 && iter_bothUnalignedPairs < bothUnalignedPairs->totalNum )
        {
            // Sorting out unaligned reads from Single-end DP
            // 1. skip reads with odd readID
            // 2. skip read pairs with at least one pass Single-end DP
            while ( iter_unalignedFromSingleEndDP < ( unalignedReadsFromSingleEndDP->end() - 1 ) && 
                    ( ( (*iter_unalignedFromSingleEndDP) & 1 ) || (*iter_unalignedFromSingleEndDP) + 1 != (*(iter_unalignedFromSingleEndDP+1)) ) )
                { ++iter_unalignedFromSingleEndDP; }
            if ( iter_unalignedFromSingleEndDP >= unalignedReadsFromSingleEndDP->end() - 1 )
                { break; }
            if ( (*iter_unalignedFromSingleEndDP) < bothUnalignedPairs->readIDs[iter_bothUnalignedPairs] )
            {
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = (*iter_unalignedFromSingleEndDP);
                ++iter_unalignedFromSingleEndDP;
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = (*iter_unalignedFromSingleEndDP);
                ++iter_unalignedFromSingleEndDP;
            }
            else
            {
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = bothUnalignedPairs->readIDs[iter_bothUnalignedPairs];
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = bothUnalignedPairs->readIDs[iter_bothUnalignedPairs]+1;
                ++iter_bothUnalignedPairs;
            }
        }
        for ( ; iter_unalignedFromSingleEndDP < unalignedReadsFromSingleEndDP->end() - 1; ++iter_unalignedFromSingleEndDP )
        {
            if ( !( (*iter_unalignedFromSingleEndDP) & 1 ) && (*iter_unalignedFromSingleEndDP) + 1 == (*(iter_unalignedFromSingleEndDP+1)) )
            {
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = (*iter_unalignedFromSingleEndDP);
                ++iter_unalignedFromSingleEndDP;
                unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = (*iter_unalignedFromSingleEndDP);
            }
            else ++cnt;
        }
    }
    for ( ; iter_bothUnalignedPairs < bothUnalignedPairs->totalNum; ++iter_bothUnalignedPairs )
    {
        unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = bothUnalignedPairs->readIDs[iter_bothUnalignedPairs];
        unalignedSingleReads->readIDs[unalignedSingleReads->totalNum++] = bothUnalignedPairs->readIDs[iter_bothUnalignedPairs]+1;
    }
}


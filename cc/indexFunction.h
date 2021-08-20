#ifndef _INDEX_FUNCTION_H_
#define _INDEX_FUNCTION_H_

#define GRID_SAMPLING_FACTOR_2_POWER 18
#define LOOP_BOUND 1000

#include "readIndex.h"

inline char getBaseBit ( unsigned int * packedDNA, unsigned long long position )
{
    return ( ( packedDNA[position >> 4] >> ( 30 - ( ( position & 15 ) << 1 ) ) ) & 3 );
}

void getChrAndPos ( unsigned int * ambiguityMap, Translate * translate,
                    unsigned long long ambPos,
                    unsigned long long * tp, unsigned * chr_id );

unsigned long long getAmbPos ( unsigned chr_id, unsigned long long offset,
                         unsigned int * ambiguityMap, Translate * translate,
                         unsigned long long dnaLength );

int getUpperTranslateIndex ( unsigned chr_id, unsigned long long offset,
                             unsigned int * ambiguityMap, Translate * translate, int segmentCount,
                             unsigned long long dnaLength );

int getLowerTranslateIndex ( unsigned chr_id, unsigned long long offset,
                             unsigned int * ambiguityMap, Translate * translate, int segmentCount,
                             unsigned long long dnaLength );

char * getChrName ( Annotation * annotation, unsigned int numOfSeq, unsigned int chrID );

unsigned int getChrIDFromName ( Annotation * annotation, unsigned int numOfSeq, const char * chrName );

#endif

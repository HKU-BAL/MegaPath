#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "indexFunction.h"

void getChrAndPos ( unsigned int * ambiguityMap, Translate * translate,
                    unsigned long long ambPos,
                    unsigned long long * tp, unsigned * chr_id )
{
    unsigned long long correctPosition;
    unsigned long long approxIndex, approxValue;
    correctPosition = ambPos;
    approxIndex = ambPos >> GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = ambiguityMap[approxIndex];

    while ( translate[approxValue].startPos > ambPos )
    {
        approxValue--;
    }

    correctPosition -= translate[approxValue].correction;
    * tp = correctPosition;
    * chr_id = translate[approxValue].chrID;
}

unsigned long long getAmbPos ( unsigned chr_id, unsigned long long offset,
                         unsigned int * ambiguityMap, Translate * translate,
                         unsigned long long dnaLength )
{
    long long low = 0;
    long long high = dnaLength - 1;
    int loop_count = 0;
    while ( low <= high )
    {
        loop_count++;
        if ( loop_count == LOOP_BOUND )
        {
            break;
        }
        unsigned long long ambPos = ( low + high ) / 2;
        unsigned long long approxIndex, approxValue;
        unsigned long long correctPosition = ambPos;
        approxIndex = ambPos >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = ambiguityMap[approxIndex];
        while ( translate[approxValue].startPos > ambPos )
        {
            approxValue--;
        }
        correctPosition -= translate[approxValue].correction;
        unsigned chr = translate[approxValue].chrID;

        if ( chr < chr_id || ( chr == chr_id && correctPosition < offset ) )
        {
            low = ambPos + 1;
        }
        else if ( chr > chr_id || ( chr == chr_id && correctPosition > offset ) )
        {
            high = ambPos - 1;
        }
        else
        {
            return ambPos;
        }
    }

    return dnaLength + 2;
}

int getUpperTranslateIndex( unsigned chr_id, unsigned long long offset,
                            unsigned int * ambiguityMap, Translate * translate, int segmentCount,
                            unsigned long long dnaLength )
{
    int low = 0;
    int high = segmentCount - 1;
    int mid;
    while (low < high) {
        mid = (low + high) / 2;
        unsigned chr = translate[mid].chrID;
        unsigned long long correction = translate[mid].correction;
        unsigned long long correctedLastBase = translate[mid+1].startPos-1 - correction;
        if ( chr > chr_id || ( chr == chr_id && correctedLastBase >= offset ) ) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    return low;
}

int getLowerTranslateIndex( unsigned  chr_id, unsigned long long offset,
                            unsigned int * ambiguityMap, Translate * translate, int segmentCount,
                            unsigned long long dnaLength )
{
    int low = 0;
    int high = segmentCount - 1;
    int mid;
    while (low < high) {
        mid = (low + high) / 2;
        unsigned chr = translate[mid].chrID;
        unsigned long long correction = translate[mid].correction;
        unsigned long long correctedFirstBase = translate[mid].startPos - correction;
        // printf("corrected: %d %d %u\n", mid, chr, correctedFirstBase);
        if ( chr > chr_id || ( chr == chr_id && correctedFirstBase > offset ) ) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    return low - 1;
}

char * getChrName ( Annotation * annotation, unsigned int numOfSeq, unsigned int chrID )
{
    if ( chrID > numOfSeq )
    {
        return NULL;
    }
    return annotation[chrID - 1].text;
}

unsigned int getChrIDFromName ( Annotation * annotation, unsigned int numOfSeq, const char * chrName )
{
    unsigned int i;
    for ( i = 0; i < numOfSeq; ++i )
    {
        if ( strcmp ( annotation[i].text, chrName ) == 0 )
        {
            return i + 1;
        }
    }
    if ( strlen(chrName) >= 3 && chrName[0] == 'c' && chrName[1] == 'h' && chrName[2] == 'r' )
    {
        for ( i = 0; i < numOfSeq; ++i )
        {
            if ( strcmp ( annotation[i].text, chrName+3 ) == 0 )
            {
                return i + 1;
            }
        }
    }
    else
    {
        char tmp[1000];
        tmp[0] = 0;
        strcat ( tmp, "chr" );
        strcat ( tmp, chrName );
        for ( i = 0; i < numOfSeq; ++i )
        {
            if ( strcmp ( annotation[i].text, chrName ) == 0 )
            {
                return i + 1;
            }
        }
    }
    return i + 1;
}

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "readIndex.h"
#include "indexFunction.h"

void convertWordToPacked ( unsigned int & packedByte )
{
    unsigned char tempChar[4];
    * ( unsigned int * ) tempChar = packedByte;
    packedByte = ( tempChar[0] << 24 ) | tempChar[1] << 16 | tempChar[2] << 8 | tempChar[3];
}

int main ( int argc, char ** argv )
{
    fprintf ( stderr, "****************************************************************************************************\n");
    fprintf ( stderr, "Only supports 64-bit index. You can use soap4/2bwt-lib/2bwt-builder to build 64-bit index\n");
    fprintf ( stderr, "****************************************************************************************************\n");
    if ( argc != 5 )
    {
        fprintf ( stderr, "Usage:\n\t%s [Index] [chr] [start pos] [end pos]\n", argv[0] );
        return 1;
    }

    char * translateFileName;
    char * annotationFileName;
    char * pacFileName;

    unsigned int indexLength = strlen ( argv[1] );
    translateFileName = ( char * ) malloc ( indexLength + 5 );
    annotationFileName = ( char * ) malloc ( indexLength + 5 );
    pacFileName = ( char * ) malloc ( indexLength + 5 );
    sprintf ( translateFileName, "%s.tra", argv[1] );
    sprintf ( annotationFileName, "%s.ann", argv[1] );
    sprintf ( pacFileName, "%s.pac", argv[1] );

    //printf ( "Loading Index\n" );
    unsigned long long dnaLength = 0;
    unsigned int * packedDNA;
    unsigned int * ambiguityMap;
    Translate * translate;
    Annotation * annotation;
    unsigned int numOfSeq = 0;

    loadTranslate ( translateFileName, dnaLength, &ambiguityMap, &translate );
    //printf ( "Loaded Translate : %u\n", dnaLength );
    loadSeqInfo ( annotationFileName, dnaLength, &annotation, NULL, numOfSeq );
    //printf ( "Loaded Annotation : %u, #seqeunce : %u\n", dnaLength, numOfSeq );
    unsigned long long startAmb, endAmb;
    if ( strcmp ( argv[2], "-1" ) ) {
        unsigned chrId = getChrIDFromName ( annotation, numOfSeq, argv[2] );
        startAmb = getAmbPos ( chrId, strtoull( argv[3], NULL, 0 ), ambiguityMap, translate, dnaLength );
        endAmb = getAmbPos ( chrId, strtoull( argv[4], NULL, 0 ), ambiguityMap, translate, dnaLength );
    }
    else
    {
        startAmb = strtoull ( argv[3], NULL, 0 );
        endAmb = strtoull ( argv[4], NULL, 0 );
    }
    if ( endAmb - startAmb > 1000 ) return 0;

    printf ( "%llu %llu\n", startAmb, endAmb );
    FILE * pacFile = fopen64 ( pacFileName, "rb" );
    unsigned long long nByte = ( startAmb >> 4 );
    unsigned int packedByte;
    fseek ( pacFile, nByte * 4, SEEK_SET );
    fread ( &packedByte, 1, 4, pacFile );
    convertWordToPacked ( packedByte );
    #define MC_OldDnaUnpack(X,i) ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
    const char dnaMap[] = {'A','C','G','T'};
    
    for ( ;startAmb<endAmb;++startAmb )
    {
        if ( nByte != ( startAmb >> 4 ) )
        {
            fread ( &packedByte, 1, 4, pacFile );   
            convertWordToPacked ( packedByte );
            nByte = ( startAmb >> 4 );
        }
        // ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
        printf ( "%c", dnaMap[(packedByte >> ((15-((startAmb)&0xF))<<1)) & 3] );
        
    }
    printf ( "\n" );
    fclose ( pacFile );

    freeTranslate ( ambiguityMap, translate );
    freeSeqInfo ( annotation, NULL );
    free ( packedDNA );

    free ( translateFileName );
    free ( annotationFileName );
}

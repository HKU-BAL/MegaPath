#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "readIndex.h"

void * loadPackedDNA ( const char * inputFileName, unsigned long long * textLength,
                       const unsigned int convertToWordPacked )
{
    FILE * inputFile;
    unsigned char tempChar[4];
    unsigned int * packedText;
    unsigned long long packedFileLen;
    unsigned char lastByteLength;
    unsigned long long wordToProcess;
    unsigned long long i;
    unsigned int trailerBufferIn128;
    trailerBufferIn128 = 4;
    inputFile = ( FILE * ) fopen64 ( inputFileName, "rb" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "Cannot Open File : %s\n", inputFileName );
        exit ( 1 );
    }

    fseek ( inputFile, -1, SEEK_END );
    packedFileLen = ftell ( inputFile );

    if ( ( long long ) packedFileLen < 0 )
    {
        fprintf ( stderr, "Cannot determine file length : %s\n", inputFileName );
        exit ( 1 );
    }

    fread ( &lastByteLength, sizeof ( unsigned char ), 1, inputFile );
    *textLength = ( packedFileLen - 1 ) * 4 + lastByteLength;
    wordToProcess = ( ( *textLength + 63 ) >> 6 ) << 2 + trailerBufferIn128 * 4;
    packedText = ( unsigned int * ) malloc ( wordToProcess * sizeof( unsigned int ) );

    for ( i = ( *textLength ) / 16; i < wordToProcess; i++ )
    {
        packedText[i] = 0;
    }

    fseek ( inputFile, 0, SEEK_SET );
    fread ( packedText, 1, packedFileLen, inputFile );
    fclose ( inputFile );

    if ( convertToWordPacked )
    {
        for ( i = 0; i < wordToProcess; i++ )
        {
            * ( unsigned int * ) tempChar = packedText[i];
            packedText[i] = ( tempChar[0] << 24 ) | tempChar[1] << 16 | tempChar[2] << 8 | tempChar[3];
        }
    }

    return ( void * ) packedText;
}

void freePackedDNA ( void * packedDNA )
{
    free ( packedDNA );
}

size_t loadTranslateWithTranslateSize ( const char * inputFileName, unsigned long long & dnaLength,
                                     unsigned int ** ambiguityMap, Translate ** translate )
{
    FILE * inputFile;
    inputFile = ( FILE * ) fopen64 ( inputFileName, "r" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "Cannot open file %s\n", inputFileName );
        exit ( 1 );
    }

    unsigned int i;
    unsigned int gridEntries;
    unsigned int removedSegmentCount;
    fscanf ( inputFile, "%llu %d %u %u\n", &dnaLength, &i, &removedSegmentCount, &gridEntries );

    unsigned int * ambMap = ( unsigned int * ) malloc ( ( gridEntries ) * sizeof ( unsigned int ) );
    Translate * tran = ( Translate * ) malloc ( ( i + removedSegmentCount + 1 ) * sizeof ( Translate ) );

    unsigned int j = 0;

    while ( !feof ( inputFile ) && j < gridEntries )
    {
        fscanf ( inputFile, "%u\n", & ( ambMap[j] ) );
        j++;
    }

    if ( j < gridEntries )
    {
        fprintf ( stderr, "Translate missing entries!\n" );
        exit ( 1 );
    }

    j = 0;

    while ( !feof ( inputFile ) && j < i + removedSegmentCount )
    {
        fscanf ( inputFile, "%llu %u %llu\n", & ( tran[j].startPos ), & ( tran[j].chrID ), & ( tran[j].correction ) );
        j++;
    }

    if ( j < i + removedSegmentCount )
    {
        fprintf ( stderr, "Translate missing entries\n" );
        exit ( 1 );
    }
    // append senitel value for convenience
    j = i + removedSegmentCount;
    tran[j].chrID = UINT_MAX;
    tran[j].startPos = ULLONG_MAX;
    tran[j].correction = ULLONG_MAX;

    fclose ( inputFile );

    *ambiguityMap = ambMap;
    *translate = tran;

    return i + removedSegmentCount + 1; // +1 for the sentinel value
}

void loadTranslate( const char * inputFileName, unsigned long long & dnaLength,
                    unsigned int ** ambiguityMap, Translate ** translate )
{
    loadTranslateWithTranslateSize( inputFileName, dnaLength, ambiguityMap, translate );
}

void freeTranslate ( unsigned int * ambiguityMap, Translate * translate )
{
    free ( ambiguityMap );
    free ( translate );
}

void loadSeqInfo ( const char * inputFileName, unsigned long long & dnaLength,
                   Annotation ** annotation, SeqOffset ** seqOffset, unsigned int & numOfSeq )
{
    FILE * inputFile;
    inputFile = ( FILE * ) fopen64 ( inputFileName, "r" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "Cannot open file %s\n", inputFileName );
        exit ( 1 );
    }

    unsigned int randomSeed;

    fscanf ( inputFile, "%llu %d %u\n", &dnaLength, &numOfSeq, &randomSeed );

    if ( numOfSeq == 0 )
    {
        fprintf ( stderr, "Annotation empty entry\n" );
        exit ( 1 );
    }

    Annotation * ann = ( Annotation * ) malloc ( ( numOfSeq + 1 ) * sizeof ( Annotation ) );
    memset ( ann, 0, ( numOfSeq + 1 ) * sizeof ( Annotation ) );
    SeqOffset * seq = ( SeqOffset * ) malloc ( ( numOfSeq + 1 ) * sizeof ( SeqOffset ) );
    memset ( seq, 0, ( numOfSeq + 1 ) * sizeof ( SeqOffset ) );

    int i = 0;
    int j, k;
    while ( !feof ( inputFile ) && i < numOfSeq )
    {
        fscanf ( inputFile, "%u ", &( ann[i].gi ) );
        fgets ( ann[i].text, MAX_SEQ_NAME_LENGTH, inputFile );
        fscanf ( inputFile, "%llu %llu %d\n", &( seq[i].startPos ), &( seq[i].endPos ), &( seq[i + 1].firstAmbiguityIndex ) );
        seq[i].lastAmbiguityIndex = seq[i + 1].firstAmbiguityIndex;
        seq[i].endPos = seq[i].startPos + seq[i].endPos - 1;

        j = 0;
        while ( j < MAX_SEQ_NAME_LENGTH )
        {
            if ( ann[i].text[j] == '\n'
                 || ann[i].text[j] == '\t'
                 || ann[i].text[j] == ' '
                 || ann[i].text[j] == '\0' )
            {
                break;
            }
            j++;
        }
        ann[i].text[j] = '\0';

        j++;
        if ( ann[i].text[j] != '\0' )
        {
            k = 0;
            while ( j < MAX_SEQ_NAME_LENGTH )
            {
                ann[i].decoratedText[k] = ann[i].text[j];
                k++;
                j++;
            }
        }

        i++;
    }

    if ( i < numOfSeq )
    {
        fprintf ( stderr, "Annotation missing entries\n" );
        exit ( 1 );
    }
    ann[i].gi = 0;
    seq[0].firstAmbiguityIndex = 1;

    if ( annotation )
    {
        *annotation = ann;
    }
    else
    {
        free ( ann );
    }
    if ( seqOffset )
    {
        *seqOffset = seq;
    }
    else
    {
        free ( seq );
    }
}

void freeSeqInfo ( Annotation * annotation, SeqOffset * seqOffset )
{
    free ( annotation );
    free ( seqOffset );
}

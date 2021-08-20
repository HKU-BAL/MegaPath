#ifndef _READ_INDEX_H_
#define _READ_INDEX_H_

#define MAX_SEQ_NAME_LENGTH 256

typedef struct Translate
{
    unsigned long long startPos;
    unsigned int chrID;
    unsigned long long correction;
} Translate;

typedef struct Annotation
{
    int gi;
    char text[MAX_SEQ_NAME_LENGTH + 1];
    char decoratedText[MAX_SEQ_NAME_LENGTH + 1];
} Annotation;

typedef struct SeqOffset
{
    unsigned long long startPos;
    unsigned long long endPos;
    int firstAmbiguityIndex;    // The index for the first ambiguity that starts on or after the sequence
    int lastAmbiguityIndex;     // The index for the last ambiguity that ends on or before the sequence
} SeqOffset; 

void * loadPackedDNA ( const char * inputFileName, unsigned long long * textLength,
                       const unsigned int convertToWordPacked );

void freePackedDNA ( void * packedDNA );

void loadTranslate ( const char * inputFileName, unsigned long long &dnaLength,
                     unsigned int ** ambiguityMap, Translate ** translate );

void freeTranslate ( unsigned int * ambiguityMap, Translate * translate );

void loadSeqInfo ( const char * inputFileName, unsigned long long & dnaLength,
                   Annotation ** annotation, SeqOffset ** seqOffset, unsigned int & numOfSeq );

void freeSeqInfo ( Annotation * annotation, SeqOffset * seqOffset );

#endif

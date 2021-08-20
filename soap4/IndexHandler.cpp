/*
 *
 *    IndexHandler.c
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

#include "IndexHandler.h"

void INDEXFillCharMap(unsigned char charMap[255]) {

    int i;

    for (i=0; i<255; i++) {
        charMap[i] = 0;
    }

    //Replacing DNA_CHAR_SIZE with ALPHABET_SIZE as SOAP3 has been developed as
    // a restrictive tool..
    for (i=0; i<ALPHABET_SIZE; i++) {
        charMap[dnaChar[i]] = (unsigned char)i;
        charMap[dnaChar[i] - 'A' + 'a'] = (unsigned char)i;
    }
    
    charMap['U'] = charMap['T'];
    charMap['N'] = charMap['G']; // N -> G
    
    charMap['u'] = charMap['t'];
    charMap['n'] = charMap['g']; // N -> G
}

// To load the index
Soap3Index * INDEXLoad ( IniParams * ini_params, char * indexName, char isShareIndex )
{
    MMMasterInitialize ( 3, 0, FALSE, NULL );
    fprintf (stderr, "[Main] loading index into host...\n" );
    double startTime = setStartTime();
    Soap3Index * index = ( Soap3Index * ) malloc ( sizeof ( Soap3Index ) );
    SRAIndex * sraIndex = ( SRAIndex * ) malloc ( sizeof ( SRAIndex ) );

    index->sraIndex = sraIndex;
    index->mmPool = MMPoolCreate ( INDEX_MMPOOL_SIZE );
    IndexFileNames indexFilenames;
    INDEXProcessFilenames ( &indexFilenames, indexName, ini_params );
    INDEXLoad ( index->mmPool, indexFilenames, & ( sraIndex->bwt ), & ( sraIndex->rev_bwt ), & ( sraIndex->hsp ), & ( sraIndex->lookupTable ), & ( sraIndex->rev_lookupTable ),  isShareIndex );

    
    /////////////////////////////////
    // Allocate and Fill CharMap
    index->charMap = (unsigned char*) malloc(sizeof(unsigned char) * 256);
    INDEXFillCharMap(index->charMap);
    
    /////////////////////////////////
    // Allocate HSPAuxilliary
    index->sraIndex->hspaux = ( HSPAux* ) malloc( sizeof ( HSPAux ) );
    
    double loadIndexTime = getElapsedTime ( startTime );
    fprintf (stderr, "[Main] Finished loading index into host.\n" );
    fprintf (stderr, "[Main] Loading time : %9.4f seconds\n\n", loadIndexTime );
    return index;
}


// To free the index
void INDEXFree ( Soap3Index * index, char isShareIndex )
{
    SRAIndex * sraIndex = index->sraIndex;
    uint numOfOccValue = ( sraIndex->bwt->textLength + GPU_OCC_INTERVAL - 1 ) / GPU_OCC_INTERVAL + 1;
    BWTFree ( index->mmPool, sraIndex->bwt, 0);
    // BWTFree ( index->mmPool, sraIndex->rev_bwt,0 );
    HSPFree ( index->mmPool, sraIndex->hsp, 1);

    {
        LTFree ( sraIndex->lookupTable );
        // LTFree ( sraIndex->rev_lookupTable );
    }

    MMPoolFree ( index->mmPool );
    free ( index->sraIndex->hspaux );
    free ( index->charMap );
    free ( sraIndex );
    free ( index );
}

void INDEXLoad ( MMPool * mmPool, IndexFileNames indexFilenames, BWT ** bwt, BWT ** revBwt, HSP ** hsp, LT ** lkt, LT ** revLkt, 
                 char isShareIndex )
{
    { ( *bwt ) = BWTLoad ( mmPool, 0, indexFilenames.bwtCodeFileName, indexFilenames.occValueFileName, indexFilenames.saCodeFileName ); }

    free ( indexFilenames.bwtCodeFileName );
    free ( indexFilenames.occValueFileName );
    free ( indexFilenames.saCodeFileName );
    free ( indexFilenames.mmapOccValueFileName );

    // { ( *revBwt ) = BWTLoad ( mmPool, 0, indexFilenames.revBwtCodeFileName, indexFilenames.revOccValueFileName, NULL ); }

    free ( indexFilenames.revBwtCodeFileName );
    free ( indexFilenames.revOccValueFileName );
    free ( indexFilenames.mmapRevOccValueFileName );

    { ( *hsp ) = HSPLoad ( mmPool, indexFilenames.packedDnaFileName, indexFilenames.annotationFileName, indexFilenames.ambiguityFileName, indexFilenames.translateFileName, 1); }
    
    free ( indexFilenames.packedDnaFileName );
    free ( indexFilenames.mmapPackedDnaFileName );
    free ( indexFilenames.annotationFileName );
    free ( indexFilenames.ambiguityFileName );
    free ( indexFilenames.translateFileName );
    unsigned int lookupWordSize = 1 << ( LOOKUP_SIZE * 2 );
    FILE *lookupTableFile, *revLookupTableFile;

    {
        ( *lkt ) = LTLoad ( indexFilenames.lookupTableFileName );
    }

    if ( ( *lkt )->tableSize != LOOKUP_SIZE )
    {
        fprintf ( stderr, "SOAP3 only supports lookup table size of %u, please re-build index.\n", LOOKUP_SIZE );
        exit ( 1 );
    }

    free ( indexFilenames.lookupTableFileName );
    /////////////////// READ REVERSE ////////////////////////////////////

    // {
    //     ( *revLkt ) = LTLoad ( indexFilenames.revLookupTableFileName );
    // }

    // if ( ( *revLkt )->tableSize != LOOKUP_SIZE )
    // {
    //     fprintf ( stderr, "SOAP3 only supports lookup table size of %u, please re-build index.\n", LOOKUP_SIZE );
    //     exit ( 1 );
    // }
    
    free ( indexFilenames.revLookupTableFileName );
    // SHOW THE MEMORY USAGE ON HOST
#ifdef BGS_OUTPUT_MEMORY_USAGE
    // for BWT
    fprintf (stderr, "BWT code size    : %lu bytes (%9.4f M)\n", ( *bwt )->bwtSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->bwtSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    fprintf (stderr, "Occ value size   : %lu bytes (%9.4f M)\n", ( ( *bwt )->occSizeInWord + ( *bwt )->occMajorSizeInWord ) * sizeof ( unsigned int ), ( float ) ( ( *bwt )->occSizeInWord + ( *bwt )->occMajorSizeInWord ) * sizeof ( unsigned int ) / 1024.0 / 1024.0 );

    if ( ( *bwt )->saValueSizeInWord > 0 )
    {
        fprintf (stderr, "SA value size    : %lu bytes (%9.4f M)\n", ( *bwt )->saValueSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->saValueSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    // for REVBWT
    fprintf (stderr, "Rev BWT code size    : %lu bytes (%9.4f M)\n", ( *revBwt )->bwtSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->bwtSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    fprintf (stderr, "Rev Occ value size   : %lu bytes (%9.4f M)\n", ( ( *revBwt )->occSizeInWord + ( *revBwt )->occMajorSizeInWord ) * sizeof ( unsigned int ), ( float ) ( ( *revBwt )->occSizeInWord + ( *revBwt )->occMajorSizeInWord ) * sizeof ( unsigned int ) / 1024.0 / 1024.0 );

    if ( ( *revBwt )->saValueSizeInWord > 0 )
    {
        fprintf (stderr, "Rev SA value size    : %lu bytes (%9.4f M)\n", ( *revBwt )->saValueSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->saValueSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    // for HSP
    uint trailerBufferInWord = 1;
    uint trailerBufferIn128 = ( trailerBufferInWord + 3 ) / 4 * 4;
    uint textLength = ( *hsp )->dnaLength;
    // hsp->packedDNA
    uint hsp_mem_usage = ( ( textLength + 64 - 1 ) / 64 * 4 + trailerBufferIn128 * 4 ) * sizeof ( unsigned int );
    // hsp->seqOffset
    hsp_mem_usage += ( ( *hsp )->numOfSeq + 1 ) * sizeof ( SeqOffset );
    // hsp->annotation
    hsp_mem_usage += ( ( *hsp )->numOfSeq + 1 ) * sizeof ( Annotation );
    // hsp->ambiguity
    hsp_mem_usage += ( ( *hsp )->numOfAmbiguity + 2 ) * sizeof ( Ambiguity );
    fprintf (stderr, "[Memory usage in Host] hsp: %u bytes (%9.4f M)\n", hsp_mem_usage, ( float ) hsp_mem_usage / 1024.0 / 1024.0 );
    fprintf (stderr, "[Memory usage in Host] occValue: %i bytes (%9.4f M)\n", numOfOccValue * ALPHABET_SIZE * sizeof ( unsigned int ), ( float ) numOfOccValue * ALPHABET_SIZE * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    fprintf (stderr, "[Memory usage in Host] lkt: %i bytes (%9.4f M)\n", sizeof ( unsigned int ) * lookupWordSize, ( float ) sizeof ( unsigned int ) * lookupWordSize / 1024.0 / 1024.0 );
    fprintf (stderr, "[Memory usage in Host] revOccValue: %i bytes (%9.4f M)\n", numOfOccValue * ALPHABET_SIZE * sizeof ( uint ), ( float ) numOfOccValue * ALPHABET_SIZE * sizeof ( uint ) / 1024.0 / 1024.0 );
    fprintf (stderr, "[Memory usage in Host] revLkt: %i bytes (%9.4f M)\n", sizeof ( uint ) * lookupWordSize, ( float ) sizeof ( uint ) * lookupWordSize / 1024.0 / 1024.0 );
    fprintf (stderr, "[Memory usage in Host] occ for all cpu threads: %i bytes (%9.4f M)\n", sizeof ( OCC ) *ini_params.Ini_NumOfCpuThreads, ( float ) sizeof ( OCC ) *ini_params.Ini_NumOfCpuThreads / 1024.0 / 1024.0 );
#endif
}



void INDEXProcessFilenames ( IndexFileNames * indexFilenames, char * indexName, IniParams * ini_params )
{
    char saExtension[10] = ".sa";

    if ( ini_params != NULL )
    {
        strcpy ( saExtension, ini_params->Ini_SaValueFileExt );
    }

    int indexNameLen = strlen ( indexName );
    // indexFilenames->iniFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->bwtCodeFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->occValueFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->lookupTableFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->revBwtCodeFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->revOccValueFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->revLookupTableFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->saCodeFileName = ( char * ) malloc ( indexNameLen + 5 );
    // indexFilenames->memControlFileName = ( char * ) malloc ( indexNameLen + 5 );
    //For HSP
    indexFilenames->packedDnaFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->annotationFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->ambiguityFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->translateFileName = ( char * ) malloc ( indexNameLen + 5 );
    // For mmap
    indexFilenames->mmapOccValueFileName = ( char * ) malloc ( indexNameLen + 25 );
    indexFilenames->mmapRevOccValueFileName = ( char * ) malloc ( indexNameLen + 35 );
    indexFilenames->mmapPackedDnaFileName = ( char * ) malloc ( indexNameLen + 25 );
    strcpy ( indexFilenames->bwtCodeFileName, indexName );
    strcpy ( indexFilenames->bwtCodeFileName + indexNameLen, ".bwt" );
    strcpy ( indexFilenames->occValueFileName, indexName );
    strcpy ( indexFilenames->occValueFileName + indexNameLen, ".fmv" );
    strcpy ( indexFilenames->lookupTableFileName, indexName );
    strcpy ( indexFilenames->lookupTableFileName + indexNameLen, ".lkt" );
    strcpy ( indexFilenames->revBwtCodeFileName, indexName );
    strcpy ( indexFilenames->revBwtCodeFileName + indexNameLen, ".rev.bwt" );
    strcpy ( indexFilenames->revOccValueFileName, indexName );
    strcpy ( indexFilenames->revOccValueFileName + indexNameLen, ".rev.fmv" );
    strcpy ( indexFilenames->revLookupTableFileName, indexName );
    strcpy ( indexFilenames->revLookupTableFileName + indexNameLen, ".rev.lkt" );
    strcpy ( indexFilenames->saCodeFileName, indexName );
    strcpy ( indexFilenames->saCodeFileName + indexNameLen, saExtension );
    strcpy ( indexFilenames->packedDnaFileName, indexName );
    strcpy ( indexFilenames->packedDnaFileName + indexNameLen, ".pac" );
    strcpy ( indexFilenames->annotationFileName, indexName );
    strcpy ( indexFilenames->annotationFileName + indexNameLen, ".ann" );
    strcpy ( indexFilenames->ambiguityFileName, indexName );
    strcpy ( indexFilenames->ambiguityFileName + indexNameLen, ".amb" );
    strcpy ( indexFilenames->translateFileName, indexName );
    strcpy ( indexFilenames->translateFileName + indexNameLen, ".tra" );
    strcpy ( indexFilenames->mmapOccValueFileName, indexName );
    strcpy ( indexFilenames->mmapOccValueFileName + indexNameLen, ".fmv.mmap" );
    strcpy ( indexFilenames->mmapRevOccValueFileName, indexName );
    strcpy ( indexFilenames->mmapRevOccValueFileName + indexNameLen, ".rev.fmv.mmap" );
    strcpy ( indexFilenames->mmapPackedDnaFileName, indexName );
    strcpy ( indexFilenames->mmapPackedDnaFileName + indexNameLen, ".pac.mmap" );
#ifdef BGS_OUTPUT_FILENAMES_TO_SCREEN
    fprintf (stderr, "[Main] Using BWT index ....                %s\n", bwtCodeFileName );
    fprintf (stderr, "[Main] Using Occ index ....                %s\n", occValueFileName );
    fprintf (stderr, "[Main] Using Lookup index ....             %s\n", lookupTableFileName );
    fprintf (stderr, "[Main] Using Rev-BWT index ....            %s\n", revBwtCodeFileName );
    fprintf (stderr, "[Main] Using Rev-Occ index ....            %s\n", revOccValueFileName );
    fprintf (stderr, "[Main] Using Rev-Lookup index ....         %s\n", revLookupTableFileName );
    fprintf (stderr, "[Main] Using Rev-Gpu-Occ index ....        %s\n", revGpuOccValueFileName );
    fprintf (stderr, "[Main] Using SA index ....                 %s\n", saCodeFileName );
    fprintf (stderr, "[Main] Using PackedDna index ....          %s\n", packedDnaFileName );
    fprintf (stderr, "[Main] Using Annotation index ....         %s\n", annotationFileName );
    fprintf (stderr, "[Main] Using Ambiguity index ....          %s\n", ambiguityFileName );
    fprintf (stderr, "[Main] Using Translation index ....        %s\n", translateFileName );
#endif
}

/*
 *
 *    SOAP3-DP.cu
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <zlib.h>
#include <algorithm>
#include <assert.h>

#include "VERSION.h"
#include "Release.h"
#include "CPUfunctions.h"
#include "alignment.h"
#include "SAM.h"
#include "samtools-0.1.18/bam.h"

#include "BGS-IO.h"
#include "2bwt-flex/SRAArguments.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"
#include "OutputDPResult.h"
#include "DV-SemiDP.h"
#include "DV-DPForBothUnalign.h"
#include "DV-DPForSingleReads.h"
#include "SeedPool.h"
#include "aio_thread.h"



int main ( int argc, char ** argv )
{
    // ======================================================================================
    // | VARIABLE DECLARATION                                                               |
    // ======================================================================================
    // local variables used in main, like BWT indexes and such.
    int i;
    double startTime, indexLoadTime, readLoadTime;
    double lastEventTime;
    double totalReadLoadTime = 0.0;
    double totalAlignmentTime = 0.0;
    double totalTrimAlignmentTime = 0.0;
    char * queryFileName = "";
    char * queryFileName2 = "";
    gzFile gzQueryFile;
    gzFile gzQueryFile2;
    bamFile bamQueryFile;
    bam_header_t * bamHeader;
    bam1_t * bam;
    char isFastq = 0;
    // user-specified maximum read length
    // and number of words per query
    uint maxReadLength;
    uint wordPerQuery;
    // uint numQueries;
    ullint roundUp;
    ullint totalQueryLength;
    uint * queries;
    uint * readLengths;
    uint * readIDs;
    char * upkdQualities;
    char ** queryNames;
    char ** queryComments;
    DPInfoForReads * dpInfoForReads;
    uint * queries0;
    uint * readLengths0;
    uint * readIDs0;
    char * upkdQualities0;
    DPInfoForReads * dpInfoForReads0;
    uint * queries1;
    uint * readLengths1;
    uint * readIDs1;
    char * upkdQualities1;
    char queryFileBuffer[INPUT_BUFFER_SIZE];
    char queryFileBuffer2[INPUT_BUFFER_SIZE];
    unsigned long long numOfAnswer;
    unsigned int numOfAlignedRead;
    unsigned int numOfUnAlignedPairs;
    SeedPool * seedPool;
    uint detected_read_length = 0;
    // for single-end reads
    // it represents the max read length for the first ten reads (i.e. 0, 1, ..., 9)
    // for paired-end reads
    // the max read length for the first ten reads with even readIDs (i.e. 0,2,...,18)
    uint detected_read_length2 = 0;
    // for paired-end reads only
    // the max read length for the first ten reads with odd readIDs (i.e. 1,3,...,19)
    // Declare variables and set up preference for device.
    // #ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    // cudaEvent_t start, stop;
    // float time1, totalDeviceTime;
    // #endif
    // Accumulated number of reads aligned and alignments
    ullint totalReadsAlignedForInputReads = 0;
    ullint totalAnsForInputReads = 0;
    // Input parameters
    InputOptions input_options;
    // The inputs for multi option
    MultiInputItem * multiInput = NULL;
    int expNum = 0;
    int currExp = 0;
    // ======================================================================================
    // | Configuration on GPU functions                                                     |
    // ======================================================================================
    // not much effect, also this is not supported by old version of cuda library
    // like the machines in TJ
    // thus these are depreciated
    // cudaDeviceSetCacheConfig ( cudaFuncCachePreferL1 ); 
    // cudaFuncSetCacheConfig(kernel, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_1, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_2, cudaFuncCachePreferShared);
    
    fprintf (stderr, "\n[Main] %s v%d.%d.%d (%s), Commit Hash: %s\n", PROJECT_NAME, PROJECT_MAJOR, PROJECT_MINOR, PROJECT_REV, PROJECT_SPECIAL, COMMIT_HASH );
    // printf("[Main] Finished parsing ini file %s.\n\n", iniFileName);
    // ======================================================================================
    // | CHECK THE INPUT ARGUMENTS                                                          |
    // ======================================================================================
    bool inputValid = parseInputArgs ( argc, argv, input_options );

    if ( !inputValid )
    { exit ( 1 ); }

    // ======================================================================================
    // | PARSING CONFIGURATION FILE                                                         |
    // ======================================================================================
    char * iniFileName = input_options.iniFileName;
    IniParams ini_params;
    if ( iniFileName == NULL ) {
        iniFileName = ( char * ) malloc ( strlen ( argv[0] ) + 5 );
        strcpy ( iniFileName, argv[0] );
        strcpy ( iniFileName + strlen ( argv[0] ), ".ini" );
    }

    if ( ParseIniFile ( iniFileName, ini_params ) != 0 )
    {
        fprintf ( stderr, "Failed to parse config file ... %s\n", iniFileName );
        return 1;
    }

    // TODO merge ini parsing and option parsing
    ini_params.Ini_NumOfCpuThreads = input_options.numCpuThreads;
    fprintf(stderr, "Number of CPU threads: %d\n", ini_params.Ini_NumOfCpuThreads);

    if ( input_options.iniFileName == NULL) {
        free ( iniFileName );
    }

    // for multi mode
    if ( input_options.isReadList == 1 )
    {
        multiInput = loadMultiInputFile ( input_options.queryFileName, ( input_options.readType == PAIR_END_READ ) ? 1 : 0,
                                          input_options.isReadBAM, expNum );
        updateInputOption ( ( &input_options ), multiInput, currExp++ );
    }

    // get the name of the query file (and the second query file for pair-ended reads)
    queryFileName = input_options.queryFileName;

    if ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 )
    {
        queryFileName2 = input_options.queryFileName2;
    }
    // printParameters( input_options , ini_params);
    maxReadLength = input_options.maxReadLength;

    if ( input_options.readType == SINGLE_READ )
    { fprintf (stderr, "[Main] Loading read file %s\n", queryFileName ); }
    else if ( input_options.isReadBAM == 0 )
    { fprintf (stderr, "[Main] Loading read files %s and %s\n", queryFileName, queryFileName2 ); }

    // restriction on the number of hits of each end for pairing
#ifdef NO_CONSTRAINT_SINGLE_READ_NUM_FOR_PAIRING

    if ( input_options.readType == PAIR_END_READ )
    { ini_params.Ini_MaxOutputPerRead = 0xFFFFFFFF; }

#endif
    
    // ======================================================================================
    // | VARIABLES SETTING AND INITIALISATION                                               |
    // ======================================================================================
    // determine number of words per query. rounded up to power of 2
    wordPerQuery = getWordPerQuery(maxReadLength);
    uint maxNumQueries = MAX_NUM_BATCH * NUM_BLOCKS * THREADS_PER_BLOCK / 6; // 2M
    // maxNumQueries has to be divible by 2
    maxNumQueries = maxNumQueries / 2 * 2;
    // For Output filenames
    char * outputFileName[MAX_NUM_CPU_THREADS];

    for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
    { outputFileName[i] = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 10 ); }

    // For Output of the DP results
    char * outputDPFileName;
    outputDPFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 9 );
    // For Output of the unpaired results
    char * outputUnpairFileName;
    outputUnpairFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 8 );

    for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
    {
        sprintf ( outputFileName[i], "%s.gout.%d", input_options.outputPrefix, i + 1 );
    }

    sprintf ( outputDPFileName, "%s.dpout.1", input_options.outputPrefix );
    sprintf ( outputUnpairFileName, "%s.unpair", input_options.outputPrefix );
    // ======================================================================================
    // | STRUCTURES FOR SEMI-GLOBAL ALIGNMENT DP                                            |
    // ======================================================================================
    // Declare the structure for storing the read IDs of the first end of the pairs
    // which both ends cannot be aligned (if the input read is paired-end)
    BothUnalignedPairsArrays * bothUnalignedPairsArrays = constructBothUnalignedPairsArrays ( ini_params.Ini_NumOfCpuThreads + 1 );
    // OR for single end cannot be aligned (if the input read is single-end)
    UnalignedSinglesArrays * unalignedSinglesArrays = bothUnalignedPairsArrays; // same array but just different name
    // Declare the structure for storing the alignment results
    // for the pairs of reads with both ends have hits but do not have valid insert size or proper strands.
    // Parameters for DP
    DPParameters dpParameters;
    getParameterForAllDP ( dpParameters, input_options, ini_params );
    // ======================================================================================
    // | INDEX LOADING                                                                      |
    // ======================================================================================
    //Start measuring runtime..
    startTime = setStartTime();
    lastEventTime = startTime;
    Soap3Index * index = INDEXLoad ( &ini_params, input_options.indexName, ini_params.Ini_shareIndex );
    HSP * hsp = index->sraIndex->hsp;
    HSPAux * hspaux = index->sraIndex->hspaux;
    hspaux->megapathMode = input_options.megapathMode; // turn on megapath mode
    hspaux->outputBAM = input_options.outputBAM; 
    hspaux->top_percentage = std::min(100, std::max(0, input_options.top)) / 100.0;
    fprintf(stderr, "[Main] top_percentage: %f\n", hspaux->top_percentage);
    indexLoadTime = getElapsedTime ( startTime );
    fprintf (stderr, "[Main] Reference sequence length : %llu\n\n", index->sraIndex->bwt->textLength );
    lastEventTime = indexLoadTime;
    roundUp = ( maxNumQueries + 31 ) / 32 * 32;
    totalQueryLength = roundUp * wordPerQuery;
    // ==============================================================
    // | QUALITY CONSTANT and DP MATCH SCORE and alignment type     |
    // ==============================================================
    int quality_constant = 0;

    if ( input_options.isIlluminaQual == 1 )
    { quality_constant = ILLUMINA_QUAL_CONST;assert(0 && "illumina quality is not supported yet"); }

    hspaux->quality_constant = quality_constant;

    hspaux->dpMatchScore = ini_params.Ini_MatchScore;
    hspaux->dpMisMatchScore = ini_params.Ini_MismatchScore;
    hspaux->alignmentType = input_options.alignmentType;
    hspaux->readType = input_options.readType;
    hspaux->peMaxOutputPerRead = ini_params.Ini_PEMaxOutputPerPair;

//    if ( input_options.readType == SINGLE_READ )
    { hspaux->ProceedDPForTooManyHits = ini_params.Ini_proceedDPForTooManyHits; }
//    else
//    { hspaux->ProceedDPForTooManyHits = 0; }

    // For Mapping Quality Score Calculation
    hspaux->x0_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->x1_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->mismatch_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->minMAPQ = ini_params.Ini_minMAPQ;
    hspaux->maxMAPQ = ini_params.Ini_maxMAPQ;
    hspaux->bwaLikeScore = ini_params.Ini_bwaLikeScore;
    hspaux->maxLenReadName = ini_params.Ini_maxReadNameLen;

    if ( hspaux->bwaLikeScore )
    { bwase_initialize ( hspaux->g_log_n ); }

    // print MD string and NM tag? 
    hspaux->isPrintMDNM = input_options.isPrintMDNM; 
    
    // For the SAM output information
    hspaux->readGroup = input_options.readGroup;

    if ( strlen ( hspaux->readGroup ) == 0 )
    { hspaux->readGroup = input_options.queryFileName; }

    hspaux->sampleName = input_options.sampleName;

    if ( strlen ( hspaux->sampleName ) == 0 )
    { hspaux->sampleName = DEFAULT_SAMPLE_NAME; }

    hspaux->readGrpOption = input_options.readGrpOption;
    // ==================================================================
    // | For construction of arrays to store the unaligned read IDs
    // | for proceeding single deep-dp on them
    // ==================================================================
    hspaux->readsIDForSingleDP = ( BothUnalignedPairsArrays * ) constructBothUnalignedPairsArrays ( 1 );
    hspaux->allHits = ( AllHits * ) constructAllHits(); // for storing the corresponding algnments
    // ======================================================================================
    // | ALLOCATE MEMORY FOR THE ARRAYS                                                     |
    // ======================================================================================
    queries0 = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
    readLengths0 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    readIDs0 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    upkdQualities0 = ( char * ) malloc ( roundUp * maxReadLength * sizeof ( char ) );
    memset ( upkdQualities0, 0, roundUp * maxReadLength );
    dpInfoForReads0 = constructDPInfoForReads ( roundUp );
    queries1 = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
    readLengths1 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    readIDs1 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    upkdQualities1 = ( char * ) malloc ( roundUp * maxReadLength * sizeof ( char ) );
    memset ( upkdQualities1, 0, roundUp * maxReadLength );

    char ** queryNames0 = (char **) calloc(roundUp, sizeof(char*));
    char ** queryNames1 = (char **) calloc(roundUp, sizeof(char*));
    char ** queryComments0 = (char **) calloc(roundUp, sizeof(char*));
    char ** queryComments1 = (char **) calloc(roundUp, sizeof(char*));

    uint maxBatchSize = NUM_BLOCKS * THREADS_PER_BLOCK * QUERIES_PER_THREAD; // queries processed in one kernel call
    // maxBatchSize has to be divible by 2
    maxBatchSize = maxBatchSize / 2 * 2;
    // #ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    // totalDeviceTime=0;
    // #endif
    // initialize the output files
    bam_header_t samOutputHeader;
    samfile_t * samOutputFilePtr[MAX_NUM_CPU_THREADS] = {NULL};
    samfile_t * samOutputDPFilePtr = NULL;
    samfile_t * samOutputUnpairFilePtr = NULL;

    if (input_options.outputBAM) {
        SAMOutputHeaderConstruct ( &samOutputHeader, hsp, hspaux, maxReadLength );
        for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ ) 
            assert(samOutputFilePtr[i] = samopen ( outputFileName[i],  "wb", &samOutputHeader )); 
        // For DP
        assert( samOutputDPFilePtr = samopen ( outputDPFileName, "wb", &samOutputHeader ));
        // For Unpaired reads
        assert(samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wb" , &samOutputHeader ));
    }
    // ======================================================================================
    // | LOADING INPUT SHORT READ FILE                                                      |
    // ======================================================================================
    size_t bufferSize;
    uint bufferIndex;
    char queryChar;
    size_t bufferSize2;
    uint bufferIndex2;
    char queryChar2;

    kseq_t *ks1, *ks2;

    if ( input_options.isReadBAM )
    {
        // the query file is in BAM format
        bamQueryFile = bam_open ( queryFileName, "r" );
        bamHeader = bam_header_init();
        bamHeader = bam_header_read ( bamQueryFile );
        bam = bam_init1();
    }
    else
    {
        //gzQueryFile = ( gzFile * ) gzopen ( queryFileName, "r" );
        gzQueryFile = gzopen ( queryFileName, "r" );

        if ( gzQueryFile == NULL ) { fprintf ( stderr, "Cannot open queryFile\n" ); exit ( 1 );}
        ks1 = kseq_init(gzQueryFile);

        if ( input_options.readType == PAIR_END_READ )
        {
            // pair-ended reads
            //gzQueryFile2 = ( gzFile * ) gzopen ( queryFileName2, "r" );
            gzQueryFile2 = gzopen ( queryFileName2, "r" );
            if ( gzQueryFile2 == NULL ) { fprintf ( stderr, "Cannot open queryFile2\n" ); exit ( 1 );}
            ks2 = kseq_init(gzQueryFile2);
        }
    }

    uint accumReadNum = 0;
    uint numQueries;
    // create buffers

    InputReadsBuffer * buffer0 = InputReadsBufferFullCreate ( maxReadLength,
                                 maxNumQueries, wordPerQuery, quality_constant, queries0,
                                 readLengths0, readIDs0, upkdQualities0, queryNames0, queryComments0,
                                 isFastq, ini_params.Ini_maxReadNameLen );
    InputReadsBuffer * buffer1 = InputReadsBufferFullCreate ( maxReadLength,
                                 maxNumQueries, wordPerQuery, quality_constant, queries1,
                                 readLengths1, readIDs1, upkdQualities1, queryNames1, queryComments1,
                                 isFastq, ini_params.Ini_maxReadNameLen );
    AIOInputBuffer * aiob = AIOInputBufferCreate ( buffer0, buffer1 );
    InputFilePointers * ifp = InputFilePointersCreate();

    if ( input_options.isReadBAM )
    {
        InputFilePointersSetBam ( ifp, bamQueryFile, bamHeader, bam );
    }
    else if ( input_options.readType == SINGLE_READ )
    {
        InputFilePointersSetSingle ( ifp, gzQueryFile );
    }
    else
    {
        InputFilePointersSetPair ( ifp, ks1, ks2 );
    }

    aiob->reads = ifp;
    InputReadsBuffer * readyReadsBuffer;
    //create io thread
    AIOInputThreadCreate ( aiob, bufferSize, bufferIndex, index->charMap, queryChar,
                           queryFileBuffer, bufferSize2, bufferIndex2, queryChar2,
                           queryFileBuffer2 );


    seedPool = constructSeedPool ( maxNumQueries, ini_params.Ini_NumOfCpuThreads );

    //load reads
    readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
    queries = readyReadsBuffer->queries;
    readLengths = readyReadsBuffer->readLengths;
    readIDs = readyReadsBuffer->readIDs;
    upkdQualities = readyReadsBuffer->upkdQualities;
    queryNames = readyReadsBuffer->queryNames;
    queryComments = readyReadsBuffer->queryComments;
    isFastq = readyReadsBuffer->isFastq;
    hspaux->isFastq = isFastq;
    numQueries = readyReadsBuffer->filledNum;

    //logically no use here
    dpInfoForReads = dpInfoForReads0;
    while ( numQueries > 0 )
    {
        fprintf (stderr, "[Main] Loaded %u short reads from the query file.\n", numQueries );
        readLoadTime = getElapsedTime ( startTime );
        fprintf (stderr, "[Main] Elapsed time on host : %9.4f seconds\n\n", readLoadTime - lastEventTime );
        totalReadLoadTime += readLoadTime - lastEventTime;
        lastEventTime = readLoadTime;
        
        numOfAnswer = 0;
        numOfAlignedRead = 0;
        numOfUnAlignedPairs = 0;
        uint origNumQueries = numQueries;
        // ini_params.Ini_skipSOAP3Alignment = 0;
        if ( detected_read_length == 0 )
        {
            // printParameters(input_options, ini_params);

            // ==================================================================
            // | DETECT THE READ LENGTH                                         |
            // ==================================================================
            if ( input_options.readType == PAIR_END_READ )
            {
                detected_read_length = GetReadLength ( readLengths, numQueries, 2 );
                detected_read_length2 = GetReadLength ( readLengths+1, numQueries, 2 );
                // the minimum insert size cannot be smaller than detected_read_length2
                if (input_options.insert_low < detected_read_length2)
                    input_options.insert_low = detected_read_length2;
                if (input_options.insert_low < detected_read_length)
                    input_options.insert_low = detected_read_length;
            }
            else
            {
                detected_read_length = GetReadLength ( readLengths, numQueries, 1 );
            }

            // ==================================================================
            // | FOR DP IS ENABLED                                              |
            // | IF READ LENGTH < MIN_READ_LEN_FOR_DP (i.e. 30)                 |
            // |    THEN DP IS DISABLE.,                                        |
            // | IF READ LENGTH > 150, SKIP SOAP3 MODULE.                       |
            // | IF READ LENGTH <= 50,  ONLY ALLOW 1 MISMATCH IN SOAP3          |
            // ==================================================================
            if ( input_options.readType == SINGLE_READ && detected_read_length < MIN_READ_LEN_FOR_DP && input_options.enableDP == 1 )
            {
                input_options.enableDP = 0;
                fprintf (stderr, "Dynamic programming is disabled because read length < %i\n", MIN_READ_LEN_FOR_DP );
            }
            else if ( input_options.readType == PAIR_END_READ && input_options.enableDP == 1 )
            {
                //ini_params.Ini_skipSOAP3Alignment = 0;
                fprintf (stderr, "All reads are directly processed by DP\n" );
            }
            else if ( input_options.readType == PAIR_END_READ && ( detected_read_length <= 50 || detected_read_length2 <= 50 ) && input_options.enableDP == 1 )
            {
                input_options.numMismatch = 1;
            }
            
            // ==================================================================
            // | IF DP IS DISABLE AND USER DOES NOT SPECIFY # OF MISMATCHES     |
            // | THEN SET THE DEFAULT # OF MISMATCHES AS:                       |
            // |   - IF READ LENGTH < 50, DEFAULT_NUM_MISMATCH_FOR_SHORT_READ   |
            // |   - IF READ LENGTH >= 50, DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ |
            // ==================================================================
            if ( input_options.enableDP == 0 && input_options.numMismatch == -1 )
            {
                // user does not specify # of mismatches
                input_options.numMismatch = getDefaultMismatchNum ( detected_read_length );
                fprintf (stderr, "Maximum number of mismatches allowed: %i\n",  input_options.numMismatch );
            }

            // get the max hit # for default DP
            if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
            {
                input_options.maxHitNum = getMaxHitNumForDefaultDP ( detected_read_length );
                input_options.maxHitNum2 = getMaxHitNumForDefaultDP ( detected_read_length2 );
            }
        }
        // Reset the array for storing the single alignment for those unaligned paired-end reads
        memset ( hspaux->x0_array, 0, roundUp * sizeof ( int ) );
        memset ( hspaux->x1_array, 0, roundUp * sizeof ( int ) );
        memset ( hspaux->mismatch_array, 0, roundUp * sizeof ( int ) );
        resetAllHits ( ( AllHits * ) hspaux->allHits );

        // =====================================================================
        // | Reset the arrays for storing alignment results for semi-global DP |
        // =====================================================================

        if ( input_options.enableDP == 1 )
        {
            resetBothUnalignedPairsArrays ( bothUnalignedPairsArrays );
        }

        //*******************************//
        // Hot fix: remove comments      //
        //*******************************//
        if (input_options.ignoreComments) {
            for (int i = 0; i < numQueries; ++i) {
                if (queryComments[i]) {
                    free(queryComments[i]);
                    queryComments[i] = NULL;
                }
            }
        }

        //*******************************//
        // Perform Alignment             //
        //*******************************//


        assert ( input_options.readType == PAIR_END_READ );
        soap3_dp_pair_align ( queries, readLengths, input_options.numMismatch, wordPerQuery,
                          maxBatchSize, numQueries, accumReadNum,
                          index,
                          ini_params, input_options,
                          maxReadLength, detected_read_length, detected_read_length2,
                          upkdQualities,
                          readIDs, queryNames, queryComments,
                          dpInfoForReads, samOutputFilePtr, samOutputDPFilePtr,
                          samOutputUnpairFilePtr,
                          numOfAnswer, numOfAlignedRead,
                          bothUnalignedPairsArrays,
                          seedPool,
                          startTime, lastEventTime, totalAlignmentTime);
        resetSeedPool(seedPool);

        // ========================================
        // | GET NEXT BATCH OF READ               |
        // ========================================
        totalReadsAlignedForInputReads += numOfAlignedRead;
        totalAnsForInputReads += numOfAnswer;
        accumReadNum += origNumQueries;
        ResetBufferStatusToUnfilled ( aiob );
        readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
        queries = readyReadsBuffer->queries;
        readLengths = readyReadsBuffer->readLengths;
        readIDs = readyReadsBuffer->readIDs;
        upkdQualities = readyReadsBuffer->upkdQualities;
        queryNames = readyReadsBuffer->queryNames;
        queryComments = readyReadsBuffer->queryComments;
        isFastq = readyReadsBuffer->isFastq;
        numQueries = readyReadsBuffer->filledNum;

        // If the current opened file still have queries returned
        // skip the code to process the result / open next file; and continue 
        if ( numQueries > 0 )
        {
            // Skip
            continue;
        }
    
        // show the summary of the result and then load another pair of read files
        // ======================================================================================
        // | SHOW THE SUMMARY                                                                   |
        // ======================================================================================
        if ( input_options.readType == PAIR_END_READ )
        {
            // fprintf(stderr,"[Main] Overall number of pairs of reads aligned: %llu (number of alignments: %llu)\n", totalReadsAlignedForInputReads/2, totalAnsForInputReads);
            fprintf (stderr, "[Main] Overall number of pairs of reads aligned: %llu\n", totalReadsAlignedForInputReads / 2 );
        }
        else
        {
            // fprintf("[Main] Overall number of reads aligned: %llu (number of alignments: %llu)\n", totalReadsAlignedForInputReads, totalAnsForInputReads);
            // fprintf("[Main] Overall number of unaligned reads: %llu\n", accumReadNum-totalReadsAlignedForInputReads);
            fprintf (stderr, "[Main] Overall number of reads aligned: %llu\n", totalReadsAlignedForInputReads );
            fprintf (stderr, "[Main] Overall number of unaligned reads: %llu\n", accumReadNum - totalReadsAlignedForInputReads );
        }

        fprintf (stderr, "[Main] Overall read load time : %9.4f seconds\n", totalReadLoadTime );
        fprintf (stderr, "[Main] Overall alignment time (excl. read loading) : %9.4f seconds\n", totalAlignmentTime + totalTrimAlignmentTime );

        // update the output files
        if (samOutputDPFilePtr) samclose ( samOutputDPFilePtr );
        for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ ) if (samOutputFilePtr[i]) samclose ( samOutputFilePtr[i] );
        if (samOutputUnpairFilePtr) samclose ( samOutputUnpairFilePtr );

        // release memory and close files
        for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
        { free ( outputFileName[i] ); }

        free ( outputDPFileName );
        free ( outputUnpairFileName );

        if ( input_options.isReadBAM )
        {
            bam_close ( bamQueryFile );
            bam_destroy1 ( bam );
        }
        else
        {
            gzclose ( gzQueryFile );

            if ( input_options.readType == PAIR_END_READ )
            {
                gzclose ( gzQueryFile2 );
            }
        }

        // ======================================================================================
        // | Load the next set of files                                                         |
        // ======================================================================================

        if ( input_options.isReadList == 1 && currExp < expNum )
        {
            fprintf (stderr, "\n[Main] Load the next set of files\n" );
            updateInputOption ( ( &input_options ), multiInput, currExp++ );
            queryFileName = input_options.queryFileName;

            if ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 )
            {
                queryFileName2 = input_options.queryFileName2;
            }

            // Update the HSP information
            hspaux->readGroup = input_options.readGroup;

            if ( strlen ( hspaux->readGroup ) == 0 )
            { hspaux->readGroup = input_options.queryFileName; }

            hspaux->sampleName = input_options.sampleName;

            if ( strlen ( hspaux->sampleName ) == 0 )
            { hspaux->sampleName = DEFAULT_SAMPLE_NAME; }

            hspaux->readGrpOption = input_options.readGrpOption;

            // update the output files
            for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
            { outputFileName[i] = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 10 ); }

            outputDPFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 9 );
            outputUnpairFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 8 );

            for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
            {
                sprintf ( outputFileName[i], "%s.gout.%d", input_options.outputPrefix, i + 1 );
            }

            sprintf ( outputDPFileName, "%s.dpout.1", input_options.outputPrefix );
            sprintf ( outputUnpairFileName, "%s.unpair", input_options.outputPrefix );

            if (input_options.outputBAM) {
                SAMOutputHeaderConstruct ( &samOutputHeader, hsp, hspaux, maxReadLength );
                for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
                    assert(samOutputFilePtr[i] = samopen ( outputFileName[i], "wb", &samOutputHeader ));
                assert(samOutputDPFilePtr = samopen ( outputDPFileName, "wb", &samOutputHeader ));
                // For Unpaired reads
                assert(samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wb", &samOutputHeader ));
            }
            // reset the variables
            accumReadNum = 0;
            totalReadsAlignedForInputReads = 0;
            totalAnsForInputReads = 0;
            detected_read_length = 0; // the max read length for the first ten reads
            totalAlignmentTime = 0;
            totalReadLoadTime = 0;

            // load reads
            if ( input_options.readType == SINGLE_READ )
            { fprintf (stderr, "\n\n[Main] Loading read file %s\n", queryFileName ); }
            else if ( input_options.isReadBAM == 0 )
            { fprintf (stderr, "\n\n[Main] Loading read files %s and %s\n", queryFileName, queryFileName2 ); }

            if ( input_options.isReadBAM )
            {
                // the query file is in BAM format
                bamQueryFile = bam_open ( queryFileName, "r" );
                bamHeader = bam_header_init();
                bamHeader = bam_header_read ( bamQueryFile );
                bam = bam_init1();
            }
            else
            {
                gzQueryFile = gzopen ( queryFileName, "r" );
                if ( gzQueryFile == NULL ) { fprintf ( stderr, "Cannot open queryFile\n" ); exit ( 1 );}
                if (ks1 != NULL) { kseq_destroy(ks1); }
                ks1 = kseq_init(gzQueryFile);

                if ( input_options.readType == PAIR_END_READ )
                {
                    // pair-ended reads
                    gzQueryFile2 = gzopen ( queryFileName2, "r" );
                    if ( gzQueryFile2 == NULL ) { fprintf ( stderr, "Cannot open queryFile2\n" ); exit ( 1 );}
                    if (ks2 != NULL) { kseq_destroy(ks2); }
                    ks2 = kseq_init(gzQueryFile2);
                }
            }

            if ( input_options.isReadBAM )
            {
                InputFilePointersSetBam ( ifp, bamQueryFile, bamHeader, bam );
            }
            else if ( input_options.readType == SINGLE_READ )
            {
                InputFilePointersSetSingle ( ifp, gzQueryFile );
            }
            else
            {
                InputFilePointersSetPair ( ifp, ks1, ks2 );
            }

            aiob->reads = ifp;
            //clear buffer' status
            AIOInputBufferClear ( aiob );
            //create io thread
            AIOInputThreadCreate ( aiob, bufferSize, bufferIndex, index->charMap, queryChar,
                                   queryFileBuffer, bufferSize2, bufferIndex2, queryChar2,
                                   queryFileBuffer2 );
            //load reads
            readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
            queries = readyReadsBuffer->queries;
            readLengths = readyReadsBuffer->readLengths;
            readIDs = readyReadsBuffer->readIDs;
            upkdQualities = readyReadsBuffer->upkdQualities;
            queryNames = readyReadsBuffer->queryNames;
            queryComments = readyReadsBuffer->queryComments;
            isFastq = readyReadsBuffer->isFastq;
            numQueries = readyReadsBuffer->filledNum;
            
            // Update the isFastq variable inside HSP
            hspaux->isFastq = isFastq;
        }
    }

    // ======================================================================================
    // | CLEAN UP                                                                           |
    // ======================================================================================
    fprintf(stderr, "Number of reads' seed alignment free: %u\n",freeSeedPool ( seedPool ));

    freeBothUnalignedPairsArrays ( ( BothUnalignedPairsArrays * ) hspaux->readsIDForSingleDP );
    releaseAllHits ( ( AllHits * ) hspaux->allHits );

    fprintf (stderr, "[Main] Free index from host memory..\n" );
    INDEXFree ( index, ini_params.Ini_shareIndex );
    fprintf (stderr, "[Main] Free query & comments..\n");
    for (int i = 0; i < roundUp; ++i) {
        if (queryNames0[i]) free(queryNames0[i]);
        if (queryNames1[i]) free(queryNames1[i]);
        if (queryComments0[i]) free(queryComments0[i]);
        if (queryComments1[i]) free(queryComments1[i]);
    }

    fprintf (stderr, "[Main] Free host memory..\n" );
    
    free ( hspaux->x0_array );
    free ( hspaux->x1_array );
    free ( hspaux->mismatch_array );
    free ( queries0 );
    free ( readLengths0 );
    free ( readIDs0 );
    free ( upkdQualities0 );
    free ( queryNames0 );
    free ( queryComments0 );

    fprintf (stderr, "[Main] freeDPInfoForReads..\n" );
    freeDPInfoForReads ( dpInfoForReads0 );
    fprintf (stderr, "[Main] freeDPInfoForReads.. Done\n" );

    free ( queries1 );
    free ( readLengths1 );
    free ( readIDs1 );
    free ( upkdQualities1 );
    free ( queryNames1 );
    free ( queryComments1 );
    AIOInputBufferFree ( aiob );

    freeBothUnalignedPairsArrays ( bothUnalignedPairsArrays ); // for both-unaligned pairs OR single-unaligned reads
    freeIniParamsSeedingProperties ( ini_params );
    freeIniParamsFilenames ( ini_params );
    if ( multiInput != NULL )
    {
        free ( multiInput );
    }

    if (input_options.outputBAM) SAMOutputHeaderDestruct ( &samOutputHeader );

    fprintf (stderr, "[Main] Overall running time: %f\n", getElapsedTime( startTime ) );

    return 0;
}

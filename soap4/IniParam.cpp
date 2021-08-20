/*
 *
 *    IniParam.c
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

 // For file checking
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "IniParam.h"


// This function is to load the input file for multi mode
MultiInputItem * loadMultiInputFile ( char * fileName, int isPair, int isReadBAM, int & lineNum )
{
    // return the array of MultiInputItem, and
    // the number of lines
    // get the number of '\n' inside the file
    lineNum = FUGetNumOfLines ( fileName );
    // create the array
    MultiInputItem * inputItemsArray = ( MultiInputItem * ) malloc ( ( lineNum + 1 ) * sizeof ( MultiInputItem ) );
    memset ( inputItemsArray, 0, ( lineNum + 1 ) * sizeof ( MultiInputItem ) );
    // read the file
    FILE * filein;
    char buffer[INPUT_BUFFER_SIZE];
    filein = ( FILE * ) fopen ( fileName, "r" );
    size_t bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
    int bufferIndex = 0;
    int currline = 0;
    int isEndOfLine = 0;
    char * tmpStr = ( char * ) malloc ( sizeof ( char ) * MAX_FIELD_LEN );

    while ( bufferSize > 0 && currline < lineNum )
    {
        // not end of file
        // get the query file 1
        FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                         inputItemsArray[currline].queryFile1,
                         isEndOfLine );

        if ( strlen ( inputItemsArray[currline].queryFile1 ) == 0 )
        {
            fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the name of query file 1.\n", fileName );
            exit ( 1 );
        }

        if ( isPair )
        {
            if ( isReadBAM == 0 )
            {
                // get the query file 2
                if ( !isEndOfLine )
                {
                    FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                     inputItemsArray[currline].queryFile2,
                                     isEndOfLine );
                }

                if ( strlen ( inputItemsArray[currline].queryFile2 ) == 0 )
                {
                    fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the name of query file 2.\n", fileName );
                    exit ( 1 );
                }
            }

            // get minimum value of insert size
            memset ( tmpStr, 0, sizeof ( char ) * MAX_FIELD_LEN ); // reset the array

            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 tmpStr, isEndOfLine );
            }

            if ( strlen ( tmpStr ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the minimum value of insert size.\n", fileName );
                exit ( 1 );
            }

            inputItemsArray[currline].insert_low = atoi ( tmpStr );
            // get maximum value of insert size
            memset ( tmpStr, 0, sizeof ( char ) * MAX_FIELD_LEN ); // reset the array

            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 tmpStr, isEndOfLine );
            }

            if ( strlen ( tmpStr ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the maximum value of insert size.\n", fileName );
                exit ( 1 );
            }

            inputItemsArray[currline].insert_high = atoi ( tmpStr );
        }

        // get the output prefix
        if ( !isEndOfLine )
        {
            FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                             inputItemsArray[currline].outputPrefix,
                             isEndOfLine );
        }

        if ( strlen ( inputItemsArray[currline].outputPrefix ) == 0 )
        {
            fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the output prefix.\n", fileName );
            exit ( 1 );
        }

        if ( !isEndOfLine )
        {
            // get the readGrpID
            FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                             inputItemsArray[currline].readGrpID,
                             isEndOfLine );

            if ( strlen ( inputItemsArray[currline].readGrpID ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which contains an empty field for read group ID.\n", fileName );
                exit ( 1 );
            }

            // get the sampleName
            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 inputItemsArray[currline].sampleName,
                                 isEndOfLine );
            }

            if ( strlen ( inputItemsArray[currline].sampleName ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which contains read group ID but does not indicate the sample name.\n", fileName );
                exit ( 1 );
            }

            updateUserInputText ( inputItemsArray[currline].sampleName );

            if ( !isEndOfLine )
            {
                // get the readGrpOpt
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 inputItemsArray[currline].readGrpOpt,
                                 isEndOfLine );

                if ( strlen ( inputItemsArray[currline].readGrpOpt ) == 0 )
                {
                    fprintf ( stderr, "Error. Inside %s, there exists a line which contains an empty field for read group other option.\n", fileName );
                    exit ( 1 );
                }

                updateUserInputText ( inputItemsArray[currline].readGrpOpt );
            }
        }

        if ( !isEndOfLine )
        {
            // skip all characters until the end of line
            FUSkipToEOL ( filein, buffer, bufferSize, bufferIndex );
        }

        currline++;
    }

    if ( currline == 0 )
    {
        fprintf ( stderr, "Error. The file %s is empty!", fileName );
        exit ( 1 );
    }

    fclose ( filein );
    return inputItemsArray;
}

// This function replace the special character \n or \t with a newline or a tab
void updateUserInputText ( char * str )
{
    // "\t" -> '\t'
    // "\n" => '\n'
    int i, j;

    for ( i = 0; i < strlen ( str ) - 1; i++ )
    {
        if ( str[i] == '\\' )
        {
            if ( str[i + 1] == 't' )
            {
                str[i] = '\t';

                // move all the remaining characters one-pace forward
                for ( j = i + 1; j < strlen ( str ) - 1; j++ )
                {
                    str[j] = str[j + 1];
                }

                str[strlen ( str ) - 1] = '\0';
            }
            else if ( str[i + 1] == 'n' )
            {
                str[i] = '\n';

                // move all the remaining characters one-pace forward
                for ( j = i + 1; j < strlen ( str ) - 1; j++ )
                {
                    str[j] = str[j + 1];
                }

                str[strlen ( str ) - 1] = '\0';
            }
        }
    }
}

// This function is to parse the ini file and collect the paramaters of
// 1. The file extension of SA file
// 2. The number of CPU threads
// 3. The alignment model
// 4. The memory size in the GPU card
int ParseIniFile ( char * iniFileName, IniParams & ini_params )
{
    // puts("parse ini");
    dictionary * ini;
    char tmp[10];
    ini = iniparser_load ( iniFileName, FALSE );

    if ( ini == NULL )
    {
        fprintf (stderr, "[ParseIniFile] File not found!\n" );
        return -1;
    }

    // Database:SaValueFileExt parameter
    iniparser_copystring ( ini, "Database:SaValueFileExt", ini_params.Ini_SaValueFileExt, ".sa", MAX_FILEEXT_LEN );
    // Alignment:NumOfCpuThreads parameter
    ini_params.Ini_NumOfCpuThreads = iniparser_getuint ( ini, "Alignment:NumOfCpuThreads", 3 );
    if ( ini_params.Ini_NumOfCpuThreads <= 0 || ini_params.Ini_NumOfCpuThreads > MAX_NUM_CPU_THREADS )
    {
        fprintf (stderr, "[ParseIniFile] Invalid value for Ini_NumOfCpuThreads(%u)!\n", ini_params.Ini_NumOfCpuThreads );
        return -1;
    }

    // Alignment:Soap3MisMatchAllow parameter
    ini_params.Ini_Soap3MisMatchAllow = iniparser_getint ( ini, "Alignment:Soap3MisMatchAllow", 2 );

    if ( ini_params.Ini_Soap3MisMatchAllow < 0 || ini_params.Ini_Soap3MisMatchAllow > 4 )
    {
        fprintf (stderr, "[ParseIniFile] The program only supports number of mismatches for soap3 alignment between 0 and 4.\n" );
        return -1;
    }

    // DP:DPScoreThreshold parameter
    ini_params.Ini_DPScoreThreshold = iniparser_getdouble( ini, "DP:DPScoreThreshold", 0.3 );

    ini_params.Ini_MaxOutputPerRead = iniparser_getuint ( ini, "Alignment:MaxOutputPerRead", 1000 );

    iniparser_copystring ( ini, "PairEnd:StrandArrangement", tmp, "+/-", 10 );

    if ( strcmp ( tmp, "+/+" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
    }
    else if ( strcmp ( tmp, "-/+" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_NEG_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
    }
    else if ( strcmp ( tmp, "-/-" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_NEG_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_NEG_STRAND;
    }
    else
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_NEG_STRAND;
    }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
    ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
    ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
#endif

    ini_params.Ini_PEMaxOutputPerPair = iniparser_getuint ( ini, "PairEnd:MaxOutputPerPair", 1000 );
    ini_params.Ini_MaxHitsEachEndForPairing = iniparser_getuint ( ini, "PairEnd:MaxHitsEachEndForPairing", 10000 );
   
    // iniparser_copystring ( ini, "" , ini_params.Ini_IndelDbFilename);
    
    // FOR DP MODULE
    ini_params.Ini_MatchScore = iniparser_getint ( ini, "DP:MatchScore", 1 );
    ini_params.Ini_MismatchScore = iniparser_getint ( ini, "DP:MismatchScore", -2 );
    ini_params.Ini_GapOpenScore = iniparser_getint ( ini, "DP:GapOpenScore", -3 );
    ini_params.Ini_GapExtendScore = iniparser_getint ( ini, "DP:GapExtendScore", -1 );
    // score
    ini_params.Ini_minMAPQ = iniparser_getint ( ini, "Score:MinMAPQ", 1 );
    ini_params.Ini_maxMAPQ = iniparser_getint ( ini, "Score:MaxMAPQ", 40 );
    // index to be shared among multiple copies of soap3-dp
    ini_params.Ini_shareIndex = iniparser_getint ( ini, "Index:ShareIndex", 0 );
    // max length allowed for read name
    ini_params.Ini_maxReadNameLen = iniparser_getint ( ini, "Reads:MaxLenReadName", 64 );
    // max length from the front of the read for clipping
    ini_params.Ini_maxFrontLenClipped = iniparser_getint ( ini, "Clipping:MaxFrontLenClipped", 3 );
    // max length from the end of the read for clipping
    ini_params.Ini_maxEndLenClipped = iniparser_getint ( ini, "Clipping:MaxEndLenClipped", 8 );
    // whether the seed will proceed to perform DP if there are too many hits
    ini_params.Ini_proceedDPForTooManyHits = iniparser_getint ( ini, "OtherSettings:ProceedDPForTooManyHits", 0 );
    // whether the read will perform SOAP3 module
    ini_params.Ini_skipSOAP3Alignment = iniparser_getint ( ini, "OtherSettings:SkipSOAP3Alignment", 0 );
    // whether the bwa-like MAPQ score should be reported
    ini_params.Ini_bwaLikeScore = iniparser_getint ( ini, "Score:BWALikeScore", 0 );
    // whether performing Single-end DP
    ini_params.Ini_skipSingleEndDP = iniparser_getint(ini, "OtherSettings:SkipSingleEndDP", 0);
    // whether performing Default DP
    ini_params.Ini_skipDefaultDP = iniparser_getint(ini, "OtherSettings:SkipDefaultDP", 0);
    // Single-end DP score threshold
    ini_params.Ini_SingleEndDPScoreThreshold = iniparser_getdouble(ini, "SingleEndDP:DPScoreThreshold", 0.3);
    // Default DP score threshold
    ini_params.Ini_DefaultDPScoreThreshold = iniparser_getdouble(ini, "DefaultDP:DPScoreThreshold", 0.3);

    // seeding properties
    ini_params.Ini_numberOfRoundOfDeepDpForShortReads = iniparser_getint(ini, "DP:NumberOfRoundOfDeepDPForShortReads", 0);
    if (ini_params.Ini_numberOfRoundOfDeepDpForShortReads > 0)
        ini_params.seedingPropertiesForShortReads = (SeedingProperties*)malloc(ini_params.Ini_numberOfRoundOfDeepDpForShortReads*sizeof(SeedingProperties));
    else
        ini_params.seedingPropertiesForShortReads = NULL;
    for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForShortReads;++i){
        char buffer[100];
        snprintf(buffer,100,"SeedingRound%dForShortReads:DPScoreThreshold",i+1);
        ini_params.seedingPropertiesForShortReads[i].dpScoreThreshold = iniparser_getdouble(ini, buffer, -1);
        snprintf(buffer,100,"SeedingRound%dForShortReads:mmpSeedingScheme",i+1);
        ini_params.seedingPropertiesForShortReads[i].mmpSeedingScheme = iniparser_getint(ini, buffer, 0);
        // construct buffer to be "Seeding:NumberOfSeedForDeepDPRound${i}"
        snprintf(buffer,100,"SeedingRound%dForShortReads:NumberOfSeedForDeepDPRound",i+1);
        int numberOfSeed = ini_params.seedingPropertiesForShortReads[i].numberOfSeed = iniparser_getint(ini, buffer, 0);
        ini_params.seedingPropertiesForShortReads[i].seedPos = (int *) malloc(numberOfSeed*sizeof(int));
        ini_params.seedingPropertiesForShortReads[i].minSeedLength = (int *) malloc(numberOfSeed*sizeof(int));
        ini_params.seedingPropertiesForShortReads[i].maxSeedLength = (int *) malloc(numberOfSeed*sizeof(int));

        for (int j=0;j<numberOfSeed;++j){
            // construct buffer to be "SeedingRound${i}:SeedPosition${j}"
            snprintf(buffer,100,"SeedingRound%dForShortReads:SeedPosition%d",i+1,j+1);
            ini_params.seedingPropertiesForShortReads[i].seedPos[j] = iniparser_getint(ini, buffer, 0);
            // construct buffer to be "SeedingRound${i}:SeedMinLength${j}"
            snprintf(buffer,100,"SeedingRound%dForShortReads:SeedMinLength%d",i+1,j+1);
            ini_params.seedingPropertiesForShortReads[i].minSeedLength[j] = iniparser_getint(ini, buffer, 0);
            // construct buffer to be "SeedingRound${i}:SeedMaxLength${j}"
            snprintf(buffer,100,"SeedingRound%dForShortReads:SeedMaxLength%d",i+1,j+1);
            ini_params.seedingPropertiesForShortReads[i].maxSeedLength[j] = iniparser_getint(ini, buffer, 0);
            snprintf(buffer,100,"SeedingRound%dForShortReads:MaxSeedHit",i+1);
            ini_params.seedingPropertiesForShortReads[i].maxNumberOfSeedHit = iniparser_getint(ini, buffer, MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ);
        }
        snprintf(buffer,100,"SeedingRound%dForShortReads:NumNBMismatch",i+1);
        ini_params.seedingPropertiesForShortReads[i].numNBMismatch = iniparser_getint(ini, buffer, SEEDING_NUM_NBM);
        snprintf(buffer,100,"SeedingRound%dForShortReads:RestrictedNBM",i+1);
        ini_params.seedingPropertiesForShortReads[i].restrictedNBM = iniparser_getint(ini, buffer, 0);
        snprintf(buffer,100,"SeedingRound%dForShortReads:ProceedDPForTooManyHits",i+1);
        ini_params.seedingPropertiesForShortReads[i].proceedDPForTooManyHits = iniparser_getint(ini, buffer, 0);
        snprintf(buffer,100,"SeedingRound%dForShortReads:StoreSeedAlignments",i+1);
        ini_params.seedingPropertiesForShortReads[i].isStoreSeedAlignments = iniparser_getint(ini, buffer, 1);
    }

    // seeding properties
    ini_params.Ini_numberOfRoundOfDeepDpForLongReads = iniparser_getint(ini, "DP:NumberOfRoundOfDeepDPForLongReads", 0);
    if (ini_params.Ini_numberOfRoundOfDeepDpForLongReads > 0)
        ini_params.seedingPropertiesForLongReads = (SeedingProperties*)malloc(ini_params.Ini_numberOfRoundOfDeepDpForLongReads*sizeof(SeedingProperties));
    else
        ini_params.seedingPropertiesForLongReads = NULL;
    for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForLongReads;++i){
        char buffer[100];
        snprintf(buffer,100,"SeedingRound%dForLongReads:DPScoreThreshold",i+1);
        ini_params.seedingPropertiesForLongReads[i].dpScoreThreshold = iniparser_getdouble(ini, buffer, -1);
        snprintf(buffer,100,"SeedingRound%dForLongReads:mmpSeedingScheme",i+1);
        ini_params.seedingPropertiesForLongReads[i].mmpSeedingScheme = iniparser_getint(ini, buffer, 0);
        // construct buffer to be "Seeding:NumberOfSeedForDeepDPRound${i}"
        snprintf(buffer,100,"SeedingRound%dForLongReads:NumberOfSeedForDeepDPRound",i+1);
        int numberOfSeed = ini_params.seedingPropertiesForLongReads[i].numberOfSeed = iniparser_getint(ini, buffer, SEEDING_NUM_NBM);
        ini_params.seedingPropertiesForLongReads[i].seedPos = (int *) malloc(numberOfSeed*sizeof(int));
        ini_params.seedingPropertiesForLongReads[i].minSeedLength = (int *) malloc(numberOfSeed*sizeof(int));
        ini_params.seedingPropertiesForLongReads[i].maxSeedLength = (int *) malloc(numberOfSeed*sizeof(int));

        for (int j=0;j<numberOfSeed;++j){
            // construct buffer to be "SeedingRound${i}:SeedPosition${j}"
            snprintf(buffer,100,"SeedingRound%dForLongReads:SeedPosition%d",i+1,j+1);
            ini_params.seedingPropertiesForLongReads[i].seedPos[j] = iniparser_getint(ini, buffer, 0);
            // construct buffer to be "SeedingRound${i}:SeedMinLength${j}"
            snprintf(buffer,100,"SeedingRound%dForLongReads:SeedMinLength%d",i+1,j+1);
            ini_params.seedingPropertiesForLongReads[i].minSeedLength[j] = iniparser_getint(ini, buffer, 0);
            // construct buffer to be "SeedingRound${i}:SeedMaxLength${j}"
            snprintf(buffer,100,"SeedingRound%dForLongReads:SeedMaxLength%d",i+1,j+1);
            ini_params.seedingPropertiesForLongReads[i].maxSeedLength[j] = iniparser_getint(ini, buffer, 0);
            snprintf(buffer,100,"SeedingRound%dForLongReads:MaxSeedHit",i+1);
            ini_params.seedingPropertiesForLongReads[i].maxNumberOfSeedHit = iniparser_getint(ini, buffer, MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ);
        }
        snprintf(buffer,100,"SeedingRound%dForLongReads:NumNBMismatch",i+1);
        ini_params.seedingPropertiesForLongReads[i].numNBMismatch = iniparser_getint(ini, buffer, SEEDING_NUM_NBM);
        snprintf(buffer,100,"SeedingRound%dForLongReads:RestrictedNBM",i+1);
        ini_params.seedingPropertiesForLongReads[i].restrictedNBM = iniparser_getint(ini, buffer, 0);
        snprintf(buffer,100,"SeedingRound%dForLongReads:ProceedDPForTooManyHits",i+1);
        ini_params.seedingPropertiesForLongReads[i].proceedDPForTooManyHits = iniparser_getint(ini, buffer, 0);
        snprintf(buffer,100,"SeedingRound%dForLongReads:StoreSeedAlignments",i+1);
        ini_params.seedingPropertiesForLongReads[i].isStoreSeedAlignments = iniparser_getint(ini, buffer, 1);
    }
    
    ini_params.mmpProperties.seedSAsizeThreshold = iniparser_getint(ini, "MMP:mmpSeedSAsizeThreshold", 30);
    ini_params.mmpProperties.seedMinLength = iniparser_getint(ini, "MMP:mmpSeedMinLength", 17);
    ini_params.mmpProperties.uniqThreshold = iniparser_getint(ini, "MMP:mmpUniqThreshold", 6);
    ini_params.mmpProperties.indelFuzz = iniparser_getint(ini, "MMP:mmpIndelFuzz", 5);
    ini_params.mmpProperties.goodSeedLen = iniparser_getint(ini, "MMP:mmpGoodSeedLen", 27);
    ini_params.mmpProperties.reseedLen = iniparser_getint(ini, "MMP:mmpReseedLen", 18);
    ini_params.mmpProperties.reseedRLTratio = iniparser_getdouble(ini, "MMP:mmpReseedRLTratio", 0.85);
    ini_params.mmpProperties.reseedAbsDiff = iniparser_getint(ini, "MMP:mmpReseedAbsDiff", 4);
    ini_params.mmpProperties.shortSeedRatio = iniparser_getdouble(ini, "MMP:mmpShortSeedRatio", 0.5);
    iniparser_freedict ( ini );
    return 0;
}


void INIPrintAlignmentUsage ( InputOptions * inputOptions, char * programName )
{
    if ( inputOptions->readType == SINGLE_READ )
    {
        if ( !inputOptions->isReadList )
        {
            UIprintUsageSingle ( programName );
        }
        else
        {
            UIprintUsageSingleList ( programName );
        }
    }
    else if ( inputOptions->readType == PAIR_END_READ )
    {
        if ( !inputOptions->isReadList )
        {
            UIprintUsagePair ( programName );
        }
        else
        {
            UIprintUsagePairList ( programName );
        }
    }
    else
    {
        UIPrintUsageOverview ( programName );
    }
}

void INIPrintBalsaUsage ( InputOptions * inputOptions, char * programName )
{
    if ( inputOptions->readType == PAIR_END_READ )
    {
        if ( !inputOptions->isReadList )
        {
            UIprintBalsaUsagePair ( programName );
        }
        else
        {
            UIprintBalsaUsagePairList ( programName );
        }
    }
    else
    {
        UIPrintUsageOverview ( programName );
    }
}


void checkExistIni ( const char * path, char isDir )
{
    if ( !path )
    {
        fprintf ( stderr, "Error : path is NULL\n" );
        exit ( 1 );
    }
    struct stat s;
    int err = stat ( path, &s );
    if ( -1 == err )
    {
        fprintf ( stderr, "Error : %s does not exist! Please check your ini setting or input arguments.\n", path );
        exit ( 1 );
    }
    else
    {
        if ( isDir && ! S_ISDIR ( s.st_mode ) )
        {
            fprintf ( stderr, "Error : %s is not a directory! Please check your ini setting or input arguments.\n", path );
            exit ( 1 );
        }
        if ( !isDir && S_ISDIR ( s.st_mode ) )
        {
            fprintf ( stderr, "Error : %s is a directory! Please check your ini setting or input arguments.\n", path );
            exit ( 1 );
        }
    }
}

void checkValidOutputPrefix ( const char * path )
{
    if ( !path )
    {
        fprintf ( stderr, "Error : path is NULL\n" );
        exit ( 1 );
    }
    if ( path[strlen ( path ) - 1] == '/' )
    {
        fprintf ( stderr, "Error : %s should be a prefix not a directory\n", path );
        exit ( 1 );
    }
    FILE * testPath = ( FILE * ) fopen ( path, "w" );
    if ( !testPath )
    {
        fprintf ( stderr, "Error : %s not accessible\n", path );
        exit ( 1 );
    }
    fclose ( testPath );
    remove ( path );
}

bool parseInputArgs ( int argc, char ** argv, InputOptions & input_options )
{
    char * programName = argv[0];
    bool readGroupSupplied = false; // whether an option "-D" is specified
    bool sampleNameSupplied = false; // whether an option "-A" is specified
    bool readGroupOptionSupplied = false; // whether an option "-R" is specified

    // This function is to check and parse the input arguments
    // arguments: <program> single bwtCodeIndex queryFileName numQueries maxReadLength [options]
    // OR         <program> pair bwtCodeIndex queryFileName1 queryFileName2 numQueries maxReadLength [options]
    // If the argument is not correct, then it returns FALSE. Else it returns TRUE
    if ( argc < 2 )
    {
        // fprintf(stderr,"No argument has been provided.\n");
        // printUsage(argv[0]);
        UIPrintUsageOverview ( programName );
        return false;
    }

    if ( strcmp ( argv[1], "single" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "pair" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "single-multi" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = TRUE;
    }
    else if ( strcmp ( argv[1], "pair-multi" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = TRUE; // list of read files
    }
    else
    {
        fprintf ( stderr, "Miss the keyword 'single', 'single-multi', 'pair' or 'pair-multi'.\n" );
        return false;
    }

    if ( argc == 2 )
    {
        INIPrintAlignmentUsage ( &input_options, programName );
        return false;
    }

    // check whether the input read type is BAM
    input_options.isReadBAM = 0;

    for ( int i = 2; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-bam" ) == 0 || strcmp ( argv[i], "-BAM" ) == 0 )
        {
            input_options.isReadBAM = 1;
            break;
        }
    }

    // check the number of arguments for the specific read type
    int min_num_args = ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 && input_options.isReadList == 0 ) ? 5 : 4;

    if ( argc < min_num_args )
    {
        fprintf ( stderr, "Invalid number of command-line arguments.\n" );
        INIPrintAlignmentUsage ( &input_options, programName );
        return false;
    }

    input_options.indexName = argv[2];
    // get the query file name(s) and the maximum read length
    input_options.queryFileName = argv[3];

    // get the second query file name
    // only for "pair" mode and not a BAM input
    if ( input_options.readType == PAIR_END_READ && !input_options.isReadList && !input_options.isReadBAM )
    {
        input_options.queryFileName2 = argv[4];
    }

    // set default values
    input_options.maxReadLength = DEFAULT_MAX_READ_LEN;
    input_options.numMismatch = -1;
    input_options.insert_low = 1;
    input_options.insert_high = 500;
    input_options.outputBAM = 0;
    input_options.outputPrefix = argv[3];
    input_options.enableDP = 1;
    input_options.isIlluminaQual = 0;
    input_options.GPUDeviceID = -1;
    input_options.numCpuThreads = 4;

    input_options.readGroup = ( char * ) "";
    input_options.sampleName = ( char * ) "";
    input_options.readGrpOption = ( char * ) ""; // read group option
    input_options.iniFileName = NULL;
    // set default alignment type = all-best alignment
    input_options.alignmentType = OUTPUT_ALL_BEST;
    input_options.isPrintMDNM = false;
    input_options.isExome = 0; 
    input_options.megapathMode = 0;
    input_options.top = 95;
    input_options.ignoreComments = 0;

    // get the options
    for ( int i = min_num_args; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-h" ) == 0 )
        {
            // output option
            // 1: all valid alignments
            // 2: all best alignments (default)
            // 3: unique best alignments
            // 4: random best alignments
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the output option after '-h' ( 1, 2, 3 or 4 )\n" );
                return false;
            }

            int input_type = atoi ( argv[++i] );

            if ( input_type == 1 )
            {
                input_options.alignmentType = OUTPUT_ALL_VALID;
            }
            else if ( input_type == 2 )
            {
                input_options.alignmentType = OUTPUT_ALL_BEST;
            }
            else if ( input_type == 3 )
            {
                input_options.alignmentType = OUTPUT_UNIQUE_BEST;
            }
            else if ( input_type == 4 )
            {
                input_options.alignmentType = OUTPUT_RANDOM_BEST;
            }
            else
            {
                fprintf ( stderr, "The output option should be 1, 2, 3 or 4\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-l" ) == 0 || strcmp ( argv[i], "-L" ) == 0 )
        {
            // the length of the longest read in the input
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the length after '-L'\n" );
                return false;
            }

            input_options.maxReadLength = atoi ( argv[++i] );

            if ( input_options.maxReadLength < 0 )
            {
                fprintf ( stderr, "The length should not be less than 0\n" );
                return false;
            }
            else if ( input_options.maxReadLength > MAX_READ_LENGTH )
            {
                fprintf ( stderr, "The length should not be greater than %u\n", MAX_READ_LENGTH );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-u" ) == 0 )
        {
            // the upper bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the maximum value of insert size after '-u'\n" );
                return false;
            }

            input_options.insert_high = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-v" ) == 0 )
        {
            // the lower bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the minimum value of insert size after '-v'\n" );
                return false;
            }

            input_options.insert_low = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-b" ) == 0 )
        {
            input_options.outputBAM = 1;
        }
        else if ( strcmp ( argv[i], "-s" ) == 0 )
        {
            // disable dynamic programming for the unmapped mate
            input_options.enableDP = 0;

            if ( i + 1 < argc && isdigit ( argv[i + 1][0] ) )
            {
                int tmp = atoi ( argv[i + 1] );

                if ( tmp >= 0 && tmp <= 4 )
                {
                    input_options.numMismatch = tmp;
                    i++;
                }
                else
                {
                    fprintf ( stderr, "The value after '-s' should be between 0 and 4\n" );
                }
            }
        }
        else if ( strcmp ( argv[i], "-o" ) == 0 )
        {
            // specify the output file prefix
            if ( input_options.isReadList == 1 )
            {
                fprintf ( stderr, "For multiple sets of input files, the output file prefix could not be specified." );
                fprintf ( stderr, "The option '-o' is thus ignored." );
            }
            else
            {
                if ( i + 1 >= argc )
                {
                    fprintf ( stderr, "Please specify the output file prefix after '-o'\n" );
                    return false;
                }

                input_options.outputPrefix = argv[++i];
            }
        }
        else if ( strcmp ( argv[i], "-I" ) == 0 )
        {
            // Quality is in Phred+64 format
            input_options.isIlluminaQual = 1;
        }
        else if ( strcmp ( argv[i], "-D" ) == 0 )
        {
            // to assign the read group ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the read group ID after '-D'\n" );
                return false;
            }

            input_options.readGroup = argv[++i];
            readGroupSupplied = true;
        }
        else if ( strcmp ( argv[i], "-A" ) == 0 )
        {
            // to assign the sample name
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-A'\n" );
                return false;
            }

            input_options.sampleName = argv[++i];
            // update the input text
            updateUserInputText ( input_options.sampleName );
            sampleNameSupplied = true;
        }
        else if ( strcmp ( argv[i], "-R" ) == 0 )
        {
            // to assign the read group option
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-R'\n" );
                return false;
            }

            input_options.readGrpOption = argv[++i];
            // update the input text
            updateUserInputText ( input_options.readGrpOption );
            readGroupOptionSupplied = true;
        }
        else if ( strcmp ( argv[i], "-c" ) == 0 )
        {
            // to specific the GPU device ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the GPU device ID after '-c'\n" );
                return false;
            }

            input_options.GPUDeviceID = atoi ( argv[++i] );

            if ( input_options.GPUDeviceID < 0 )
            {
                fprintf ( stderr, "The GPU device ID should not be less than 0\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-p" ) == 0 )
        {
            // print MD string and NM tag
            input_options.isPrintMDNM = true;
        }
        else if ( strcmp ( argv[i], "-e" ) == 0 )
        {
            // to specifiy the Exome Region
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the Exome Region Index after '-e'\n" );
                return false;
            }

            input_options.exomeRegionFileName = argv[++i];
            input_options.isExome = 1;
        }
        else if ( strcmp ( argv[i], "-T" ) == 0 ) {
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the number of CPU threads after '-T'\n" );
                return false;
            }

            input_options.numCpuThreads = atoi ( argv[++i] );

            if ( input_options.numCpuThreads <= 0 ) {
                fprintf ( stderr, "Please specify a positive number of CPU threads after '-T'\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-C" ) == 0 ) {
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify ini file name after '-C'\n" );
                return false;
            }

            input_options.iniFileName = argv[++i];
        }
        else if ( strcmp(argv[i], "-F") == 0) {
            input_options.megapathMode = 1;
        }
        else if ( strcmp(argv[i], "-P") == 0) {
            input_options.megapathMode = 2;
        }
        else if ( strcmp(argv[i], "-top") == 0 ) {
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the value of '-top'\n" );
                return false;
            }
            input_options.top = atoi(argv[i+1]);
        }
        else if ( strcmp(argv[i], "-nc") == 0 ) {
            input_options.ignoreComments = 1;
        }
    }

    if ( input_options.readType == PAIR_END_READ )
    {
        if ( input_options.insert_high == -1 )
        {
            fprintf ( stderr, "Please specify the maximum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low == -1 )
        {
            fprintf ( stderr, "Please specify the minimum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low > input_options.insert_high )
        {
            fprintf ( stderr, "The minimum value of insert size should not be greater than the maximum value of insert size.\n" );
            return false;
        }
    }

    // only allow three cases:
    // case 1: user does not specify "-D", "-A" or "-R"
    // case 2: user only specifies "-D" and "-A"
    // case 3: user only specifies "-D", "-A" and "-R"
    if ( readGroupOptionSupplied && ( ( !sampleNameSupplied ) || ( !readGroupSupplied ) ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, both options '-A' and '-R' have to be specified too.\n" );
        return false;
    }

    if ( sampleNameSupplied && ( !readGroupSupplied ) )
    {
        fprintf ( stderr, "Since an option '-A' is specified, the options '-D' has to be specified too.\n" );
        return false;
    }

    if ( readGroupSupplied && ( !sampleNameSupplied ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, the options '-A' has to be specified too.\n" );
        return false;
    }

    return true;
}

bool parseBalsaInputArgs ( int argc, char ** argv, InputOptions & input_options )
{
    char * programName = argv[0];
    bool readGroupSupplied = false; // whether an option "-D" is specified
    bool sampleNameSupplied = false; // whether an option "-A" is specified
    bool readGroupOptionSupplied = false; // whether an option "-R" is specified

    // This function is to check and parse the input arguments
    // arguments: <program> single bwtCodeIndex queryFileName numQueries maxReadLength [options]
    // OR         <program> pair bwtCodeIndex queryFileName1 queryFileName2 numQueries maxReadLength [options]
    // If the argument is not correct, then it returns FALSE. Else it returns TRUE
    if ( argc < 2 )
    {
        // fprintf(stderr,"No argument has been provided.\n");
        // printUsage(argv[0]);
        UIPrintUsageOverview ( programName );
        return false;
    }

    if ( strcmp ( argv[1], "single" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "pair" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "single-multi" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = TRUE;
    }
    else if ( strcmp ( argv[1], "pair-multi" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = TRUE; // list of read files
    }
    else
    {
        fprintf ( stderr, "Miss the keyword 'single', 'single-multi', 'pair' or 'pair-multi'.\n" );
        return false;
    }

    if ( argc == 2 )
    {
        INIPrintBalsaUsage ( &input_options, programName );
        return false;
    }

    // check whether the input read type is BAM
    input_options.isReadBAM = 0;

    for ( int i = 2; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-bam" ) == 0 || strcmp ( argv[i], "-BAM" ) == 0 )
        {
            input_options.isReadBAM = 1;
            break;
        }
    }

    // check the number of arguments for the specific read type
    int min_num_args = ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 && input_options.isReadList == 0 ) ? 9 : 8;

    if ( argc < min_num_args )
    {
        fprintf ( stderr, "Invalid number of command-line arguments.\n" );
        INIPrintBalsaUsage ( &input_options, programName );
        return false;
    }

    input_options.indexName = argv[2];
    input_options.dbSnpIndexFileName = argv[3];
    char * dbSnpBVname = ( char * ) malloc ( strlen ( input_options.dbSnpIndexFileName ) + 4 );
    char * dbSnpListname = ( char * ) malloc ( strlen ( input_options.dbSnpIndexFileName ) + 7 );
    sprintf ( dbSnpBVname, "%s.bv", input_options.dbSnpIndexFileName );
    sprintf ( dbSnpListname, "%s.1list", input_options.dbSnpIndexFileName );
    checkExistIni ( dbSnpBVname, 0 );
    checkExistIni ( dbSnpListname, 0 );
    free ( dbSnpBVname );
    free ( dbSnpListname );
    input_options.indelDBIndexFileName = argv[4];
    checkExistIni ( input_options.indelDBIndexFileName, 0 );
    input_options.geneRegionFileName = argv[5];
    checkExistIni ( input_options.geneRegionFileName, 0 );
    // get the query file name(s) and the maximum read length
    input_options.queryFileName = argv[6];

    // get the second query file name
    // only for "pair" mode and not a BAM input
    if ( input_options.readType == PAIR_END_READ && !input_options.isReadList && !input_options.isReadBAM )
    {
        input_options.queryFileName2 = argv[7];
        input_options.resultPrefix = argv[8];
    }
    else
    {
        input_options.resultPrefix = argv[7];
    }

    // set default values
    input_options.maxReadLength = DEFAULT_MAX_READ_LEN;
    input_options.numMismatch = -1;
    input_options.insert_low = 1;
    input_options.insert_high = 500;
    input_options.outputBAM = 0;
    input_options.outputPrefix = argv[6];
    input_options.enableDP = 1;
    input_options.isIlluminaQual = 0;
    input_options.GPUDeviceID = -1;
    input_options.readGroup = ( char * ) "";
    input_options.sampleName = ( char * ) "";
    input_options.readGrpOption = ( char * ) ""; // read group option
    // set default alignment type = all-best alignment
    input_options.alignmentType = OUTPUT_ALL_BEST;
    input_options.isPrintMDNM = false;
    input_options.default_tempPrefix = 1;
    input_options.tempPrefix = NULL;
    input_options.isExome = 0; 

    // get the options
    for ( int i = min_num_args; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-h" ) == 0 )
        {
            // output option
            // 1: all valid alignments
            // 2: all best alignments (default)
            // 3: unique best alignments
            // 4: random best alignments
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the output option after '-h' ( 1, 2, 3 or 4 )\n" );
                return false;
            }

            int input_type = atoi ( argv[++i] );

            if ( input_type == 1 )
            {
                input_options.alignmentType = OUTPUT_ALL_VALID;
            }
            else if ( input_type == 2 )
            {
                input_options.alignmentType = OUTPUT_ALL_BEST;
            }
            else if ( input_type == 3 )
            {
                input_options.alignmentType = OUTPUT_UNIQUE_BEST;
            }
            else if ( input_type == 4 )
            {
                input_options.alignmentType = OUTPUT_RANDOM_BEST;
            }
            else
            {
                fprintf ( stderr, "The output option should be 1, 2, 3 or 4\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-l" ) == 0 || strcmp ( argv[i], "-L" ) == 0 )
        {
            // the length of the longest read in the input
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the length after '-L'\n" );
                return false;
            }

            input_options.maxReadLength = atoi ( argv[++i] );

            if ( input_options.maxReadLength < 0 )
            {
                fprintf ( stderr, "The length should not be less than 0\n" );
                return false;
            }
            else if ( input_options.maxReadLength > MAX_READ_LENGTH )
            {
                fprintf ( stderr, "The length should not be greater than %u\n", MAX_READ_LENGTH );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-u" ) == 0 )
        {
            // the upper bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the maximum value of insert size after '-u'\n" );
                return false;
            }

            input_options.insert_high = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-v" ) == 0 )
        {
            // the lower bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the minimum value of insert size after '-v'\n" );
                return false;
            }

            input_options.insert_low = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-b" ) == 0 )
        {
            input_options.outputBAM = 1;
        }
        else if ( strcmp ( argv[i], "-s" ) == 0 )
        {
            // disable dynamic programming for the unmapped mate
            input_options.enableDP = 0;

            if ( i + 1 < argc && isdigit ( argv[i + 1][0] ) )
            {
                int tmp = atoi ( argv[i + 1] );

                if ( tmp >= 0 && tmp <= 4 )
                {
                    input_options.numMismatch = tmp;
                    i++;
                }
                else
                {
                    fprintf ( stderr, "The value after '-s' should be between 0 and 4\n" );
                }
            }
        }
        else if ( strcmp ( argv[i], "-o" ) == 0 )
        {
            // specify the output file prefix
            if ( input_options.isReadList == 1 )
            {
                fprintf ( stderr, "For multiple sets of input files, the output file prefix could not be specified." );
                fprintf ( stderr, "The option '-o' is thus ignored." );
            }
            else
            {
                if ( i + 1 >= argc )
                {
                    fprintf ( stderr, "Please specify the output file prefix after '-o'\n" );
                    return false;
                }

                input_options.outputPrefix = argv[++i];
            }
        }
        else if ( strcmp ( argv[i], "-I" ) == 0 )
        {
            // Quality is in Phred+64 format
            input_options.isIlluminaQual = 1;
        }
        else if ( strcmp ( argv[i], "-D" ) == 0 )
        {
            // to assign the read group ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the read group ID after '-D'\n" );
                return false;
            }

            input_options.readGroup = argv[++i];
            readGroupSupplied = true;
        }
        else if ( strcmp ( argv[i], "-A" ) == 0 )
        {
            // to assign the sample name
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-A'\n" );
                return false;
            }

            input_options.sampleName = argv[++i];
            // update the input text
            updateUserInputText ( input_options.sampleName );
            sampleNameSupplied = true;
        }
        else if ( strcmp ( argv[i], "-R" ) == 0 )
        {
            // to assign the read group option
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-R'\n" );
                return false;
            }

            input_options.readGrpOption = argv[++i];
            // update the input text
            updateUserInputText ( input_options.readGrpOption );
            readGroupOptionSupplied = true;
        }
        else if ( strcmp ( argv[i], "-c" ) == 0 )
        {
            // to specific the GPU device ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the GPU device ID after '-c'\n" );
                return false;
            }

            input_options.GPUDeviceID = atoi ( argv[++i] );

            if ( input_options.GPUDeviceID < 0 )
            {
                fprintf ( stderr, "The GPU device ID should not be less than 0\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-p" ) == 0 )
        {
            // print MD string and NM tag
            input_options.isPrintMDNM = true;
        }
        else if ( strcmp ( argv[i], "-e" ) == 0 )
        {
            // to specifiy the Exome Region
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the Exome Region Index after '-e'\n" );
                return false;
            }

            input_options.exomeRegionFileName = argv[++i];
            input_options.isExome = 1;
        }
        else if ( strcmp (argv[i], "-t") == 0 )
        {
            if ( i + 1 >= argc )
            {
                    fprintf ( stderr, "Please specify the temporary directory after '-t'\n" );
                    return false;
            }

            input_options.tempPrefix = argv[++i];
            checkValidOutputPrefix ( input_options.tempPrefix );
            input_options.default_tempPrefix = 0;
        }
        else if ( strcmp ( argv[i], "-snapshot" ) == 0 )
        {
            input_options.outputSnapshot = 1;
        }
    }

    if ( input_options.tempPrefix == NULL )
    {
        checkValidOutputPrefix ( "./balsa_" );
        input_options.tempPrefix = ( char * ) malloc ( sizeof ( char ) * 7 );
        sprintf ( input_options.tempPrefix, "balsa_" );
    }

    if ( input_options.readType == PAIR_END_READ )
    {
        if ( input_options.insert_high == -1 )
        {
            fprintf ( stderr, "Please specify the maximum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low == -1 )
        {
            fprintf ( stderr, "Please specify the minimum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low > input_options.insert_high )
        {
            fprintf ( stderr, "The minimum value of insert size should not be greater than the maximum value of insert size.\n" );
            return false;
        }
    }

    // only allow three cases:
    // case 1: user does not specify "-D", "-A" or "-R"
    // case 2: user only specifies "-D" and "-A"
    // case 3: user only specifies "-D", "-A" and "-R"
    if ( readGroupOptionSupplied && ( ( !sampleNameSupplied ) || ( !readGroupSupplied ) ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, both options '-A' and '-R' have to be specified too.\n" );
        return false;
    }

    if ( sampleNameSupplied && ( !readGroupSupplied ) )
    {
        fprintf ( stderr, "Since an option '-A' is specified, the options '-D' has to be specified too.\n" );
        return false;
    }

    if ( readGroupSupplied && ( !sampleNameSupplied ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, the options '-A' has to be specified too.\n" );
        return false;
    }

    return true;
}



// print out the parameters
void printParameters ( InputOptions input_options, IniParams ini_params )
{
    fprintf (stderr, "\n----------Parameters------------------\n" );
    fprintf (stderr, "Number of mismatches: %i\n", input_options.numMismatch );

    if ( input_options.alignmentType == OUTPUT_ALL_VALID )
    { fprintf (stderr, "Output all valid alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_ALL_BEST )
    { fprintf (stderr, "Ouput all best alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_UNIQUE_BEST )
    { fprintf (stderr, "Ouput unique best alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_RANDOM_BEST )
    { fprintf (stderr, "Ouput random best alignments\n" ); }

    fprintf (stderr, "Number of CPU threads: %u\n", ini_params.Ini_NumOfCpuThreads );

    if ( input_options.insert_high > 0 )
    { fprintf (stderr, "Upper bound of insert size: %u\n", input_options.insert_high ); }

    if ( input_options.insert_low > 0 )
    { fprintf (stderr, "Lower bound of insert size: %u\n", input_options.insert_low ); }

    fprintf (stderr, "enableDP: %i\n", ( int ) input_options.enableDP );
    fprintf (stderr, "outputPrefix: %s\n", input_options.outputPrefix );
    fprintf (stderr, "readGroup: %s\n", input_options.readGroup );
    fprintf (stderr, "readGrpOption: %s\n", input_options.readGrpOption );
    fprintf (stderr, "--------------------------------------\n\n" );
    // seeding properties
    fprintf(stderr, "Number Of Round Of Deep DP For short reads: %d\n",ini_params.Ini_numberOfRoundOfDeepDpForShortReads);
    if (ini_params.Ini_numberOfRoundOfDeepDpForShortReads > 0){
        for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForShortReads;++i){
            fprintf(stderr, "Seed For Deep DP Round %d\n", i+1);
            fprintf(stderr, "DP Threshold: %d\n", ini_params.seedingPropertiesForShortReads[i].dpScoreThreshold);
            fprintf(stderr, "Number Of Seed For Deep DP Round %d\n",ini_params.seedingPropertiesForShortReads[i].numberOfSeed);
         
            for (int j=0;j<ini_params.seedingPropertiesForShortReads[i].numberOfSeed;++j){
                fprintf(stderr, "Seeding Round %d: Seed Position %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForShortReads[i].seedPos[j]);
                fprintf(stderr, "Seeding Round %d: Seed MinLength %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForShortReads[i].minSeedLength[j]);
                fprintf(stderr, "Seeding Round %d: Seed MaxLength %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForShortReads[i].maxSeedLength[j]);
            }
        }
        fprintf(stderr, "\n");
    }
    // seeding properties
    fprintf(stderr, "Number Of Round Of Deep DP For Long reads: %d\n",ini_params.Ini_numberOfRoundOfDeepDpForLongReads);
    if (ini_params.Ini_numberOfRoundOfDeepDpForLongReads > 0){
        for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForLongReads;++i){
            fprintf(stderr, "Seed For Deep DP Round %d\n", i+1);
            fprintf(stderr, "DP Threshold: %d\n", ini_params.seedingPropertiesForLongReads[i].dpScoreThreshold);
            fprintf(stderr, "Number Of Seed For Deep DP Round %d\n",ini_params.seedingPropertiesForLongReads[i].numberOfSeed);
         
            for (int j=0;j<ini_params.seedingPropertiesForLongReads[i].numberOfSeed;++j){
                fprintf(stderr, "Seeding Round %d: Seed Position %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForLongReads[i].seedPos[j]);
                fprintf(stderr, "Seeding Round %d: Seed MinLength %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForLongReads[i].minSeedLength[j]);
                fprintf(stderr, "Seeding Round %d: Seed MaxLength %d: %d\n",i+1,j+1,ini_params.seedingPropertiesForLongReads[i].maxSeedLength[j]);
            }
        }
        fprintf(stderr, "\n");
    }
}

// print out the parameters
void printDPParameters ( DPParameters dp_params )
{
#define BGS_OUTPUT_PARAMETERS
#ifdef BGS_OUTPUT_PARAMETERS
    fprintf (stderr, "\n----------DP Parameters------------------\n" );
    fprintf (stderr, "Match score: %i\n", dp_params.matchScore );
    fprintf (stderr, "Mismatch score: %i\n", dp_params.mismatchScore );
    fprintf (stderr, "Open gap score: %i\n", dp_params.openGapScore );
    fprintf (stderr, "Extend gap score: %i\n", dp_params.extendGapScore );
    fprintf (stderr, "numOfCPUThreads: %i\n", dp_params.numOfCPUThreads );
    fprintf (stderr, "numOfCPUForSeeding: %i\n", dp_params.numOfCPUForSeeding );
    fprintf (stderr, "softClipLeft: %i\n", dp_params.softClipLeft );
    fprintf (stderr, "softClipRight: %i\n", dp_params.softClipRight );
    fprintf (stderr, "tailTrimLen: %i\n", dp_params.tailTrimLen );
    fprintf (stderr, "singleDPSeedNum: %i\n", dp_params.singleDPSeedNum );

    for ( int i = 0; i < 3; i++ )
    {
        fprintf (stderr, "singleDPSeedPos[%i]: %i\n", i, dp_params.singleDPSeedPos[i] );
    }

    for ( int i = 0; i < 2; i++ )
    {
        fprintf (stderr, "maxHitNum[%i]: %i\n", i, dp_params.paramRead[i].maxHitNum );
//        fprintf (stderr, "sampleDist[%i]: %i\n", i, dp_params.paramRead[i].sampleDist );
//        fprintf (stderr, "seedLength[%i]: %i\n", i, dp_params.paramRead[i].seedLength );
    }
    
    fprintf (stderr, "--------------------------------------\n\n" );
#endif
}


// This function is to update the values inside InputOptions
// according to the i-th set of values inside the array of MultiInputItem
void updateInputOption ( InputOptions * input_options, MultiInputItem * multiInputArray, int i )
{
    input_options->queryFileName = multiInputArray[i].queryFile1;
    input_options->queryFileName2 = multiInputArray[i].queryFile2;
    input_options->insert_low = multiInputArray[i].insert_low;
    input_options->insert_high = multiInputArray[i].insert_high;
    input_options->outputPrefix = multiInputArray[i].outputPrefix;
    input_options->readGroup = multiInputArray[i].readGrpID;
    input_options->sampleName = multiInputArray[i].sampleName;
    input_options->readGrpOption = multiInputArray[i].readGrpOpt;
}

void freeIniParamsSeedingProperties(IniParams ini_params)
{
    for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForLongReads;++i)
    {
        free ( ini_params.seedingPropertiesForLongReads[i].maxSeedLength );
        free ( ini_params.seedingPropertiesForLongReads[i].seedPos );
        free ( ini_params.seedingPropertiesForLongReads[i].minSeedLength );
    }
    free ( ini_params.seedingPropertiesForLongReads );
    
    for (int i=0;i<ini_params.Ini_numberOfRoundOfDeepDpForShortReads;++i){
        free ( ini_params.seedingPropertiesForShortReads[i].maxSeedLength );
        free ( ini_params.seedingPropertiesForShortReads[i].seedPos );
        free ( ini_params.seedingPropertiesForShortReads[i].minSeedLength );
    }
    free ( ini_params.seedingPropertiesForShortReads );
}

void freeIniParamsFilenames ( IniParams ini_params )
{

}

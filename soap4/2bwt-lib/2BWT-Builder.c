/*

   BWTFormatdb.c        Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

   Copyright (C) 2006, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Thus, changing all references to 2bwt lib to subdirectory.

   Date   : 23rd October 2011
   Author : Edward MK Wu
   Change : Fix a rounding error when building reverse packed sequence.
   
*/

#include <stdio.h>
#include <stdlib.h>

#include "TypeNLimit.h"
#include "BWTConstruct.h"
#include "LTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "iniparser.h"
#include "HSP.h"
#include "Timing.h"

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);

// Parameters
char IniFileName[MAX_FILENAME_LEN+1];
int Confirmation;

// BuildTasks parameters
int ParseFASTA = TRUE;
int BuildBWT = TRUE;
int BuildSaValue = TRUE;
int BuildLookup = TRUE;
int BuildReverse=TRUE;
// Memory parameters
unsigned long long PoolSize = 2097152;                // 2M  - fixed; not configurable through ini

// Display parameters
int ShowProgress = TRUE;

// Database parameters
char FASTAFileName[MAX_FILENAME_LEN+1] = "";
char DatabaseName[MAX_FILENAME_LEN+1] = "";
char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.index.ann";
char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.index.amb";
char TranslateFileName[MAX_FILENAME_LEN+1] = "*.index.tra";
char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.pac";
char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.bwt";
char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.fmv";
char SaValueFileName[MAX_FILENAME_LEN+1] = "*.index.sa";
char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.index.sai";

char RevPackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.rev.pac";
char RevBWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.rev.bwt";
char RevBWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.rev.fmv";

char LookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.lkt";
char RevLookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.rev.lkt";


// Parse FASTA parameters
unsigned long long FASTARandomSeed = 0;
int MaskLowerCase = FALSE;

// Build BWT parameters
unsigned int OccValueFreq = 256;
float TargetNBit = 8;
unsigned long long InitialMaxBuildSize = 10000000;
unsigned long long IncMaxBuildSize = 10000000;

// Build SA value parameters
unsigned long long SaValueFreq = 16;

// Build Lookup Table parameters
unsigned int LookUpTableSize = 13;
    
void BuildReversePacked(const char *inputFileName, unsigned long long *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord) {

    FILE *inputFile;
    FILE *outputFile;
    unsigned char * packedText;
    unsigned char * revPackedText;
    off64_t packedFileLen;
    unsigned char lastByteLength;
    long long i,j;
    int k,l;

    inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");
    outputFile = (FILE*)(FILE*)fopen64(RevPackedDNAFileName, "wb");

    if (inputFile == NULL) {
        fprintf(stderr, "BuildReversePacked() : Cannot open inputFileName!\n");
        exit(1);
    }

    fseek(inputFile, -1, SEEK_END);
    packedFileLen = ftello64(inputFile);
    if (packedFileLen == -1) {
        fprintf(stderr, "BuildReversePacked(): Cannot determine file length!\n");
        exit(1);
    }
    
    
    fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);
    *textLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

    if (ShowProgress) {
        printf("Packed file size = %llu\n",(unsigned long long) packedFileLen);
        printf("Text Length = %llu\n",*textLength);
    }
    
    unsigned long long byteToProcess = (*textLength+CHAR_PER_BYTE-1) / CHAR_PER_BYTE;
    packedText = (unsigned char*) malloc(byteToProcess+1);
    revPackedText = (unsigned char*) malloc(byteToProcess+1);

    fseek(inputFile, 0, SEEK_SET);
    fread(packedText, 1, packedFileLen, inputFile);
    fclose(inputFile);
    i=byteToProcess-1;
    j=0;
    revPackedText[j]=0;
    k=0;
    
    unsigned char allOneChar = (1<<BIT_PER_CHAR) - 1;
    if (lastByteLength>0) {
        unsigned char lastByte = packedText[i];
        lastByte >>= (CHAR_PER_BYTE - lastByteLength)*BIT_PER_CHAR;
        for (k=0;k<lastByteLength;k++) {
            revPackedText[j]<<=BIT_PER_CHAR;
            revPackedText[j]|=(lastByte & allOneChar);
            lastByte>>=BIT_PER_CHAR;
        }
        i--;
    }

    for (;i>=0;i--) {
        unsigned char lastByte = packedText[i];
        for (l=0;l<CHAR_PER_BYTE;l++) {
            revPackedText[j]<<=BIT_PER_CHAR;
            revPackedText[j]|=(lastByte & allOneChar);
            k++;
            if (k>=CHAR_PER_BYTE) {
                j++;
                k=0;
                revPackedText[j]=0;
            }
            lastByte>>=BIT_PER_CHAR;
        }
    }
    
    if (k!=0) {
        revPackedText[j]<<=(CHAR_PER_BYTE - k)*BIT_PER_CHAR;
    }

    fwrite(revPackedText,sizeof(unsigned char),byteToProcess,outputFile);
    if (lastByteLength==0) {
        fwrite(&lastByteLength,sizeof(unsigned char),1,outputFile);
    }
    fwrite(&lastByteLength,sizeof(unsigned char),1,outputFile);
    
    free(revPackedText);
    free(packedText);
    fclose(outputFile);
    //*/
}

int main(int argc, char** argv) {
    

    char c;
    MMPool *mmPool;
    dictionary *programInput;
    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0;

    char filename[MAX_FILENAME_LEN+1];
    BWT *bwt = NULL;
    BWT *rev_bwt = NULL;
    HSP *hsp = NULL;
    unsigned long long textLength = 0;
    unsigned long long numSeq;

    BWTInc *bwtInc = NULL;
    BWTInc *rev_bwtInc = NULL;

    // Program input
    programInput = ParseInput(argc, argv);
    // Ini
    if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
        *(argv[0] + strlen(argv[0]) - 4) = '\0';
    }
    sprintf(filename, "%s.ini", argv[0]);
    ParseIniFile(filename);
    ProcessIni();
    ValidateIni();
    
    if (Confirmation == TRUE) {
        printf("Press Y to go or N to cancel. ");
        c = (char)getchar();
        while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
            c = (char)getchar();
        }
        if (c == 'n' || c == 'N') {
            exit(0);
        }
    }

    startTime = setStartTime();

    MMMasterInitialize(1, 0, FALSE, NULL);
    mmPool = MMPoolCreate(PoolSize);

    // Parse FASTA file to produce packed DNA and annotation file
    if (ParseFASTA == TRUE) {

        printf("Parsing FASTA file..\n");
        numSeq = HSPParseFASTAToPacked(FASTAFileName, AnnotationFileName, PackedDNAFileName, AmbiguityFileName, TranslateFileName, FASTARandomSeed, MaskLowerCase);
        
        printf("Finished. Parsed %llu sequences.\n", numSeq);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");
        
        //Parse packed DNA to construct the packed reversed DNA
        if (BuildReverse == TRUE) {
            printf("Parsing FASTA file reverse..\n");
            unsigned long long textLen;
            BuildReversePacked(PackedDNAFileName,&textLen,TRUE,1);
            printf("Reversed Packed DNA generated..\n");
        }

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

    }
    
    if (BuildLookup == TRUE) {
        printf("Building Look-Up..\n");
        BuildLookupTable(PackedDNAFileName,LookupTableFileName,LookUpTableSize);
        if (BuildReverse == TRUE) {
            BuildLookupTable(RevPackedDNAFileName,RevLookupTableFileName,LookUpTableSize);
        }
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Finished.\nElapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");
    }
    
    // Construct BWTInc from text
    if (BuildBWT == TRUE) {

        printf("Building BWT..\n");

        bwtInc = BWTIncConstructFromPacked(mmPool, PackedDNAFileName, ShowProgress, 
                                           TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

        printf("Finished constructing BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        printf("Saving BWT..\n");
        BWTSaveBwtCodeAndOcc(bwtInc->bwt, BWTCodeFileName, BWTOccValueFileName);
        printf("Finished saving BWT.  ");
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = bwtInc->bwt->textLength;
        BWTIncFree(mmPool, bwtInc);


        if (BuildReverse == TRUE) {
        //Building Reversed BWT
            printf("Building Reversed BWT..\n");

            rev_bwtInc = BWTIncConstructFromPacked(mmPool, RevPackedDNAFileName, ShowProgress, 
                                               TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

            printf("Finished constructing Reversed BWT in %u iterations.  ", rev_bwtInc->numberOfIterationDone);
            
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");

            printf("Saving BWT..\n");
            BWTSaveBwtCodeAndOcc(rev_bwtInc->bwt, RevBWTCodeFileName, RevBWTOccValueFileName);
            printf("Finished saving BWT.  ");
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");

            textLength = rev_bwtInc->bwt->textLength;
            BWTIncFree(mmPool, rev_bwtInc);
        }
    }

    // Load BWT
    if (BuildSaValue) {

        printf("Loading BWT...\n");

        bwt = BWTLoad(mmPool, 0, BWTCodeFileName, BWTOccValueFileName, NULL);
        //Use BWT to build the hash table

        printf("Finished loading BWT.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = bwt->textLength;

        printf("Building SA value...\n");
        
        if (ShowProgress) {
            BWTGenerateSaValue(bwt, SaValueFreq, bwt->textLength / SaValueFreq / 10);
        } else {
            BWTGenerateSaValue(bwt, SaValueFreq, 0);
        }
        BWTSaveSaValue(bwt, SaValueFileName);

        printf("Finished building SA value.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

    }

    // Free BWT
    if (BuildSaValue) {
        BWTFree(mmPool, bwt, false);
    }
    
    // Finished all construction tasks
    printf("Index building is completed.\n");
    totalElapsedTime = getElapsedTime(startTime);
    printf("Total elapsed time = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, totalElapsedTime);
    printf("\n");

    MMPoolFree(mmPool);
    iniparser_freedict(programInput);
    return 0;
}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;
    char t1[3] = "-c";    // specify that this is a boolean type parameter
    char t2[3] = "-U";    // specify that this is a boolean type parameter
    char *d[2];

    d[0] = t1;
    d[1] = t2;
    
    programInput = paraparser_load(argc, argv, 2, d);    // 2 boolean type parameters

    // Get database name
    if (!iniparser_find_entry(programInput, "argument:1")) {
        printf("Usage: ./%s <sequence file>\n",argv[0]);
        exit(1);
    }
    iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
    if (strlen(DatabaseName) + 4 > MAX_FILENAME_LEN) {
        printf("Usage: ./%s <sequence file>\n",argv[0]);
        exit(1);
    }

    // Get FASTA file name
    iniparser_copystring(programInput, "argument:2", FASTAFileName, DatabaseName, MAX_FILENAME_LEN);
    if (strlen(FASTAFileName) > MAX_FILENAME_LEN) {
        printf("Usage: ./%s <sequence file>\n",argv[0]);
        exit(1);
    }


    // Whether confirmation is needed
    Confirmation = iniparser_find_entry(programInput, "parameter:-c");

    MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");

    return programInput;

}

void ParseIniFile(char *iniFileName) {

    dictionary *ini;

    printf("Loading %s ..", iniFileName);
    ini = iniparser_load(iniFileName, FALSE);
    if (ini == NULL) {
        printf("not found.\n");
        return;
    }
    //printf("done.\n");

    // BuildTasks parameters
    ParseFASTA = iniparser_getboolean(ini, "BuildTasks:ParseFASTA", ParseFASTA);
    BuildBWT = iniparser_getboolean(ini, "BuildTasks:BuildBWT", BuildBWT);
    BuildSaValue = iniparser_getboolean(ini, "BuildTasks:BuildSaValue", BuildSaValue);
    BuildLookup = iniparser_getboolean(ini, "BuildTasks:BuildLookup", BuildLookup);
    BuildReverse = iniparser_getboolean(ini, "BuildTasks:BuildReverse", BuildReverse);

    // Display parameters
    ShowProgress = iniparser_getboolean(ini, "Display:ShowProgress", ShowProgress);

    // Parse FASTA parameters
    FASTARandomSeed = iniparser_getint(ini, "ParseFASTA:RandomSeed", FASTARandomSeed);
    if (FASTARandomSeed == 0) {
        FASTARandomSeed = getRandomSeed();
    }

    // Build BWT parameters
    OccValueFreq = iniparser_getint(ini, "BuildBWT:OccValueFreq", OccValueFreq);
    TargetNBit = (float)iniparser_getdouble(ini, "BuildBWT:TargetNBit", TargetNBit);
    InitialMaxBuildSize = iniparser_getint(ini, "BuildBWT:InitialMaxBuildSize", InitialMaxBuildSize);
    IncMaxBuildSize = iniparser_getint(ini, "BuildBWT:IncMaxBuildSize", IncMaxBuildSize);

    // Build SA value parameters
    SaValueFreq = iniparser_getint(ini, "BuildSAValue:SaValueFreq", SaValueFreq);

    // Build Look Up parameters
    LookUpTableSize = iniparser_getint(ini, "BuildLookUp:TableSize", LookUpTableSize);

    // Database parameters
    iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:TranslateFileName", TranslateFileName, TranslateFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);
    
    iniparser_copystring(ini, "Database:RevPackedDNAFileName", RevPackedDNAFileName, RevPackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevBWTCodeFileName", RevBWTCodeFileName, RevBWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevBWTOccValueFileName", RevBWTOccValueFileName, RevBWTOccValueFileName, MAX_FILENAME_LEN);

    iniparser_copystring(ini, "Database:LookupTableFileName", LookupTableFileName, LookupTableFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevLookupTableFileName", RevLookupTableFileName, RevLookupTableFileName, MAX_FILENAME_LEN);

    
    iniparser_freedict(ini);

}

void ProcessIni() {

    ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseName);
    ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseName);
    ProcessFileName(TranslateFileName, TranslateFileName, DatabaseName);
    ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseName);
    ProcessFileName(RevPackedDNAFileName, RevPackedDNAFileName, DatabaseName);
    ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseName);
    ProcessFileName(RevBWTCodeFileName, RevBWTCodeFileName, DatabaseName);
    ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseName);
    ProcessFileName(RevBWTOccValueFileName, RevBWTOccValueFileName, DatabaseName);
    ProcessFileName(SaValueFileName, SaValueFileName, DatabaseName);
    ProcessFileName(SaIndexFileName, SaIndexFileName, DatabaseName);
    ProcessFileName(LookupTableFileName, LookupTableFileName, DatabaseName);
    ProcessFileName(RevLookupTableFileName, RevLookupTableFileName, DatabaseName);
    
}

void ValidateIni() {

    if (!ParseFASTA && !BuildBWT && !BuildSaValue && !BuildLookup) {
        fprintf(stderr, "No action is specified!\n");
        exit(1);
    }
    if (ParseFASTA) {
        if (PackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Packed DNA file name is not specified!\n");
            exit(1);
        }
        if (AnnotationFileName[0] == '\0') {
            fprintf(stderr, "Annotation file name is not specified!\n");
            exit(1);
        }
        if (AmbiguityFileName[0] == '\0') {
            fprintf(stderr, "Ambiguity file name is not specified!\n");
            exit(1);
        }
    }
    if (BuildBWT) {
        if (PackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Packed DNA file is not specified!\n");
            exit(1);
        }
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file name is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file name is not specified!\n");
            exit(1);
        }
        if (TargetNBit < 2.5) {
            fprintf(stderr, "Target NBit should be at least 2.5!\n");
            exit(1);
        }
    }
    if (BuildSaValue) {
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file is not specified!\n");
            exit(1);
        }
        if (SaValueFileName[0] == '\0') {
            fprintf(stderr, "SA value file name is not specified!\n");
            exit(1);
        }
        if (SaValueFreq <= 0) {
            fprintf(stderr, "SA value frequency must > 0!\n");
            exit(1);
        }
    }
}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName) {

    char tempChar[MAX_FILENAME_LEN];
    unsigned long long i;

    if (inputFileName == NULL) {
        if (outputFileName != inputFileName) {
            outputFileName[0] = '\0';
        }
        return;
    }

    if (strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
        fprintf(stderr, "File length is too long!\n");
        exit(1);
    }

    strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

    // locate the *
    for (i=0; i<MAX_FILENAME_LEN; i++) {
        if (tempChar[i] == '*') {
            break;
        }
    }
    if (i<MAX_FILENAME_LEN) {
        tempChar[i] = '\0';
        sprintf(outputFileName, "%s%s%s", tempChar, databaseName, tempChar + i + 1);
    } else {
        sprintf(outputFileName, "%s", tempChar);
    }

}

/*

   HSP.c        BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

   Copyright (C) 2004, Wong Chi Kwong.

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

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <ctype.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "Socket.h"
#include "r250.h"
#include "HSP.h"
#include "HSPstatistic.h"


extern  double stat_expectationValue;

void HSPFillCharMap(unsigned char charMap[255]) {

    int i;

    for (i=0; i<255; i++) {
        charMap[i] = nonMatchDnaCharIndex;
    }
    for (i=0; i<DNA_CHAR_SIZE; i++) {
        charMap[dnaChar[i]] = (unsigned char)i;
        charMap[dnaChar[i] - 'A' + 'a'] = (unsigned char)i;
    }

}

HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName, const char *AmbiguityFileName, const char * TranslateFileName, const unsigned int trailerBufferInWord) {

    HSP *hsp;
    unsigned long long dnaLength;
    unsigned int randomSeed;
    int i;
    char c;
    unsigned char charMap[255];

    FILE *annotationFile = NULL, *ambiguityFile = NULL, *translateFile = NULL;

    hsp = (HSP*) MMPoolDispatch(mmPool, sizeof(HSP));

    // Load packed DNA
    if (PackedDNAFileName != NULL && PackedDNAFileName[0] != '\0' && PackedDNAFileName[0] != '-') {
        hsp->packedDNA = (unsigned int *) DNALoadPacked(PackedDNAFileName, &hsp->dnaLength, TRUE, trailerBufferInWord);
    } else {
        hsp->packedDNA = NULL;
        hsp->dnaLength = 0;
    }

    // Load annotation
    if (AnnotationFileName != NULL && AnnotationFileName[0] != '\0' && AnnotationFileName[0] != '-') {

        annotationFile = (FILE*)fopen64(AnnotationFileName, "r");
        if (annotationFile == NULL) {
            fprintf(stderr, "Cannot open annotation file!\n");
            exit(1);
        }

        fscanf(annotationFile, "%llu %d %u\n", &dnaLength, &hsp->numOfSeq, &randomSeed);
        if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
            fprintf(stderr, "Annotation database length not match!\n");
            exit(1);
        }
        hsp->dnaLength = dnaLength;
        if (hsp->numOfSeq == 0) {
            fprintf(stderr, "Annotation number of sequence = 0!\n");
            exit(1);
        }
        hsp->annotation = (Annotation*) MMUnitAllocate((hsp->numOfSeq+1) * sizeof(Annotation));
        hsp->seqOffset = (SeqOffset*) MMUnitAllocate((hsp->numOfSeq+1) * sizeof(SeqOffset));

        i = 0;
        hsp->minSeqLength = UINT32_MAX;
        while (!feof(annotationFile) && i < hsp->numOfSeq) {
            fscanf(annotationFile, "%u ", &hsp->annotation[i].gi);
            fgets(hsp->annotation[i].text, MAX_SEQ_NAME_LENGTH, annotationFile);

            if (strlen(hsp->annotation[i].text) != 0 && hsp->annotation[i].text[strlen(hsp->annotation[i].text)-1] == '\n') {
                hsp->annotation[i].text[strlen(hsp->annotation[i].text)-1] = '\0';
            }

            fscanf(annotationFile, "%llu %llu %d\n", &hsp->seqOffset[i].startPos, &hsp->seqOffset[i].endPos, &hsp->seqOffset[i+1].firstAmbiguityIndex);
            hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].firstAmbiguityIndex;
            if (hsp->seqOffset[i].endPos < hsp->minSeqLength) {
                hsp->minSeqLength = hsp->seqOffset[i].endPos;
            }
            hsp->seqOffset[i].endPos = hsp->seqOffset[i].startPos + hsp->seqOffset[i].endPos - 1;    // length of sequence is stored

            i++;
        }
        if (i < hsp->numOfSeq) {
            fprintf(stderr, "Annotation missing entries!\n");
            exit(1);
        }
        fclose(annotationFile);

        hsp->annotation[i].gi = 0;
        hsp->annotation[i].text[0] = '\0';

        hsp->seqOffset[i].startPos = UINT32_MAX;
        hsp->seqOffset[i].endPos = UINT32_MAX;
        hsp->seqOffset[0].firstAmbiguityIndex = 1;    // ambiguity[0] and ambiguity[numOfAmbiguity+1] are dummy
        for (i=1; i<=hsp->numOfSeq; i++) {
            hsp->seqOffset[i].firstAmbiguityIndex += hsp->seqOffset[i-1].firstAmbiguityIndex;    // number of ambiguity is stored
        }
        // hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex = total number of ambiguity + 1 now
        hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1;
        // hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = total number of ambiguity now
        for (i=hsp->numOfSeq; i>1; i--) {
            hsp->seqOffset[i-1].lastAmbiguityIndex = hsp->seqOffset[i].lastAmbiguityIndex - hsp->seqOffset[i-1].lastAmbiguityIndex;    // number of ambiguity is stored
        }
        for (i=0; i<hsp->numOfSeq; i++) {
            hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].lastAmbiguityIndex;
        }

    } else {

        // Make up a dummy annotation and a dummy sequence offset
        hsp->annotation = (Annotation*) MMUnitAllocate((1+1) * sizeof(Annotation));
        hsp->seqOffset = (SeqOffset*) MMUnitAllocate((1+1) * sizeof(SeqOffset));

        hsp->numOfSeq = 1;
        hsp->numOfAmbiguity = 0;
        
        hsp->annotation[0].gi = 0;
        hsp->annotation[0].text[0] = '\0';
        hsp->annotation[1].gi = 0;
        hsp->annotation[1].text[0] = '\0';
        hsp->seqOffset[0].startPos = 0;
        hsp->seqOffset[0].endPos = hsp->dnaLength - 1;
        hsp->seqOffset[0].firstAmbiguityIndex = 1;
        hsp->seqOffset[0].lastAmbiguityIndex = 0;
        hsp->seqOffset[1].startPos = UINT32_MAX;
        hsp->seqOffset[1].endPos = UINT32_MAX;
        hsp->seqOffset[1].firstAmbiguityIndex = UINT32_MAX;    // should not be referred
        hsp->seqOffset[1].lastAmbiguityIndex = UINT32_MAX;    // should not be referred

    }

    // Load ambigity
    if (AmbiguityFileName != NULL && AmbiguityFileName[0] != '\0' && AmbiguityFileName[0] != '-') {

        // Load ambigity
        ambiguityFile = (FILE*)fopen64(AmbiguityFileName, "r");
        if (ambiguityFile == NULL) {
            fprintf(stderr, "Cannot open ambiguity file!\n");
            exit(1);
        }

        fscanf(ambiguityFile, "%llu %d %d\n", &dnaLength, &i, &hsp->numOfAmbiguity);
        if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
            fprintf(stderr, "Ambiguity database length not match!\n");
            exit(1);
        }
        hsp->dnaLength = dnaLength;
        if (i != hsp->numOfSeq) {
            fprintf(stderr, "Ambiguity database number of sequence not match!\n");
            exit(1);
        }
        if (hsp->numOfAmbiguity != hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1) {
            fprintf(stderr, "Ambiguity number of ambiguity not match!\n");
            exit(1);
        }

        HSPFillCharMap(charMap);

        hsp->ambiguity = (Ambiguity*) MMUnitAllocate((hsp->numOfAmbiguity + 2) * sizeof(Ambiguity));
        hsp->ambiguity[0].startPos = 0;
        hsp->ambiguity[0].rightOfEndPos = 0;
        hsp->ambiguity[0].symbol = 0;
        i = 1;
        while (!feof(ambiguityFile) && i-1 < hsp->numOfAmbiguity) {
            fscanf(ambiguityFile, "%llu %llu %c\n", &hsp->ambiguity[i].startPos, &hsp->ambiguity[i].rightOfEndPos, &c);
            hsp->ambiguity[i].rightOfEndPos = hsp->ambiguity[i].startPos + hsp->ambiguity[i].rightOfEndPos;    // number of character is stored
            hsp->ambiguity[i].symbol = charMap[c];
            i++;
        }
        hsp->ambiguity[i].startPos = UINT32_MAX;
        hsp->ambiguity[i].rightOfEndPos = UINT32_MAX;
        hsp->ambiguity[i].symbol = 0;

        if (i-1 < hsp->numOfAmbiguity) {
            fprintf(stderr, "Ambiguity missing entries!\n");
            exit(1);
        }
        fclose(ambiguityFile);

    } else {

        if (hsp->numOfAmbiguity > 0) {
            fprintf(stderr, "Ambiguity file missing!\n");
            exit(1);
        }

        hsp->ambiguity = (Ambiguity*) MMUnitAllocate((0 + 2) * sizeof(Ambiguity));
        hsp->ambiguity[0].startPos = 0;
        hsp->ambiguity[0].rightOfEndPos = 0;
        hsp->ambiguity[0].symbol = 0;
        hsp->ambiguity[1].startPos = UINT32_MAX;
        hsp->ambiguity[1].rightOfEndPos = UINT32_MAX;
        hsp->ambiguity[1].symbol = 0;

    }
    
    // Load translate
    if (TranslateFileName != NULL && TranslateFileName[0] != '\0' && TranslateFileName[0] != '-') {

        // Load translate, ambiguity map
        translateFile = (FILE*)fopen64(TranslateFileName, "r");
        if (translateFile == NULL) {
            fprintf(stderr, "Cannot open translate file!\n");
            exit(1);
        }
        unsigned int gridEntries;
        unsigned int removedSegmentCount;
        fscanf(translateFile, "%llu %d %u %u\n", &dnaLength, &i, &removedSegmentCount, &gridEntries);
        if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
            fprintf(stderr, "Translate database length not match!\n");
            exit(1);
        }
        hsp->dnaLength = dnaLength;
        if (i != hsp->numOfSeq) {
            fprintf(stderr, "Translate database number of sequence not match!\n");
            exit(1);
        }
        /*if (hsp->numOfAmbiguity != hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1) {
            fprintf(stderr, "Translate number of ambiguity not match!\n");
            exit(1);
        }*/
        hsp->numOfRemovedSegment = removedSegmentCount;
        hsp->numOfGridEntry = gridEntries;
        hsp->seqActualOffset = (SeqActualOffset*) MMUnitAllocate((hsp->numOfSeq+1) * sizeof(SeqActualOffset));
        hsp->ambiguityMap = (unsigned int*) MMUnitAllocate((gridEntries) * sizeof(unsigned int));//MMUnitAllocate((gridEntries) * sizeof(unsigned int));
        hsp->translate = (Translate*) MMUnitAllocate((hsp->numOfSeq+removedSegmentCount) * sizeof(Translate));//MMUnitAllocate((hsp->numOfSeq+hsp->numOfAmbiguity) * sizeof(Translate));
        
        unsigned long long j=0;
        while (!feof(translateFile) && j<gridEntries) {
            unsigned int tmp;
            fscanf(translateFile, "%u\n", &tmp);
            hsp->ambiguityMap[j]=tmp;j++;
        }
        if (j < gridEntries) {
            fprintf(stderr, "Translate missing entries!\n");
            exit(1);
        }
        
        j=0;
        while (!feof(translateFile) && j<hsp->numOfSeq+removedSegmentCount) {
            fscanf(translateFile, "%llu %u %llu\n", &hsp->translate[j].startPos, &hsp->translate[j].chrID, &hsp->translate[j].correction);
            j++;
        }
        if (j < hsp->numOfSeq+removedSegmentCount) {
            fprintf(stderr, "Translate missing entries!\n");
            exit(1);
        }
        
        j=0;
        while (!feof(translateFile) && j<hsp->numOfSeq) {
            fscanf(translateFile, "%llu %llu\n", &hsp->seqActualOffset[j].startPos, &hsp->seqActualOffset[j].endPos);
            hsp->seqActualOffset[j].endPos = hsp->seqActualOffset[j].startPos + hsp->seqActualOffset[j].endPos;
            j++;
        }
        // Backward compatibility
        if (j < hsp->numOfSeq) {
            fprintf(stderr, "[Translation] Missing actual chromosome lengths.\n");
            fprintf(stderr, "Re-build Packed indexes for more accurate chromosome length in output.\n");
            
            //Translation file is missing, initialise the value with inaccurate chromosome length
            j=0;
            while (j<hsp->numOfSeq) {
                hsp->seqActualOffset[j].startPos = hsp->seqOffset[j].startPos;
                hsp->seqActualOffset[j].endPos  = hsp->seqOffset[j].endPos;
                j++;
            }
        }
        fclose(translateFile);

    } else {

        if (hsp->numOfAmbiguity > 0 || hsp->numOfSeq > 0) {
            fprintf(stderr, "Translate file missing!\n");
            exit(1);
        }

        hsp->numOfRemovedSegment = 0;
        unsigned int gridEntries = (hsp->dnaLength/GRID_SAMPLING_FACTOR)+1;
        hsp->numOfGridEntry = gridEntries;
        hsp->seqActualOffset = (SeqActualOffset*) MMUnitAllocate((hsp->numOfSeq+1) * sizeof(SeqActualOffset));
        hsp->ambiguityMap = (unsigned int*) MMUnitAllocate((gridEntries) * sizeof(unsigned int));
        hsp->translate = (Translate*) MMUnitAllocate((hsp->numOfSeq) * sizeof(Translate));//MMUnitAllocate((hsp->numOfSeq+hsp->numOfAmbiguity) * sizeof(Translate));
        unsigned int j=0;
        while (!feof(translateFile) && j<gridEntries) {
            hsp->ambiguityMap[j]=0;
        }
        hsp->translate[0].startPos=0;
        hsp->translate[0].chrID=0;
        hsp->translate[0].correction=0;
        
        j=0;
        while (j<hsp->numOfSeq) {
            hsp->seqActualOffset[j].startPos = hsp->seqOffset[j].startPos;
            hsp->seqActualOffset[j].endPos  = hsp->seqOffset[j].endPos;
            j++;
        }
    }
    return hsp;

}

void HSPFree(MMPool *mmPool, HSP *hsp, const unsigned int trailerBufferInWord) {

    if (hsp->packedDNA != NULL) {
        DNAFreePacked(hsp->packedDNA, hsp->dnaLength, trailerBufferInWord,0);
    }
    MMUnitFree(hsp->seqOffset, (hsp->numOfSeq+1) * sizeof(SeqOffset));
    MMUnitFree(hsp->annotation, (hsp->numOfSeq+1) * sizeof(Annotation));
    MMUnitFree(hsp->ambiguity, (hsp->numOfAmbiguity+2) * sizeof(Ambiguity));
    MMUnitFree(hsp->seqActualOffset, (hsp->numOfSeq+1) * sizeof(SeqActualOffset));
    MMUnitFree(hsp->ambiguityMap, (hsp->numOfGridEntry) * sizeof(unsigned int));
    MMUnitFree(hsp->translate, (hsp->numOfSeq+hsp->numOfRemovedSegment) * sizeof(Translate));

    MMPoolReturn(mmPool, hsp, sizeof(hsp));

}

unsigned int HSPParseFASTAToPacked(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName, const char * translateFileName,
                      const unsigned int FASTARandomSeed, const int maskLowerCase) {

    FILE *FASTAFile, *annotationFile, *packedDNAFile, *ambiguityFile, * translateFile;
    Annotation *annotation;
    SeqOffset *seqOffset;
    SeqActualOffset *seqActualOffset;
    Ambiguity *ambiguity;
    int annotationAllocated = 256;
    int ambiguityAllocated =  256;
    unsigned long long totalDiscardedChar = 0;
    
    char c;
    int i;
    int fastaComment = 0;
    int numAmbiguity;
    unsigned long long lastAmbiguityPos;
    unsigned long long numChar, numCharInBuffer, totalNumChar, totalActualNumChar;
    unsigned int sequenceRandomSeed;
    int numSeq;
    unsigned char charMap[255];

    unsigned char buffer[PACKED_BUFFER_SIZE];
    unsigned char packedBuffer[PACKED_BUFFER_SIZE / CHAR_PER_BYTE];

    FASTAFile = (FILE*)fopen64(FASTAFileName, "r");
    if (FASTAFile == NULL) {
        fprintf(stderr, "ParseFASTToPacked() : Cannot open FASTAFileName!\n");
        exit(1);
    }

    annotationFile = (FILE*)fopen64(annotationFileName, "w");
    if (annotationFile == NULL) {
        fprintf(stderr, "ParseFASTToPacked() : Cannot open annotationFileName!\n");
        exit(1);
    }

    packedDNAFile = (FILE*)fopen64(packedDNAFileName, "wb");
    if (packedDNAFile == NULL) {
        fprintf(stderr, "ParseFASTToPacked() : Cannot open packedDNAFileName!\n");
        exit(1);
    }

    ambiguityFile = (FILE*)fopen64(ambiguityFileName, "w");
    if (ambiguityFile == NULL) {
        fprintf(stderr, "ParseFASTToPacked() : Cannot open ambiguityFileName!\n");
        exit(1);
    }

    translateFile = (FILE*)fopen64(translateFileName, "w");
    if (translateFile == NULL) {
        fprintf(stderr, "ParseFASTToPacked() : Cannot open translateFileName!\n");
        exit(1);
    }

    HSPFillCharMap(charMap);

    c = (char)getc(FASTAFile);
    if (c != '>') {
        fprintf(stderr, "ParseFASTToPacked() : FASTA file does not begin with '>'!\n");
        exit(1);
    }

    totalNumChar = 0;
    totalActualNumChar = 0;
    numSeq = 0;
    numAmbiguity = -1;
    numCharInBuffer = 0;

    annotation = (Annotation*) MMUnitAllocate(sizeof(Annotation) * annotationAllocated);
    seqOffset = (SeqOffset*) MMUnitAllocate(sizeof(SeqOffset) * annotationAllocated);
    ambiguity = (Ambiguity*) MMUnitAllocate(sizeof(Ambiguity) * ambiguityAllocated);
    seqActualOffset = (SeqActualOffset*) MMUnitAllocate(sizeof(SeqActualOffset) * annotationAllocated);

    while (!feof(FASTAFile)) {
        fastaComment = 0;
        numChar = 0;
        if (numSeq >= annotationAllocated) {
            annotation = (Annotation*) MMUnitReallocate(annotation, sizeof(Annotation) * annotationAllocated * 2, sizeof(Annotation) * annotationAllocated);
            seqOffset = (SeqOffset*) MMUnitReallocate(seqOffset, sizeof(SeqOffset) * annotationAllocated * 2, sizeof(SeqOffset) * annotationAllocated);
            seqActualOffset = (SeqActualOffset*) MMUnitReallocate(seqActualOffset, sizeof(SeqActualOffset) * annotationAllocated * 2, sizeof(SeqActualOffset) * annotationAllocated);
            annotationAllocated *= 2;
        }

        annotation[numSeq].gi = 0;

        c = (char)getc(FASTAFile);
        while (!feof(FASTAFile) && c != '\n') {
            if (!fastaComment && isspace(c)) {
                fastaComment = 1;
            }

            if (!fastaComment && numChar < MAX_SEQ_NAME_LENGTH) {
                annotation[numSeq].text[numChar] = c;
                numChar++;
            }
            c = (char)getc(FASTAFile);
        }
        annotation[numSeq].text[numChar] = '\0';

        if (numChar > 3 && annotation[numSeq].text[0] == 'g'
                             && annotation[numSeq].text[1] == 'i'
                             && annotation[numSeq].text[2] == '|') {
            sscanf(annotation[numSeq].text + 3, "%u", &(annotation[numSeq].gi));
        }

        // Set random seed for the sequence
        sequenceRandomSeed = FASTARandomSeed;
        if (annotation[numSeq].gi > 0) {
            sequenceRandomSeed += annotation[numSeq].gi;
        } else {
            for (i=0; i<(int)numChar; i++) {
                c = annotation[numSeq].text[i];
                sequenceRandomSeed += c;
            }
        }
        r250_init(sequenceRandomSeed);

        seqOffset[numSeq].startPos = totalNumChar;
        seqActualOffset[numSeq].startPos = totalActualNumChar;
        numChar = 0;
        lastAmbiguityPos = (unsigned int)-2;

        c = (char)getc(FASTAFile);
        unsigned long long nCount=0;
        while (!feof(FASTAFile) && c != '>') {
            // Get sequence
            if (c != '\n' && c != '\t') {
                if (maskLowerCase && c >= 'a' && c <= 'z') {
                    c = dnaChar[lowercaseDnaCharIndex];
                }
                
                if (charMap[c] != nonMatchDnaCharIndex) {
                    if (ambiguityCount[charMap[c]] == 1) {
                        buffer[numCharInBuffer] = c;
                    } else {
    //Char != A C G T...can be anything else on Earth!!
                        c = (char)getc(FASTAFile);
                        unsigned long long nCount=1;
                        while (!feof(FASTAFile) && c != '>') {
                            if (charMap[c] != nonMatchDnaCharIndex) {
                                if (ambiguityCount[charMap[c]] != 1) {
                                    nCount++;
                                } else {
                                    break;
                                }
                            }
                            c = (char)getc(FASTAFile);
                        }
                        if (nCount<10) {
                            int k;
                            for (k=0;k<nCount;k++) {
                                buffer[numCharInBuffer] = 'G'; //Substitute rest of the 'N' with 'G'
                                numCharInBuffer++;
                                if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
                                    ConvertTextToBytePacked(buffer, packedBuffer, charMap, ALPHABET_SIZE, PACKED_BUFFER_SIZE);
                                    fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / CHAR_PER_BYTE, packedDNAFile);
                                    numCharInBuffer = 0;
                                }
                                numChar++;
                            }
                        } else {
                            totalDiscardedChar+=nCount;
                            
                            numAmbiguity++;
                            if (numAmbiguity >= ambiguityAllocated) {
                                ambiguity = (Ambiguity*) MMUnitReallocate(ambiguity, sizeof(Ambiguity) * ambiguityAllocated * 2, sizeof(Ambiguity) * ambiguityAllocated);
                                ambiguityAllocated *= 2;
                            }
                            ambiguity[numAmbiguity].startPos = totalNumChar + numChar;
                            ambiguity[numAmbiguity].rightOfEndPos = totalNumChar + numChar + totalDiscardedChar - 1;
                            ambiguity[numAmbiguity].symbol = 0;
                        }
                        continue;
                    }
                    numCharInBuffer++;
                    if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
                        ConvertTextToBytePacked(buffer, packedBuffer, charMap, ALPHABET_SIZE, PACKED_BUFFER_SIZE);
                        fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / CHAR_PER_BYTE, packedDNAFile);
                        numCharInBuffer = 0;
                    }
                    numChar++;
                }
            }
            c = (char)getc(FASTAFile);
        }
        seqOffset[numSeq].endPos = totalNumChar + numChar - 1;
        seqOffset[numSeq].firstAmbiguityIndex = numAmbiguity + 1;
        seqActualOffset[numSeq].endPos = totalNumChar + totalDiscardedChar + numChar - 1;
        totalNumChar += numChar;
        totalActualNumChar = totalNumChar + totalDiscardedChar;
        numSeq++;

    }
    fclose(FASTAFile);


    // Finalize packed DNA file

    numAmbiguity++;

    if (numCharInBuffer > 0) {
        ConvertTextToBytePacked(buffer, packedBuffer, charMap, ALPHABET_SIZE, numCharInBuffer);
        fwrite(packedBuffer, 1, (numCharInBuffer + CHAR_PER_BYTE-1) / CHAR_PER_BYTE, packedDNAFile);
        numCharInBuffer = 0;
    }
    if (totalNumChar % CHAR_PER_BYTE == 0) {
        c = 0;
        fwrite(&c, 1, 1, packedDNAFile);
    }
    c = (char)(totalNumChar % CHAR_PER_BYTE);
    fwrite(&c, 1, 1, packedDNAFile);

    fclose(packedDNAFile);
    
    //Compute Translation Map
    unsigned int gridEntries = (totalNumChar/GRID_SAMPLING_FACTOR)+1;
    unsigned int *ambiguityMap = (unsigned int*)malloc(sizeof(unsigned int)*gridEntries);
    Translate * translate = (Translate*) malloc(sizeof(Translate)*(numAmbiguity+numSeq));
    unsigned int j=0,k=0;
    for (j=0;j<gridEntries;j++) {
        ambiguityMap[j]=0;
    }
    j=0;i=0;
    while(j<numAmbiguity && k<numSeq) {
        if (ambiguity[j].startPos<seqOffset[k].startPos) {
            //Write Ambiguous Segment Flag
            unsigned long long index = ambiguity[j].startPos;
            ambiguityMap[index/GRID_SAMPLING_FACTOR]+=1;
            translate[i].startPos=ambiguity[j].startPos;
            translate[i].chrID=k;
            
            if (k-1>0) {
                unsigned long long correction = ambiguity[seqOffset[k-2].firstAmbiguityIndex-1].rightOfEndPos -
                              ambiguity[seqOffset[k-2].firstAmbiguityIndex-1].startPos + 1;
                translate[i].correction=(seqOffset[k-1].startPos)+ 
                                        correction-
                                        (ambiguity[j].rightOfEndPos - ambiguity[j].startPos +1) - 1;
            } else {
                translate[i].correction=-ambiguity[j].rightOfEndPos + ambiguity[j].startPos-2;
            }
            j++;
            i++;
        } else if (ambiguity[j].startPos>seqOffset[k].startPos) {
            //Write Artificial Chromosome Start Flag
            unsigned long long index = seqOffset[k].startPos;
            ambiguityMap[index/GRID_SAMPLING_FACTOR]+=1;
            translate[i].startPos=seqOffset[k].startPos;
            translate[i].chrID=k+1;
            translate[i].correction=seqOffset[k].startPos-1;
            k++;
            i++;
        } else {
            k++;
        }
    }
    while (j<numAmbiguity) {
        unsigned long long index = ambiguity[j].startPos;
        ambiguityMap[index/GRID_SAMPLING_FACTOR]+=1;
        translate[i].startPos=ambiguity[j].startPos;
        translate[i].chrID=k;            

        if (k-1>0) {
            unsigned long long correction = ambiguity[seqOffset[k-2].firstAmbiguityIndex-1].rightOfEndPos -
                          ambiguity[seqOffset[k-2].firstAmbiguityIndex-1].startPos + 1;
            translate[i].correction=(seqOffset[k-1].startPos)+ 
                                    correction-
                                    (ambiguity[j].rightOfEndPos - ambiguity[j].startPos +1) - 1;
        } else {
            translate[i].correction=-ambiguity[j].rightOfEndPos + ambiguity[j].startPos-2;
        }
        //printf("%u\tAMB%u\tSP:%u\tCHR:%u\tCORRECT:%u\n",i,j,translate[i].startPos,translate[i].chrID,translate[i].correction);
        j++;i++;
    }
    while (k<numSeq) {
        unsigned long long index = seqOffset[k].startPos;
        ambiguityMap[index/GRID_SAMPLING_FACTOR]+=1;
        translate[i].startPos=seqOffset[k].startPos;
        translate[i].chrID=k+1;
        translate[i].correction=seqOffset[k].startPos-1;
        //printf("%u\tCHR%u\tSP:%u\tCHR:%u\tCORRECT:%u\n",i,k,translate[i].startPos,translate[i].chrID,translate[i].correction);    
        k++;i++;
    }
    while (i<numAmbiguity+numSeq) {
        if (i>0) {
            translate[i].startPos=translate[i-1].startPos;
            translate[i].chrID=translate[i-1].chrID;
            translate[i].correction=translate[i-1].correction;
        } i++;
    }
    
    for (j=1;j<gridEntries;j++) {
        ambiguityMap[j]+=ambiguityMap[j-1];
    }
    for (j=0;j<gridEntries;j++) {
        ambiguityMap[j]--;
    }

    // Output annotation file
    fprintf(annotationFile, "%llu %u %u\n", totalNumChar, numSeq, FASTARandomSeed);
    for (i=0; i<numSeq; i++) {
        fprintf(annotationFile, "%u %s\n", annotation[i].gi, annotation[i].text);
        fprintf(annotationFile, "%llu %llu 0\n", seqOffset[i].startPos, seqOffset[i].endPos - seqOffset[i].startPos + 1);
        // output number of ambiguity
        /*if (i > 0) {
            fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex - seqOffset[i-1].firstAmbiguityIndex);
        } else {
            fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex);
        }*/
    }
    fclose(annotationFile);

    MMUnitFree(annotation, sizeof(Annotation) * annotationAllocated);
    MMUnitFree(seqOffset, sizeof(SeqOffset) * annotationAllocated);

    // Output ambiguity file
    fprintf(ambiguityFile, "%llu %u %u\n", totalNumChar, numSeq, 0);//numAmbiguity);
    //for (i=0; i<numAmbiguity; i++) {
        // The ambiguity for the dummy length is visible
        //fprintf(ambiguityFile, "%u %u %c\n", ambiguity[i].startPos, ambiguity[i].rightOfEndPos - ambiguity[i].startPos + 1, dnaChar[ambiguity[i].symbol]);
    //}
    fclose(ambiguityFile);

    MMUnitFree(ambiguity, ambiguityAllocated * sizeof(Ambiguity));
    
    fprintf(translateFile, "%llu %u %u %u\n", totalNumChar, numSeq, numAmbiguity, gridEntries);
    for (j=0; j<gridEntries; j++) {
        // The ambiguity for the dummy length is visible
        fprintf(translateFile, "%u\n", ambiguityMap[j]);
    }
    for (j=0; j<numAmbiguity+numSeq; j++) {
        // The ambiguity for the dummy length is visible
        fprintf(translateFile, "%llu %u %llu\n", translate[j].startPos,translate[j].chrID, translate[j].correction );
    }
    for (j=0; j<numSeq; j++) {
        fprintf(translateFile, "%llu %llu\n", seqActualOffset[j].startPos, seqActualOffset[j].endPos - seqActualOffset[j].startPos + 1);
    }
    fclose(translateFile);
    
    free(translate);
    free(ambiguityMap);
    MMUnitFree(seqActualOffset, sizeof(SeqActualOffset) * annotationAllocated);

    return numSeq;

}

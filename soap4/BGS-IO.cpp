/*
 *
 *    BGS-IO.c
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

#include "BGS-IO.h"
#include "assert.h"
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

// mapping score
// dimension 1: # of mismatches (i.e. 0, 1, 2, 3, 4) OR
//              difference % from the full DP score (i.e. 0, 1-5%, 6-10%, 11-15%, 15-20%, >20%)
// dimension 2: average mismatch quality value (i.e. 1-10, 11-20, 21-30, 31-40)
//static int mapping_score[6][4] = {{40,40,40,40},{38,37,35,34},{35,33,30,28},{30,28,25,22},{25,22,19,16},{19,16,13,10}};
// static int mapping_score_orig[6][4] = {{MAPQ_MAX,MAPQ_MAX,MAPQ_MAX,MAPQ_MAX},{MAPQ_MAX*0.875,MAPQ_MAX*0.875,MAPQ_MAX*0.85,MAPQ_MAX*0.85},{MAPQ_MAX*0.75,MAPQ_MAX*0.75,MAPQ_MAX*0.7,MAPQ_MAX*0.7},{MAPQ_MAX*0.625,MAPQ_MAX*0.625,MAPQ_MAX*0.55,MAPQ_MAX*0.55},{MAPQ_MAX*0.475,MAPQ_MAX*0.475,MAPQ_MAX*0.4,MAPQ_MAX*0.4},{MAPQ_MAX*0.325,MAPQ_MAX*0.325,MAPQ_MAX*0.25,MAPQ_MAX*0.25}};
static double mapping_score[6][2] = {{1.0, 1.0}, {0.875, 0.85}, {0.75, 0.7}, {0.625, 0.55}, {0.475, 0.4}, {0.325, 0.25}};


//===========================================================//
// The following is for the updated MAPQ scoring function    //
// for single-end reads for DP module (Date: Oct 19, 2012)   //
//===========================================================//

// The penality score when consideration of the average quality value in mismatch positions
static float penalty_score_avg_mis_qual[41] = {3, 2.85, 2.71, 2.57, 2.43, 2.3, 2.17, 2.04, 1.92, 1.8, 1.69, 1.58, 1.47, 1.37, 1.27, 1.17, 1.08, 0.99, 0.91, 0.83, 0.75, 0.68, 0.61, 0.54, 0.48, 0.42, 0.37, 0.32, 0.27, 0.23, 0.19, 0.15, 0.12, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0, 0};

// The penality ratio when consideration of X1
static float penalty_ratio_x1[101] = {1, 0.5, 0.33, 0.25, 0.2, 0.17, 0.14, 0.13, 0.11, 0.1, 0.09, 0.08, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

//


//Define the below parameter to output the alignment result(text position)
// on screen instead of writing into output file
//#define DEBUG_2BWT_OUTPUT_TO_SCREEN

//Define the below parameter to output the alignment result(text position)
// with the old 2BWT-Aligner format
//#define DEBUG_2BWT_OUTPUT_32_BIT

//Define the below parameter to stop cache reported SA range and text position.
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_OUTPUT

//Define the below parameter to skip writing the alignment result(text position)
// into disk. This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_WRITE_FILE

//Define the below parameter to skip translating the alignment result(text position).
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_TRANSLATION

OCC * OCCConstruct ()
{
    OCC * occ = ( OCC * ) malloc ( sizeof ( OCC ) );
    SAMOccurrenceConstruct ( occ );
    occ->occPositionCacheCount = 0;
    return occ;
}

void OCCReset ( OCC * occ )
{
    occ->occPositionCacheCount = 0;
}

void OCCFree ( OCC * occ )
{
    SAMOccurrenceDestruct ( occ );
    free ( occ );
}

void OCCWriteOutputHeader ( HSP * hsp, FILE * outFilePtr,
                                    unsigned int maxReadLength,
                                    unsigned int numOfReads,
                                    int outputFormat )
{
    unsigned long long tp;

    switch ( outputFormat )
    {
        case SRA_OUTPUT_FORMAT_SAM_API:
            //Ad-hoc writing header is not supported by API
            break;

        case SRA_OUTPUT_FORMAT_SAM:
            fprintf ( outFilePtr, "@HD\tVN:1.4\tSO:unsorted\n" );

            for ( int i = 0; i < hsp->numOfSeq; i++ )
            {
                int j;
                for ( j = 0; j < 255; j++ )
                {
                    if ( hsp->annotation[i].text[j] == '\0' ||
                            hsp->annotation[i].text[j] == ' ' ||
                            hsp->annotation[i].text[j] == '\t' ||
                            hsp->annotation[i].text[j] == '\r' ||
                            hsp->annotation[i].text[j] == '\n' )
                    {
                        break;
                    }
                }

                hsp->annotation[i].text[j] = '\0';
                tp = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos + maxReadLength;
                fprintf ( outFilePtr, "@SQ\tSN:%s\tLN:%llu\n", hsp->annotation[i].text, tp );
            }

            fprintf ( outFilePtr, "@PG\tID:%s\tVN:v%d.%d.%d (%s)\n", PROJECT_NAME, PROJECT_MAJOR, PROJECT_MINOR, PROJECT_REV, PROJECT_SPECIAL );
            break;

        case SRA_OUTPUT_FORMAT_PLAIN:
            assert(false);
    }
}

void AssignCigarStrToSAMIU ( bam1_t * samAlgnmt, int * curSize,
                             char * cigar_string )
{
    char operCode[256];
    operCode['M'] = 0; // match or mismatch
    operCode['I'] = 1; // insertion to the reference
    operCode['D'] = 2; // deletion from the reference
    operCode['S'] = 4; // soft clipping
    samAlgnmt->core.n_cigar = 0;
    int num = 0;

    for ( int i = 0; i < strlen ( cigar_string ); i++ )
    {
        if ( cigar_string[i] >= '0' && cigar_string[i] <= '9' )
        {
            num = num * 10 + cigar_string[i] - '0';
        }
        else if ( num > 0 )
        {
            SAMIUint8ConcatUint32 ( samAlgnmt->data, curSize, num << BAM_CIGAR_SHIFT | operCode[cigar_string[i]] ); //CIGAR
            samAlgnmt->core.n_cigar++;                                //SAM: number of CIGAR operations
            num = 0;
        }
    }
}

unsigned long long getChrAndPos ( SRAQueryInput * qInput, unsigned long long ambPos,
                            unsigned long long * tp, unsigned int * chr_id )
{
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // get the chromosome and the position for the position on the packed sequence
    HSP * hsp = aIndex->hsp;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    Translate * occTranslate = hsp->translate;
    unsigned long long correctPosition;
    unsigned long long approxIndex, approxValue;
    correctPosition = ambPos;
    approxIndex = ambPos >> GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = occAmbiguityMap[approxIndex];

    while ( occTranslate[approxValue].startPos > ambPos )
    {
        approxValue--;
    }

    correctPosition -= occTranslate[approxValue].correction;
    ( *tp ) = correctPosition;
    ( *chr_id ) = occTranslate[approxValue].chrID;

    return approxValue < hsp->numOfRemovedSegment - 1 ?
      occTranslate[approxValue + 1].startPos - 1 : hsp->dnaLength;
}


long long BoundaryCheck ( unsigned long long pacPos, int chrID, unsigned long long chrEndPos, long long readLength, unsigned long long segmentEndPos,
                    HSP * hsp, unsigned long long & correctedPac, char ** newCIGAR )
{
    segmentEndPos = chrEndPos < segmentEndPos ? chrEndPos : segmentEndPos;

    if ( pacPos + readLength <= segmentEndPos + 1 ) { return 0; }

    long long actualAlignedLength = segmentEndPos - pacPos + 1;

    if ( newCIGAR ) { *newCIGAR = ( char * ) malloc ( 15 * sizeof ( char ) ); }

    if ( actualAlignedLength >= ( readLength + 1 ) / 2 )
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dM%dS", actualAlignedLength, readLength - actualAlignedLength ); }

        correctedPac = pacPos;
        return - ( readLength - actualAlignedLength ); // trim right
    }
    else
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dS%dM", actualAlignedLength, readLength - actualAlignedLength ); }

        correctedPac = segmentEndPos + 1;
        return actualAlignedLength; // trim left
    }
}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
long long BoundaryCheckDP ( unsigned long long pacPos, int chrID, unsigned long long chrEndPos, long long readLength, char * cigar, unsigned long long segmentEndPos, HSP * hsp,
                      unsigned long long & correctedPac, char ** newCIGAR )
{
    if ( pacPos + readLength * 2 <= chrEndPos + 1 && pacPos + readLength * 2 <= segmentEndPos + 1 ) { return 0; } // is this readLength * 2 bound good enough?

    segmentEndPos = chrEndPos < segmentEndPos ? chrEndPos : segmentEndPos;

    // CIGAR string parsing...
    char * buffer = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    char * leftBuf = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    char * rightBuf = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    long long leftLen, rightLen;
    long long leftOff = 0, rightOff = 0;
    long long leftS, rightS;
    char * cigarP = cigar;
    {
        unsigned long long refPos;
        refPos = pacPos;
        leftLen = rightLen = 0;
        leftBuf[0] = rightBuf[0] = 0;
        leftS = rightS = 0;

        while ( *cigarP )
        {
            long long num = 0;
            char op;

            while ( *cigarP <= '9' )
            {
                num = num * 10 + ( *cigarP - '0' );
                cigarP++;
            }

            op = *cigarP;
            cigarP++;

            if ( op == 'S' )
            {
                // refPos += num;
                sprintf ( buffer, "%dS", num );

                if ( refPos <= segmentEndPos )
                {
                    strcat ( leftBuf, buffer );
                    leftS += num;
                }
                else
                {
                    strcat ( rightBuf, buffer );
                    rightS += num;
                }
            }
            else if ( op == 'M' || op == 'm' )
            {
                if ( refPos > segmentEndPos )
                {
                    rightLen += num;
                    sprintf ( buffer, "%d%c", num, op );
                    strcat ( rightBuf, buffer );
                }
                else if ( refPos + num <= segmentEndPos + 1 )
                {
                    leftLen += num;
                    sprintf ( buffer, "%d%c", num, op );
                    strcat ( leftBuf, buffer );
                }
                else     // crossing boundary
                {
                    leftLen += segmentEndPos - refPos + 1;
                    sprintf ( buffer, "%d%c", segmentEndPos - refPos + 1, op );
                    strcat ( leftBuf, buffer );
                    rightLen += num - segmentEndPos + refPos - 1;
                    sprintf ( buffer, "%d%c", num - segmentEndPos + refPos - 1, op );
                    strcat ( rightBuf, buffer );
                }

                refPos += num;
            }
            else if ( op == 'D' )
            {
                sprintf ( buffer, "%dD", num );

                if ( refPos > segmentEndPos )
                {
                    if ( strlen ( rightBuf ) == 0 )
                    {
                        rightOff = num;
                    }
                    else
                    {
                        strcat ( rightBuf, buffer );
                    }
                }
                else if ( refPos + num <= segmentEndPos + 1 )
                {
                    strcat ( leftBuf, buffer );
                }
                else     // crossing boundary
                {
                    sprintf ( buffer, "%dD", segmentEndPos - refPos + 1 );
                    strcat ( leftBuf, buffer );
                    sprintf ( buffer, "%dD", num - segmentEndPos + refPos - 1 );

                    if ( strlen ( rightBuf ) == 0 )
                    {
                        rightOff = num - segmentEndPos + refPos - 1;
                    }
                    else
                    {
                        strcat ( rightBuf, buffer );
                    }
                }

                refPos += num;
            }
            else if ( op == 'I' )
            {
                if ( refPos == chrEndPos + 1 ) // exactly at boundary
                {
                    // ignore both sides
                }
                else if ( refPos <= segmentEndPos )
                {
                    leftLen += num;
                    sprintf ( buffer, "%dI", num );
                    strcat ( leftBuf, buffer );
                }
                else
                {
                    rightLen += num;
                    sprintf ( buffer, "%dI", num );
                    strcat ( rightBuf, buffer );
                }
            }
        }
    }

    if ( newCIGAR ) { *newCIGAR = buffer; }

    if ( !leftLen || !rightLen )
    {
        free ( leftBuf );
        free ( rightBuf );

        if ( newCIGAR )
        {
            free ( *newCIGAR );
            *newCIGAR = NULL;
        }

        return 0; // no trimming
    }
    else if ( leftLen >= rightLen )
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%s%dS", leftBuf, readLength - leftLen - leftS ); }

        correctedPac = pacPos;
        free ( leftBuf );
        free ( rightBuf );
        return - ( readLength - leftLen - leftS ); /// trim right
    }
    else
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dS%s", readLength - rightLen - rightS, rightBuf ); }

        correctedPac = segmentEndPos + 1 + rightOff;
        free ( leftBuf );
        free ( rightBuf );
        return readLength - rightLen - rightS; // trim left
    }

}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
long long getChrAndPosWithBoundaryCheck ( SRAQueryInput * qInput, unsigned long long readLength, unsigned long long ambPos,
                                    unsigned long long * tp, unsigned int * chr_id, char ** buffer )
{
    unsigned long long segmentEndPos = getChrAndPos ( qInput, ambPos, tp, chr_id );
    HSP * hsp = qInput->AlgnmtIndex->hsp;
    unsigned long long chrEndPos = hsp->seqOffset[ *chr_id - 1 ].endPos;
    unsigned long long correctedPac;
    long long ret;

    if ( ret = BoundaryCheck ( ambPos, *chr_id, chrEndPos, readLength, segmentEndPos, hsp, correctedPac, buffer ) )
    {
        if ( correctedPac > ambPos ) // trimmed left side
        {
            getChrAndPos ( qInput, correctedPac, tp, chr_id );
        }
    }

    return ret;
}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
long long getChrAndPosWithBoundaryCheckDP ( SRAQueryInput * qInput, unsigned long long readLength, unsigned long long ambPos, char * cigar,
                                      unsigned long long * tp, unsigned int * chr_id, char ** buffer )
{
    unsigned long long segmentEndPos = getChrAndPos ( qInput, ambPos, tp, chr_id );
    HSP * hsp = qInput->AlgnmtIndex->hsp;
    unsigned long long chrEndPos = hsp->seqOffset[ *chr_id - 1 ].endPos;
    unsigned long long correctedPac;
    long long ret;

    if ( ret = BoundaryCheckDP ( ambPos, *chr_id, chrEndPos, readLength, cigar, segmentEndPos, hsp, correctedPac, buffer ) )
    {
        if ( correctedPac > ambPos ) // trimmed left side
        {
            getChrAndPos ( qInput, correctedPac, tp, chr_id );
        }
    }

    return ret;
}

void initializeSAMAlgnmt ( bam1_t * samAlgnmt, int readlen, char * queryName, unsigned char * querySeq, char * queryComment,
                           char * qualities, char strand, uint8_t * xazStr, int xazlen, char * cigar_str,
                           int isUnmapped, char * readGroup )
{
    int i;
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );      //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = readlen;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( queryName ) + 1; //SAM: length of the query name

    if ( isUnmapped )
    { samAlgnmt->core.qual = 0; }     //SAM: mapping quality
    else
    { samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE; }     //SAM: mapping quality

    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), queryName, strlen ( queryName ) + 1 ); //Name

    if ( isUnmapped )
    {
        samAlgnmt->core.n_cigar = 0;                     //SAM: number of CIGAR operations
    }
    else
    {
        if ( cigar_str != NULL )
        {
            AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );
        }
        else
        {
            samAlgnmt->core.n_cigar = 1;                     //SAM: number of CIGAR operations
            SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), readlen << BAM_CIGAR_SHIFT ); //CIGAR
        }
    }

    if ( strand == QUERY_NEG_STRAND )
    {
        // for negative strand
        if ( readlen % 2 == 1 )
        {
            for ( i = ( readlen - 1 ) / 2; i > 0 ; i-- )                                                    //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 - 1]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            uint8_t biChar = bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[0]]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }
        else
        {
            for ( i = readlen / 2 - 1; i >= 0 ; i-- )                                                      //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 + 1]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }
        }

        for ( i = readlen - 1; i >= 0 ; i-- )                                                          //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] - 33);
        }
    }
    else
    {
        // for positive strand
        for ( i = 0; i < readlen / 2; i++ )                                                        //Read
        {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2]]] << 4;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2 + 1]]];
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        if ( readlen % 2 == 1 )
        {
            uint8_t biChar = bam_nt16_table[dnaChar[querySeq[readlen - 1]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        for ( i = 0; i < readlen; i++ )                                                            //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] - 33);
        }
    }

    unsigned int auxStart = samAlgnmt->data_len;
    bam_aux_append ( samAlgnmt, "RG", 'Z', strlen ( readGroup ) + 1, ( uint8_t * ) readGroup );

    if ( xazlen > 0 )
    { bam_aux_append ( samAlgnmt, "XA", 'Z', xazlen + 1, xazStr ); }
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;

}

void initializeSAMAlgnmt2 ( bam1_t * samAlgnmt, int readlen, char * queryName, unsigned char * querySeq, char * queryComment,
                            char * qualities, char strand, uint8_t * xazStr, int xazlen, char * cigar_str, int isUnmapped,
                            int mismatchNum, int editDist, int bestHitNum, int secBestHitNum, int gapOpenNum, int gapExtendNum,
                            char * mdStr, int mdlen, int map_qual_score, char * readGroup, bool isPrintMDNM, int moduleId = -1 )
{
    int i;
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );      //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = readlen;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( queryName ) + 1; //SAM: length of the query name

    if ( isUnmapped )
    { samAlgnmt->core.qual = 0; }     //SAM: mapping quality
    else
    { samAlgnmt->core.qual = map_qual_score; }     //SAM: mapping quality

    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), queryName, strlen ( queryName ) + 1 ); //Name

    if ( isUnmapped )
    {
        samAlgnmt->core.n_cigar = 0;                     //SAM: number of CIGAR operations
    }
    else
    {
        if ( cigar_str != NULL )
        {
            AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );
        }
        else
        {
            samAlgnmt->core.n_cigar = 1;                     //SAM: number of CIGAR operations
            SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), readlen << BAM_CIGAR_SHIFT ); //CIGAR
        }
    }

    if ( strand == QUERY_NEG_STRAND )
    {
        // for negative strand
        if ( readlen % 2 == 1 )
        {
            for ( i = ( readlen - 1 ) / 2; i > 0 ; i-- )                                                    //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 - 1]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            uint8_t biChar = bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[0]]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }
        else
        {
            for ( i = readlen / 2 - 1; i >= 0 ; i-- )                                                      //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 + 1]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }
        }

        for ( i = readlen - 1; i >= 0 ; i-- )                                                          //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] - 33);
        }
    }
    else
    {
        // for positive strand
        for ( i = 0; i < readlen / 2; i++ )                                                        //Read
        {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2]]] << 4;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2 + 1]]];
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        if ( readlen % 2 == 1 )
        {
            uint8_t biChar = bam_nt16_table[dnaChar[querySeq[readlen - 1]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        for ( i = 0; i < readlen; i++ )                                                            //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] - 33);
        }
    }

    unsigned int auxStart = samAlgnmt->data_len;
    bam_aux_append ( samAlgnmt, "RG", 'Z', strlen ( readGroup ) + 1, ( uint8_t * ) readGroup );

    if ( !isUnmapped )
    {
        if ( isPrintMDNM && editDist >= 0 )
        {
            bam_aux_append ( samAlgnmt, "NM", 'i', 4, ( uint8_t * ) &editDist );
        }

        if ( bestHitNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "X0", 'i', 4, ( uint8_t * ) &bestHitNum );
        }

        if ( secBestHitNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "X1", 'i', 4, ( uint8_t * ) &secBestHitNum );
        }

        if ( mismatchNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XM", 'i', 4, ( uint8_t * ) &mismatchNum );
        }

        if ( gapOpenNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XO", 'i', 4, ( uint8_t * ) &gapOpenNum );
        }

        if ( gapExtendNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XG", 'i', 4, ( uint8_t * ) &gapExtendNum );
        }

        if ( isPrintMDNM && mdlen > 0 )
        {
            bam_aux_append ( samAlgnmt, "MD", 'Z', mdlen + 1, ( uint8_t * ) mdStr );
        }

        if ( xazlen > 0 )
        {
            bam_aux_append ( samAlgnmt, "XA", 'Z', xazlen + 1, xazStr );
        }

        if ( moduleId != -1 )
        {
            bam_aux_append ( samAlgnmt, "PH", 'i', 4, ( uint8_t * ) &moduleId );
        }
    }
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
}

int getMapQualScore ( int n, int mismatchNum, int avgMismatchQual, int maxMAPQ, int minMAPQ )
{
    // not considering the second best
    if ( n == 1 )
    {
        // n is 1
        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int bwaLikeSingleQualScore ( int x0, int x1, int * g_log_n )
{
    // fprintf(stderr, "x0 : %i; x1 : %i\n", x0, x1);
    int score;

    if ( x0 > 1 )
    { score = 0; }
    else if ( x1 == 0 )
    { score = 37; }
    else
    {
        x1 = ( x1 > 255 ) ? 255 : x1;
        int n = g_log_n[x1];
        score = ( 23 < n ) ? 0 : 23 - n;
    }

    // fprintf(stderr, "score : %i\n", score);
    return score;
}

int getMapQualScoreSingle ( int mismatchNum, int avgMismatchQual, int x0, int x1, int maxMAPQ, int minMAPQ, int isBWALike, int * g_log_n )
{
    if ( isBWALike )
    {
        return bwaLikeSingleQualScore ( x0, x1, g_log_n );
    }

    if ( x0 == 1 )
    {
        if ( x1 > 0 )
        {
            return minMAPQ;
        }

        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int getMapQualScoreForSingleDP ( int maxDPScore, int avgMismatchQual, int x0, int x1_t1, int x1_t2, int bestDPScore, int secondBestDPScore, int maxMAPQ, int minMAPQ, int dpThres, int isBWALike, int * g_log_n )
{
    if ( isBWALike )
    {
        return bwaLikeSingleQualScore ( x0, x1_t1 + x1_t2, g_log_n );
    }

    // considersation of uniqueness
    if ( x0 > 1 || x1_t1 > 0 )
    { return minMAPQ; }

    float R1, R2, R3, P;

    // consideration of suboptimal score (i.e. secondBestDPScore)
    if ( x1_t2 > 0 )
    {
        R1 = 1.0 - ( ( float ) ( secondBestDPScore - dpThres ) ) / ( 0.7 * bestDPScore - dpThres );
    }
    else
    {
        R1 = 1.0;
    };

    // consideration of X1 (i.e. x1_t1 + x1_t2)
    int x1 = x1_t1 + x1_t2;

    R2 = ( x1 > 100 ) ? penalty_ratio_x1[100] : penalty_ratio_x1[x1];

    // consideration of best DP score
    R3 = ( ( float ) ( bestDPScore - dpThres ) ) / ( maxDPScore - dpThres );

    // consideration of average of quality value in mismatch positions
    if ( avgMismatchQual < 0 ) { avgMismatchQual = 0; }
    else if ( avgMismatchQual > 40 ) { avgMismatchQual = 40; }

    P = penalty_score_avg_mis_qual[avgMismatchQual];
    // calculation of the MAPQ score
    int mapqScore = ( int ) ( maxMAPQ * R1 * R2 * R3 - P );

    if ( mapqScore < minMAPQ ) { mapqScore = minMAPQ; }

    return mapqScore;
}


void bwaLikePairQualScore ( int x0_0, int x1_0, int x0_1, int x1_1, int * g_log_n, int op_score, int op_num, int subop_score, int subop_num, int readlen_0, int readlen_1, int * map_score0, int * map_score1 )
{
    // fprintf(stderr, "x0_0:%i x1_0:%i x0_1:%i x1_1:%i op_score:%i op_num:%i subop_score:%i subop_num:%i readlen_0:%i readlen_1:%i \n", x0_0, x1_0, x0_1, x1_1, op_score, op_num, subop_score, subop_num, readlen_0, readlen_1);
    int mapq0 = bwaLikeSingleQualScore ( x0_0, x1_0, g_log_n );
    int mapq1 = bwaLikeSingleQualScore ( x0_1, x1_1, g_log_n );
    // fprintf(stderr, "mapq0 : %i ; mapq1 : %i\n", mapq0, mapq1);
    op_score = op_score * 10;
    subop_score = subop_score * 10;
    int mapq_p = 0;

    if ( mapq0 > 0 && mapq1 > 0 )
    {
        mapq_p = mapq0 + mapq1;

        if ( mapq_p > 60 ) { mapq_p = 60; }

        mapq0 = mapq_p;
        mapq1 = mapq_p;
    }
    else
    {
        if ( op_num == 1 )
        {
            if ( subop_num == 0 )
            { mapq_p = 29; }
            else if ( op_score - subop_score > ( 0.3 * ( ( readlen_0 + readlen_1 ) / 2 ) ) ) { mapq_p = 23; }
            else
            {
                subop_num = ( subop_num > 255 ) ? 255 : subop_num;
                mapq_p = ( op_score - subop_score ) / 2 - g_log_n[subop_num];

                if ( mapq_p < 0 ) { mapq_p = 0; }
            }
        }

        if ( mapq0 == 0 )
        {
            mapq0 = ( mapq_p + 7 < mapq1 ) ? mapq_p + 7 : mapq1;
        }

        if ( mapq1 == 0 )
        {
            mapq1 = ( mapq_p + 7 < mapq0 ) ? mapq_p + 7 : mapq0;
        }
    }

    ( *map_score0 ) = mapq0;
    ( *map_score1 ) = mapq1;
}

int getMapQualScore2 ( int mismatchNum, int avgMismatchQual, int x0, int x1, char isBestHit, unsigned int totalNumValidPairs, int maxMAPQ, int minMAPQ )
{
    if ( x0 == 1 && totalNumValidPairs == 1 )
    {
        if ( isBestHit == 0 && ( x1 > 1 ) )
        {
            return minMAPQ;
        }

        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}


int getMapQualScoreForDP ( int n, int dpScore, int maxDPScore, int avgMismatchQual, int maxMAPQ, int minMAPQ )
{
    // does not consider second best
    if ( n == 1 )
    {
        // n = 1
        int difference_index = 0;

        if ( dpScore < maxDPScore )
        { difference_index = ( int ) ( ( 1.0 - ( double ) dpScore / maxDPScore ) * 100.0 - 1.0 ) / 5 + 1; }

        if ( difference_index > 5 )
        { difference_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[difference_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int getMapQualScoreForDP2 ( int dpScore, int maxDPScore, int avgMismatchQual, int x0, int x1, int bestDPScore, int secondBestDPScore, char isBestHit, int totalNumValidPairs, int maxMAPQ, int minMAPQ )
{
    if ( x0 == 1 && totalNumValidPairs == 1 )
    {
        if ( isBestHit == 0 && ( x1 > 1 ) )
        {
            return minMAPQ;
        }

        int difference_index = 0;

        if ( dpScore < maxDPScore )
        {
            difference_index = ( int ) ( ( 1.0 - ( double ) dpScore / maxDPScore ) * 100.0 - 1.0 ) / 4 + 1;
        }

        if ( difference_index > 5 )
        { difference_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[difference_index][avg_mismatch_qual_index] );

        if ( bestDPScore > secondBestDPScore && ( ( double ) bestDPScore - secondBestDPScore ) / maxDPScore < 0.2 )
        { mapQualScore = minMAPQ; }

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}


int getMapQualScoreForPair ( int score1, int score2 )
{
    return ( score1 > score2 ) ? ( int ) ( score1 * 0.2 + score2 * 0.8 ) : ( int ) ( score1 * 0.8 + score2 * 0.2 );
}

#ifdef VC_LEFT_ALIGNMENT
char performLeftAlignment ( unsigned int * packedSeq, unsigned long long genomePos,
                            unsigned char * query, unsigned int readLength,
                            char strand, char * cigar, unsigned int cigarStrLen,
                            char * newCigar )
{
    if ( cigar == NULL )
    {
        sprintf ( newCigar, "%u%c", readLength, 'M' );
        return 1;
    }

    char isValid = 1;
    char validOpt[] = { 1, 0, 0, 1 }; // Valid to have { M I D S } for next operation.
#define setValidOpt(_M, _I, _D, _S) do { validOpt[0] = _M; validOpt[1] = _I; validOpt[2] = _D; validOpt[3] = _S; } while ( 0 );
    for ( unsigned int i = 0; i < cigarStrLen; ++i )
    {
        while ( cigar[i] <= '9' )
        {
            ++i;
        }

        if ( cigar[i] == 'S' )
        {
            if ( !validOpt[3] )
            {
                isValid = 0;
                break;
            }
            setValidOpt ( 1, 1, 1, 0 );
        }
        else if ( cigar[i] == 'M' )
        {
            if ( !validOpt[0] )
            {
                isValid = 0;
                break;
            }
            setValidOpt ( 0, 1, 1, 1 );
        }
        else if ( cigar[i] == 'I' )
        {
            if ( !validOpt[1] )
            {
                isValid = 0;
                break;
            }
            setValidOpt ( 1, 0, 0, 1 );
        }
        else if ( cigar[i] == 'D' )
        {
            if ( !validOpt[2] )
            {
                isValid = 0;
                break;
            }
            setValidOpt ( 1, 0, 0, 1 );
        }
        else
        {
            isValid = 0;
            break;
        }
    }
    if ( !isValid )
    {
        fprintf ( stderr, "Invalid Cigar String : %s\n", cigar );
        return 0;
    }

    unsigned int range[MAX_READ_LENGTH] = { 0 };
    unsigned char operation[MAX_READ_LENGTH] = { '\0' };
    unsigned int numOpt = 0;
    for ( unsigned int i = 0; i < cigarStrLen; ++i )
    {
        while ( cigar[i] <= '9' )
        {
            range[numOpt] = range[numOpt] * 10 + ( cigar[i] - '0' );
            ++i;
        }
        
        if ( numOpt > 0 && cigar[i] == 'S' )
        {
            range[numOpt + 1] = range[numOpt];
            operation[numOpt + 1] = cigar[i];
            
            range[numOpt] = 0;
            operation[numOpt] = 'M';
            
            numOpt += 2;
        }
        else
        {
            operation[numOpt++] = cigar[i];
        }
    }
    assert ( numOpt < MAX_READ_LENGTH );
    
    unsigned long long readPos = 0;
    unsigned long long gPos = genomePos;
    unsigned int opt = 0;
    while ( opt < numOpt )
    {
        if ( operation[opt] == 'I' && opt > 0 )
        {
            unsigned long long startPos = readPos;
            int searchRange = range[opt - 1];
            unsigned long long offset = range[opt];
            /////////////////////////////////////////////////////////////////////////////
            //             Left Aligned the Insertion                                  //
            /////////////////////////////////////////////////////////////////////////////
            for ( int i=1;i<=searchRange;++i )
            {
                unsigned long long posR = startPos + offset - 1;
                unsigned long long posL = startPos - 1;
                char baseBitR = strand ? query[posR] : 3 - query[readLength - posR - 1];
                char baseBitL = strand ? query[posL] : 3 - query[readLength - posL - 1];
                if ( baseBitR == baseBitL )
                {
                    --startPos;
                }
                else 
                {
                    break;
                }
            }
            /////////////////////////////////////////////////////////////////////////////
            //    Below commented code will perform Left Alinged and                   //
            //                              kept the insert pattern unchanged          //
            /////////////////////////////////////////////////////////////////////////////
            /*
            if ( offset > searchRange )
            {
                offset = searchRange;
            }
            
            unsigned int i = 0;
            while ( i < range[opt] )
            {
                char baseBitR = strand ? query[startPos + i]
                                : 3 - query[readLength - startPos - i - 1];
                char baseBitL = strand ? query[startPos + i - offset]
                                : 3 - query[readLength - startPos - i + offset - 1];
                char baseBitT = strand ? query[startPos + i + ( range[opt] - offset )]
                                : 3 - query[readLength - startPos - i - ( range[opt] - offset ) - 1];
                if ( baseBitL != baseBitR || baseBitL != baseBitT )
                {
                    offset--;
                    if ( offset == 0 )
                    {
                        break;
                    }
                    i = 0;
                }
                else
                {
                    ++i;
                    if ( i == range[opt] )
                    {
                        startPos -= offset;
                        searchRange -= offset;
                        if ( offset == range[opt] && searchRange > 0 )
                        {
                            if ( offset > searchRange )
                            {
                                offset = searchRange;
                            }
                            i = 0;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
            */
            if ( readPos - startPos > 0 )
            {
                unsigned long long realOffset = readPos - startPos;
                readPos -= realOffset;
                range[opt - 1] -= realOffset;
                range[opt + 1] += realOffset;
                if ( opt + 1 == numOpt )
                {
                    operation[opt + 1] = operation[opt - 1];
                    numOpt++;
                }
            }
            readPos += range[opt];
        }
        else if ( operation[opt] == 'D' && opt > 0 )
        {
            unsigned long long startPos = gPos;
            int searchRange = range[opt - 1];
            unsigned long long offset = range[opt];
            /////////////////////////////////////////////////////////////////////////////
            //             Left Aligned the Deletion                                   //
            /////////////////////////////////////////////////////////////////////////////
            for ( int i=1;i<=searchRange;++i )
            {
                unsigned long long posR = startPos + offset - 1;
                unsigned long long posL = startPos - 1;
                char baseBitR = ( char ) ( ( packedSeq[posR >> 4] >> ( 30 - ( ( posR & 15 ) << 1 ) ) ) & 3 );
                char baseBitL = ( char ) ( ( packedSeq[posL >> 4] >> ( 30 - ( ( posL & 15 ) << 1 ) ) ) & 3 );
                if ( baseBitR == baseBitL )
                {
                    --startPos;
                }
                else 
                {
                    break;
                }
            }

            /////////////////////////////////////////////////////////////////////////////
            //    Below commented code will perform Left Alinged and                   //
            //                              kept the delete pattern unchanged          //
            /////////////////////////////////////////////////////////////////////////////
            /*
            if ( offset > searchRange )
            {
                offset = searchRange;
            }
            
            unsigned int i = 0;
            while ( i < range[opt] )
            {
                unsigned int posR = startPos + i;
                unsigned int posL = posR - offset;
                unsigned int posT = posR + ( range[opt] - offset );
                char baseBitR = ( char ) ( ( packedSeq[posR >> 4] >> ( 30 - ( ( posR & 15 ) << 1 ) ) ) & 3 );
                char baseBitL = ( char ) ( ( packedSeq[posL >> 4] >> ( 30 - ( ( posL & 15 ) << 1 ) ) ) & 3 );
                char baseBitT = ( char ) ( ( packedSeq[posL >> 4] >> ( 30 - ( ( posT & 15 ) << 1 ) ) ) & 3 );
                if ( baseBitL != baseBitR || baseBitL != baseBitT )
                {
                    offset--;
                    if ( offset == 0 )
                    {
                        break;
                    }
                    i = 0;
                }
                else
                {
                    ++i;
                    if ( i == range[opt] )
                    {
                        startPos -= offset;
                        searchRange -= offset;
                        if ( offset == range[opt] && searchRange > 0 )
                        {
                            if ( offset > searchRange )
                            {
                                offset = searchRange;
                            }
                            i = 0;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
            */
            if ( gPos - startPos > 0 )
            {
                unsigned long long realOffset = gPos - startPos;
                gPos -= realOffset;
                readPos -= realOffset;
                range[opt - 1] -= realOffset;
                range[opt + 1] += realOffset;
                if ( opt + 1 == numOpt )
                {
                    operation[opt + 1] = operation[opt - 1];
                    numOpt++;
                }
            }
            
            gPos += range[opt];
        }
        else if ( operation[opt] == 'M' )
        {
            gPos += range[opt];
            readPos += range[opt];
        }
        else if ( operation[opt] == 'S' )
        {
            readPos += range[opt];
        }
        else
        {
            fprintf ( stderr, "Invalid cigar string %s\n", cigar );
            exit ( 1 );
        }
        
        opt++;
    }
    assert ( numOpt < MAX_READ_LENGTH );
    
    unsigned int curRange = range[0];
    char curOperation = operation[0];
    
    opt = 1;
    while ( opt < numOpt )
    {
        if ( range[opt] == 0 )
        {
            opt++;
            continue;
        }
        
        if ( operation[opt] == curOperation )
        {
            curRange+= range[opt];
        }
        else
        {
            sprintf ( newCigar + strlen ( newCigar ), "%u%c", curRange, curOperation );
            curRange = range[opt];
            curOperation = operation[opt];
        }
        opt++;
    }
    sprintf ( newCigar + strlen ( newCigar ), "%u%c", curRange, curOperation );
    assert ( strlen ( newCigar ) < MAX_READ_LENGTH );
}

#endif

static inline int decideTargetChr(SRAQueryInput *qInput, unsigned long long ambPosition, uint readLen, char *cigar = NULL) {
    unsigned long long posStart, posEnd;
    unsigned int chrStart, chrEnd;
    getChrAndPos(qInput, ambPosition, &posStart, &chrStart);
    getChrAndPos(qInput, ambPosition + readLen - 1, &posEnd, &chrEnd);

    if (chrStart == chrEnd) return chrStart;
    else if (cigar == NULL) return -1;

    // parse cigar string to get alignment length
    long long l = 0;
    long long alignmentLen = 0;
    char op;
    char *p = cigar;
    while (*p) {
        if (isdigit(*p)) l = l * 10 + (*p - '0');
        else {
            op = *p;
            if (op != 'I' && op != 'S') {
                alignmentLen += l;
            }
            l = 0;
        }
        ++p;
    }

    getChrAndPos(qInput, ambPosition + readLen - 1, &posEnd, &chrEnd);
    if (chrStart == chrEnd) return chrStart;
    else return -1;
}

struct MappingRecordFromHeader {
    int score;
    int s, e; // start and end position the mapping reference name in the string of "comment"
};

static inline int getMappingFromHeader(char *comment, vector<MappingRecordFromHeader> &v, double top_percentage, double scoreT = 0) {
    // return maximum mapping score
    v.clear();
    if (comment == NULL || strcmp(comment, "IGNORE") == 0) return 0;

    int score = atoi(comment + 6);
    if (score < scoreT) return score;
    else if (scoreT < score * top_percentage)
        scoreT = score * top_percentage;

    char* p = strchr(comment + 6, ';');
    while (*(p + 1) != '\0') {
        MappingRecordFromHeader m;
        m.score = atoi(p + 1);
        m.s = p + 1 - comment;
        p = strchr(p + 1, ';');
        m.e = p - comment;

        if (score >= scoreT)
            v.push_back(m);
    }

    return score;
}

static inline string outputFastqReadSeq(unsigned char *query, char *quality, int readLen, int qconst) {
    string ret;
    for (int i = 0; i < readLen; ++i) {
        ret += "ACGT"[query[i]];
    }
    string tmp;
    for (int i=0;i<readLen;i++) tmp += char((int)quality[i] + qconst);
    ret += "\n+\n" + tmp + "\n";
    return ret;
}

string unproperlypairDPOutputFastqAPI(SRAQueryInput * qInput, Algnmt * algn_list, int hitNum,
                                    unsigned char * query, char * qualities, int readLen, 
                                    char * queryName, char * queryComment)
{
    int qconst = qInput->AlgnmtIndex->hspaux->quality_constant;
    string ret = "@" + string(queryName);
    int bestScore = 0;
    double top_percentage = qInput->AlgnmtIndex->hspaux->top_percentage;
    
    if (qInput->AlgnmtIndex->hspaux->megapathMode == 2) {
        // i.e. not paired == not align
        hitNum = 0;
    }

    if (queryComment && strcmp(queryComment, "IGNORE") == 0) {
        return ret + "\tIGNORE\n" + outputFastqReadSeq(query, qualities, readLen, qconst);
    }

    vector<pair<int, int> > chrHits;

    for (int i = 0; i < hitNum; ++i) {
        int chr_id = decideTargetChr(qInput, algn_list[i].algnmt, readLen, algn_list[i].cigarString);
        if (chr_id < 0) {
            // set unaligned
            algn_list[i].algnmt = NOT_ALIGNED;
            algn_list[i].score = 0;
        } else {
            chrHits.push_back(make_pair(chr_id, -algn_list[i].score)); // use -score here for sorting it acending order
        }

        if (bestScore < algn_list[i].score) {
            bestScore = algn_list[i].score;
        }
    }

    sort(chrHits.begin(), chrHits.end());

    vector<MappingRecordFromHeader> v;
    int prevScore = getMappingFromHeader(queryComment, v, top_percentage, bestScore * top_percentage);
    if (prevScore > bestScore) {
        bestScore = prevScore;
    }

    // output scores
    ret += "\tSCORE:" + to_string((long long int)bestScore) + ";";
    if (bestScore > 0) {
        for (size_t i = 0; i < chrHits.size(); ++i) {
            if (i > 0 && chrHits[i].first == chrHits[i-1].first) {
                continue;
            }
            if (-chrHits[i].second > 0 && -chrHits[i].second >= bestScore * top_percentage) {
                ret += to_string(-(long long int)chrHits[i].second) + "," + qInput->AlgnmtIndex->hsp->annotation[chrHits[i].first-1].text + ";";
            }
        }
    }

    for (int i = 0; i < v.size(); ++i) {
        if (v[i].score >= bestScore * top_percentage) {
            ret += string(queryComment + v[i].s, queryComment + v[i].e) + ";";
        }
    }
    return ret + "\n" + outputFastqReadSeq(query, qualities, readLen, qconst);
}

void unproperlypairDPOutputSAMAPI ( SRAQueryInput * qInput, unsigned int firstReadID, Algnmt * algn_list1, unsigned int secondReadID, Algnmt * algn_list2,
                                    int hitNum1, int hitNum2,
                                    unsigned char * query1, unsigned char * query2,
                                    char * qualities1, char * qualities2,
                                    int readlen1, int readlen2,
                                    char * queryName1, char * queryName2,
                                    char * queryComment1, char * queryComment2,
                                    DynamicUint8Array * xazArray, Algnmt ** bestAlgns )
{
    if (qInput->AlgnmtIndex->hspaux->megapathMode) {
        string fq1 = unproperlypairDPOutputFastqAPI(qInput, algn_list1, hitNum1, query1, qualities1, readlen1, queryName1, queryComment1);
        string fq2 = unproperlypairDPOutputFastqAPI(qInput, algn_list2, hitNum2, query2, qualities2, readlen2, queryName2, queryComment2);
        fputs((fq1 + fq2).c_str(), stdout);
    }

    //printf("%s\n",__func__);
    // output the DP alignments which are not properly paired
    // since the alignments are not properly paired, it follows the single-end scoring scheme
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    
    OCC * occ = qSetting->occ;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    if (!samFilePtr) return;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned int chr_1 = 0;
    unsigned int chr_2 = 0;
    unsigned int samFlag;
    unsigned int i;
    unsigned long long currAmbPos;
    unsigned long long currTP = 0;
    unsigned int currChr = 0;
    char currStrand;
    int currStrLen = 0;
    unsigned int currScore;
    char currOccStr[500];
    char * currSpCigar = NULL;
    int currEditDist = 0;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    long long mdStrLen1 = 0;
    long long mdStrLen2 = 0;
    char cigarStr1[MAX_READ_LENGTH + 1];
    char cigarStr2[MAX_READ_LENGTH + 1];
    int cigarStrLen1 = 0;
    int cigarStrLen2 = 0;
    int map_qual_score1 = 0;
    int map_qual_score2 = 0;
    int bestMismatchNum1 = 0;
    int bestGapOpen1 = 0;
    int bestGapExtend1 = 0;
    int avg_mismatch_qual1;
    int bestMismatchNum2 = 0;
    int bestGapOpen2 = 0;
    int bestGapExtend2 = 0;
    int avg_mismatch_qual2;
    // obtain the best hit
    Algnmt * bestAlgn1 = NULL;
    int bestScore1 = 0;
    int bestScoreNum1 = 0;
    long long boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int deletedEnd1 = 0, deletedEnd2 = 0;

    if ( hitNum1 > 0 )
    {
        bestAlgn1 = & ( algn_list1[0] );
        bestScore1 = bestAlgn1->score;
        bestScoreNum1 = 1;

        for ( i = 1; i < hitNum1; i++ )
        {
            if ( algn_list1[i].score > bestScore1 )
            {
                bestAlgn1 = & ( algn_list1[i] );
                bestScore1 = bestAlgn1->score;
                bestScoreNum1 = 1;
            }
            else if ( algn_list1[i].score == bestScore1 )
            {
                bestScoreNum1++;
            }
        }
    }

    Algnmt * bestAlgn2 = NULL;
    int bestScore2 = 0;
    int bestScoreNum2 = 0;

    if ( hitNum2 > 0 )
    {
        bestAlgn2 = & ( algn_list2[0] );
        bestScore2 = bestAlgn2->score;
        bestScoreNum2 = 1;

        for ( i = 1; i < hitNum2; i++ )
        {
            if ( algn_list2[i].score > bestScore2 )
            {
                bestAlgn2 = & ( algn_list2[i] );
                bestScore2 = bestAlgn2->score;
                bestScoreNum2 = 1;
            }
            else if ( algn_list2[i].score == bestScore2 )
            {
                bestScoreNum2++;
            }
        }
    }

    // obtain x1_t1 and x1_t2, and secondBestScore
    int x1_t1_1 = 0;
    int x1_t2_1 = 0;
    int secondBestScore1 = -9999;
    int subopt_class_thres1 = ( int ) ( 0.7 * bestScore1 );

    for ( i = 0; i < hitNum1; i++ )
    {
        if ( algn_list1[i].score < bestScore1 )
        {
            if ( algn_list1[i].score > secondBestScore1 )
            {
                secondBestScore1 = algn_list1[i].score;
            }

            if ( algn_list1[i].score >= subopt_class_thres1 )
            {
                x1_t1_1++;
            }
            else
            {
                x1_t2_1++;
            }
        }
    }

    int secondBestNum1;

    if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
    { secondBestNum1 = -1; }
    else
    { secondBestNum1 = x1_t1_1 + x1_t2_1; }

    int x1_t1_2 = 0;
    int x1_t2_2 = 0;
    int secondBestScore2 = -9999;
    int subopt_class_thres2 = ( int ) ( 0.7 * bestScore2 );

    for ( i = 0; i < hitNum2; i++ )
    {
        if ( algn_list2[i].score < bestScore2 )
        {
            if ( algn_list2[i].score > secondBestScore2 )
            {
                secondBestScore2 = algn_list2[i].score;
            }

            if ( algn_list2[i].score >= subopt_class_thres2 )
            {
                x1_t1_2++;
            }
            else
            {
                x1_t2_2++;
            }
        }
    }

    int secondBestNum2;

    if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
    { secondBestNum2 = -1; }
    else
    { secondBestNum2 = x1_t1_2 + x1_t2_2; }

    if ( ( bestAlgn1 != NULL ) && bestAlgn1->algnmt != NOT_ALIGNED && ( hspaux->alignmentType != OUTPUT_UNIQUE_BEST || bestScoreNum1 == 1 ) )
    {

        if ( bestAlgn1->isFromDP == 1 )
        {
            boundTrim1 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen1, bestAlgn1->algnmt, bestAlgn1->cigarString, &tp_1, &chr_1, &newCigar1 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen1 = boundTrim1 ? convertToCigarStr ( newCigar1, cigarStr1 ) : convertToCigarStr ( bestAlgn1->cigarString, cigarStr1, &deletedEnd1 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen1 = getMisInfoForDP ( hsp, query1, qualities1, readlen1, bestAlgn1->algnmt, bestAlgn1->strand,
                                          boundTrim1 ? newCigar1 : bestAlgn1->cigarString, mdStr1, &bestMismatchNum1, &bestGapOpen1, &bestGapExtend1,
                                          &avg_mismatch_qual1, boundTrim1 );
        }
        else
        {
            boundTrim1 = getChrAndPosWithBoundaryCheck ( qInput, readlen1, bestAlgn1->algnmt, &tp_1, &chr_1, &newCigar1 );
            // to get the md str
            mdStrLen1 = getMdStr ( hsp, query1, qualities1, readlen1, bestAlgn1->algnmt, bestAlgn1->strand, bestAlgn1->editdist, mdStr1, &avg_mismatch_qual1, boundTrim1 );
            bestMismatchNum1 = 0;

            for ( int ii = 0; ii < mdStrLen1; ++ii )
            { bestMismatchNum1 += ( mdStr1[ii] > '9' ); } // A C G or T
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score1 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {
            map_qual_score1 = getMapQualScoreForSingleDP ( readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestScoreNum1, x1_t1_1, x1_t2_1, bestScore1, secondBestScore1, hspaux->maxMAPQ, hspaux->minMAPQ, std::max(hspaux->singleDPcutoffRatio * readlen1, hspaux->singleDPcutoffLB) , hspaux->bwaLikeScore, hspaux->g_log_n );

            if ( !hspaux->bwaLikeScore )
            { map_qual_score1 = map_qual_score1 >> 1; }

            if ( map_qual_score1 < hspaux->minMAPQ )
            {
                map_qual_score1 = hspaux->minMAPQ;
            }
        }

        if ( boundTrim1 ) { map_qual_score1 = 0; }
    }
    else
    {
        bestAlgn1 = NULL;
    }

    if ( ( bestAlgn2 != NULL ) && bestAlgn2->algnmt != NOT_ALIGNED && ( hspaux->alignmentType != OUTPUT_UNIQUE_BEST || bestScoreNum2 == 1 ) )
    {

        if ( bestAlgn2->isFromDP == 1 )
        {
            boundTrim2 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen2, bestAlgn2->algnmt, bestAlgn2->cigarString, &tp_2, &chr_2, &newCigar2 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen2 = boundTrim2 ? convertToCigarStr ( newCigar2, cigarStr2 ) : convertToCigarStr ( bestAlgn2->cigarString, cigarStr2, &deletedEnd2 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen2 = getMisInfoForDP ( hsp, query2, qualities2, readlen2, bestAlgn2->algnmt, bestAlgn2->strand,
                                          boundTrim2 ? newCigar2 : bestAlgn2->cigarString, mdStr2, &bestMismatchNum2, &bestGapOpen2, &bestGapExtend2,
                                          &avg_mismatch_qual2, boundTrim2 );

        }
        else
        {
            boundTrim2 = getChrAndPosWithBoundaryCheck ( qInput, readlen2, bestAlgn2->algnmt, &tp_2, &chr_2, &newCigar2 );
            // to get the md str
            mdStrLen2 = getMdStr ( hsp, query2, qualities2, readlen2, bestAlgn2->algnmt, bestAlgn2->strand, bestAlgn2->editdist, mdStr2, &avg_mismatch_qual2, boundTrim2 );
            bestMismatchNum2 = 0;

            for ( int ii = 0; ii < mdStrLen2; ++ii )
            { bestMismatchNum2 += ( mdStr2[ii] > '9' ); } // A C G or T
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score2 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {
            map_qual_score2 = getMapQualScoreForSingleDP ( readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestScoreNum2, x1_t1_2, x1_t2_2, bestScore2, secondBestScore2, hspaux->maxMAPQ, hspaux->minMAPQ, std::max(hspaux->singleDPcutoffRatio * readlen2, hspaux->singleDPcutoffLB), hspaux->bwaLikeScore, hspaux->g_log_n );

            if ( !hspaux->bwaLikeScore )
            { map_qual_score2 = map_qual_score2 >> 1; }

            if ( map_qual_score2 < hspaux->minMAPQ )
            {
                map_qual_score2 = hspaux->minMAPQ;
            }

            if ( boundTrim2 ) { map_qual_score2 = 0; }
        }
    }
    else
    {
        bestAlgn2 = NULL;
    }

    //-------------------------------------------//
    // report the first alignment                //
    //-------------------------------------------//
    if ( hspaux->alignmentType == OUTPUT_ALL_BEST || hspaux->alignmentType == OUTPUT_ALL_VALID )
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );

        for ( i = 0; i < hitNum1; i++ )
        {
            if ( bestAlgn1 == & ( algn_list1[i] ) )
            { continue; }

            if ( hspaux->alignmentType == OUTPUT_ALL_BEST && algn_list1[i].score < bestScore1 )
            { continue; }

            currScore = algn_list1[i].score;
            currStrand = algn_list1[i].strand;
            currSpCigar = algn_list1[i].cigarString;
            currEditDist = algn_list1[i].editdist;
            currAmbPos = algn_list1[i].algnmt;
            if ( algn_list1[i].algnmt == NOT_ALIGNED) continue; 
            getChrAndPos ( qInput, algn_list1[i].algnmt,
                           &currTP, &currChr );
            char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
            memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
            int pos = strlen ( chr_name );
            currOccStr[pos++] = ',';
            currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
            pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
            currOccStr[pos++] = ',';
            pos += convertToCigarStr ( currSpCigar, & ( currOccStr[pos] ) ); // cigar string
            currOccStr[pos++] = ',';
            pos += writeNumToStr ( currEditDist, & ( currOccStr[pos] ) ); // edit distance
            currOccStr[pos++] = ';';
            currOccStr[pos] = '\0';
            currStrLen = pos;
            appendStringToUint8Array ( xazArray, currOccStr, currStrLen );
        }

        if ( bestAlgn1 != NULL )
            initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                                   qualities1, bestAlgn1->strand, xazArray->charStr, xazArray->length, bestAlgn1->isFromDP ? cigarStr1 : newCigar1, 0,
                                   bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestScoreNum1, secondBestNum1, bestGapOpen1, bestGapExtend1, mdStr1, mdStrLen1, map_qual_score1, hspaux->readGroup,
                                   hspaux->isPrintMDNM, hspaux->dpStageId );
        else
            initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                                  qualities1, QUERY_POS_STRAND, xazArray->charStr,
                                  xazArray->length, NULL, 1, hspaux->readGroup );

        // compute the value of the flag
        samFlag = 1;

        if ( bestAlgn1 == NULL )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestAlgn2 == NULL )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        samFlag |= SAM_FLAG_FIRST_IN_PAIR;

        if ( bestAlgn1 != NULL && bestAlgn1->strand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( bestAlgn2 != NULL && bestAlgn2->strand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
        samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

        if ( chr_1 > 0 && chr_1 == chr_2 )
            if ( tp_2 > tp_1 )
            { samAlgnmt->core.isize = tp_2 - deletedEnd2 + readlen2 - tp_1; }
            else
            { samAlgnmt->core.isize = - ( tp_1 - deletedEnd1 + readlen1 - tp_2 ); }
        else
        { samAlgnmt->core.isize = 0; }

        samwrite ( samFilePtr, samAlgnmt );
        
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        if ( hspaux->alignmentType == OUTPUT_ALL_BEST || hspaux->alignmentType == OUTPUT_ALL_VALID )
        {
            // Build the XA:Z tag
            DynamicUint8ArrayReset ( xazArray );

            for ( i = 0; i < hitNum2; i++ )
            {
                if ( bestAlgn2 == & ( algn_list2[i] ) )
                { continue; }

                if ( hspaux->alignmentType == OUTPUT_ALL_BEST && algn_list2[i].score < bestScore2 )
                { continue; }

                currScore = algn_list2[i].score;
                currStrand = algn_list2[i].strand;
                currSpCigar = algn_list2[i].cigarString;
                currEditDist = algn_list2[i].editdist;
                currAmbPos = algn_list2[i].algnmt;
                if ( algn_list2[i].algnmt == NOT_ALIGNED) continue;
                getChrAndPos ( qInput, algn_list2[i].algnmt, &currTP, &currChr );
                char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
                memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                currOccStr[pos++] = ',';
                currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
                currOccStr[pos++] = ',';
                pos += convertToCigarStr ( currSpCigar, & ( currOccStr[pos] ) ); // cigar string
                currOccStr[pos++] = ',';
                pos += writeNumToStr ( currEditDist, & ( currOccStr[pos] ) ); // edit distance
                currOccStr[pos++] = ';';
                currOccStr[pos] = '\0';
                currStrLen = pos;
                appendStringToUint8Array ( xazArray, currOccStr, currStrLen );
            }
        }

        if ( bestAlgn2 != NULL )
            initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                                   qualities2, bestAlgn2->strand, xazArray->charStr, xazArray->length, bestAlgn2->isFromDP ? cigarStr2 : newCigar2, 0,
                                   bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestScoreNum2, secondBestNum2, bestGapOpen2, bestGapExtend2, mdStr2, mdStrLen2, map_qual_score2, hspaux->readGroup,
                                   hspaux->isPrintMDNM, hspaux->dpStageId );
        else
            initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                                  qualities2, QUERY_POS_STRAND, xazArray->charStr,
                                  xazArray->length, NULL, 1, hspaux->readGroup );

        // compute the value of the flag
        samFlag = 1;

        if ( bestAlgn2 == NULL )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestAlgn1 == NULL )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        samFlag |= SAM_FLAG_SECOND_IN_PAIR;

        if ( bestAlgn1 != NULL && bestAlgn1->strand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        if ( bestAlgn2 != NULL && bestAlgn2->strand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
        samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

        if ( chr_1 > 0 && chr_1 == chr_2 )
            if ( tp_1 > tp_2 )
            { samAlgnmt->core.isize = tp_1 - deletedEnd1 + readlen1 - tp_2; }
            else
            { samAlgnmt->core.isize = - ( tp_2 - deletedEnd2 + readlen2 - tp_1 ); }
        else
        { samAlgnmt->core.isize = 0; }

        samwrite ( samFilePtr, samAlgnmt );
    }
    
    bestAlgns[0] = bestAlgn1;
    bestAlgns[1] = bestAlgn2;
    
    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }

}

int readLengthWithCigar ( char * cigar )
{
    char * p = cigar;
    char op = ' ';
    int x = 0;
    int len = 0;

    while ( *p )
    {
        x = 0;

        while ( *p <= '9' )
        {
            x = x * 10 + ( *p - '0' );
            p++;
        }

        op = *p;
        p++;

        if ( op == 'M' || op == 'm' || op == 'D' ) { len += x; }
    }

    if ( op == 'D' ) // last delete
    { len -= x; } // ignore => undo

    return len;
}

// if r1 and r2 are paired, then r1.score = r2.score = sum of the two score
// else r1 and r2 unchange
// if r is ummap, r.score = 0
void normalizeScore(DeepDPAlignResult * algnResult, bool paired) {
    if (algnResult->algnmt_1 == NOT_ALIGNED) {
        algnResult->score_1 = 0;
    }

    if (algnResult->algnmt_2 == NOT_ALIGNED) {
        algnResult->score_2 = 0;
    }

    if (paired && algnResult->algnmt_1 != NOT_ALIGNED && algnResult->algnmt_2 != NOT_ALIGNED) {
        int sumScore = algnResult->score_1 + algnResult->score_2;
        algnResult->score_1 = sumScore;
        algnResult->score_2 = sumScore;
    }
}

// output fastq, alignment results append to fastq header
string pairDeepDPOutputFastqAPI(SRAQueryInput * qInput, DeepDPAlignResult * algnResult,
                              DeepDPAlignResult * bestResult,
                              unsigned char * query1, unsigned char * query2,
                              char * qualities1, char * qualities2,
                              int readlen1, int readlen2,
                              char * queryName1, char * queryName2,
                              char * queryComment1, char * queryComment2,
                              unsigned int start, unsigned int num)
{
    int qconst = qInput->AlgnmtIndex->hspaux->quality_constant;
    int bestScore1 = 0, bestScore2 = 0;
    double top_percentage = qInput->AlgnmtIndex->hspaux->top_percentage;
    string ret;

    vector<pair<int, int> > chrHits1, chrHits2;

    for (int i = start; i < start + num; ++i) {
        int chr1 = -1, chr2 = -1;
        if (algnResult[i].algnmt_1 != NOT_ALIGNED) {
            chr1 = decideTargetChr(qInput, algnResult[i].algnmt_1, readlen1, algnResult[i].cigarString_1);
        }
        if (algnResult[i].algnmt_2 != NOT_ALIGNED) {
            chr2 = decideTargetChr(qInput, algnResult[i].algnmt_2, readlen2, algnResult[i].cigarString_2);
        }

        // if alignment on chr boundary, set as not aligned
        if (chr1 == -1) {
            algnResult[i].algnmt_1 = NOT_ALIGNED;
            algnResult[i].score_1 = 0;
        }
        if (chr2 == -1) {
            algnResult[i].algnmt_2 = NOT_ALIGNED;
            algnResult[i].score_2 = 0;
        }

        if (qInput->AlgnmtIndex->hspaux->megapathMode == 2) {
            // not paired == not align
            if (algnResult[i].algnmt_1 == NOT_ALIGNED || algnResult[i].algnmt_2 == NOT_ALIGNED) {
                algnResult[i].algnmt_1 = NOT_ALIGNED;
                algnResult[i].score_1 = 0;
                algnResult[i].algnmt_2 = NOT_ALIGNED;
                algnResult[i].score_2 = 0;
            }
        }

        // update bestScore
        normalizeScore(algnResult + i, chr1 == chr2);

        if (bestScore1 < algnResult[i].score_1) {
            bestScore1 = algnResult[i].score_1;
        }

        if (bestScore2 < algnResult[i].score_2) {
            bestScore2 = algnResult[i].score_2;
        }

        if (chr1 != -1) {
            chrHits1.push_back(make_pair(chr1, -algnResult[i].score_1));
        }
        if (chr2 != -1) {
            chrHits2.push_back(make_pair(chr2, -algnResult[i].score_2));
        }
    }

    sort(chrHits1.begin(), chrHits1.end());
    sort(chrHits2.begin(), chrHits2.end());

    vector<MappingRecordFromHeader> v1, v2;
    int prevScore1 = getMappingFromHeader(queryComment1, v1, top_percentage, bestScore1 * top_percentage);
    int prevScore2 = getMappingFromHeader(queryComment2, v2, top_percentage, bestScore2 * top_percentage);

    if (prevScore1 > bestScore1) bestScore1 = prevScore1;
    if (prevScore2 > bestScore2) bestScore2 = prevScore2;

    // output read name
    ret = "@" + string(queryName1);
    // output scores
    if (queryComment1 && strcmp(queryComment1, "IGNORE") == 0) {
        ret += "\tIGNORE\n";
    } else {
        ret += "\tSCORE:" + to_string((long long int)bestScore1) + ";";
        if (bestScore1 > 0) {
            for (size_t i = 0; i < chrHits1.size(); ++i) {
                if (i > 0 && chrHits1[i].first == chrHits1[i-1].first) continue;
                if (-chrHits1[i].second > 0 && -chrHits1[i].second >= bestScore1 * top_percentage) {
                    ret += to_string(-(long long int)chrHits1[i].second) + "," + qInput->AlgnmtIndex->hsp->annotation[chrHits1[i].first-1].text + ";";
                }
            }
        }
        for (int i = 0; i < v1.size(); ++i) {
            if (v1[i].score >= bestScore1 * top_percentage) {
                ret += string(queryComment1 + v1[i].s, queryComment1 + v1[i].e);
                ret += ";";
            }
        }
        ret += "\n";
    }
    // output seq & quality
    ret += outputFastqReadSeq(query1, qualities1, readlen1, qconst);

    ret += "@" + string(queryName2);
    // output scores
    if (queryComment2 && strcmp(queryComment2, "IGNORE") == 0) {
        ret += "\tIGNORE\n";
    } else {
        ret += "\tSCORE:" + to_string((long long int)bestScore2) + ";";
        if (bestScore2 > 0) {
            for (size_t i = 0; i < chrHits2.size(); ++i) {
                if (i > 0 && chrHits2[i].first == chrHits2[i-1].first) continue;
                if (-chrHits2[i].second > 0 && -chrHits2[i].second >= bestScore2 * top_percentage) {
                    ret += to_string(-(long long int)chrHits2[i].second) + "," + qInput->AlgnmtIndex->hsp->annotation[chrHits2[i].first-1].text + ";";
                }
            }
        }
        for (int i = 0; i < v2.size(); ++i) {
            if (v2[i].score >= bestScore2 * top_percentage) {
                ret += string(queryComment2 + v2[i].s, queryComment2 + v2[i].e);
                ret += ";";
            }
        }
        ret += "\n";
    }
    // output seq & quality
    ret += outputFastqReadSeq(query2, qualities2, readlen2, qconst);
    return ret;
}

void pairDeepDPOutputSAMAPI ( SRAQueryInput * qInput, DeepDPAlignResult * algnResult,
                              DeepDPAlignResult * bestResult,
                              unsigned int start, unsigned int num,
                              unsigned char * query1, unsigned char * query2,
                              char * qualities1, char * qualities2,
                              int readlen1, int readlen2,
                              char * queryName1, char * queryName2,
                              char * queryComment1, char * queryComment2,
                              DynamicUint8Array * xazArray, char ** twoSamStrings, int threadId )
{
    if (qInput->AlgnmtIndex->hspaux->megapathMode) {
        string fq2interleaved = pairDeepDPOutputFastqAPI(qInput, algnResult, bestResult, query1, query2,
                                 qualities1, qualities2, readlen1, readlen2,
                                 queryName1, queryName2,
                                 queryComment1, queryComment2,
                                 start, num);
        fputs(fq2interleaved.c_str(), stdout);
    }
    //printf("%s\n",__func__);
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    
    OCC * occ = qSetting->occ;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    if (!samFilePtr) return;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned long long curr_tp = 0;
    unsigned int chr_1 = 0;
    unsigned int chr_2 = 0;
    unsigned int curr_chr = 0;
    unsigned int samFlag;
    unsigned int i;
    char curr_strand;
    int curr_len = 0;
    char curr_occStr[500];
    char best_strand1 = 1;
    char best_strand2 = 1;
    char * best_cigar1 = NULL;
    char * best_cigar2 = NULL;
    int best_insert = 0;
    // int mapping_qual_score = 0; // for a valid paired-end alignment
    int mapq1, mapq2;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen1 = 0;
    char cigarStr1[MAX_READ_LENGTH + 1];
    int cigarStrLen1 = 0;
    int bestMismatchNum1 = 0;
    int bestGapOpen1 = 0;
    int bestGapExtend1 = 0;
    int avg_mismatch_qual1;
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen2 = 0;
    char cigarStr2[MAX_READ_LENGTH + 1];
    int cigarStrLen2 = 0;
    int bestMismatchNum2 = 0;
    int bestGapOpen2 = 0;
    int bestGapExtend2 = 0;
    int avg_mismatch_qual2;

    long long boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int r1 = readlen1, r2 = readlen2;

    if ( bestResult != NULL )
    {
        // for the reads which do not have any hit for both ends
        if ( bestResult->algnmt_1 == NOT_ALIGNED && bestResult->algnmt_2 == NOT_ALIGNED )
        {
            return;
        }

        // collect the information of the first alignment
        if ( bestResult->algnmt_1 != NOT_ALIGNED )
        {
            best_strand1 = bestResult->strand_1;
            best_cigar1 = bestResult->cigarString_1;
        }

        // collect the information of the second alignment
        if ( bestResult->algnmt_2 != NOT_ALIGNED ) {
            best_strand2 = bestResult->strand_2;
            best_cigar2 = bestResult->cigarString_2;
        }

        if ( bestResult->algnmt_1 != NOT_ALIGNED && bestResult->algnmt_2 != NOT_ALIGNED )
        {
            best_insert = bestResult->insertSize;
        }

        // obtain the corresponding normal cigar string, # of gap open, # of gap extension, # of mismatch,
        // and the average quality values in the mismatch positions
        if ( bestResult->algnmt_1 != NOT_ALIGNED )
        {
            boundTrim1 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen1, bestResult->algnmt_1, bestResult->cigarString_1, &tp_1, &chr_1, &newCigar1 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen1 = boundTrim1 ? convertToCigarStr ( newCigar1, cigarStr1 ) : convertToCigarStr ( bestResult->cigarString_1, cigarStr1 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen1 = getMisInfoForDP ( hsp, query1, qualities1, readlen1, bestResult->algnmt_1, best_strand1,
                                          boundTrim1 ? newCigar1 : bestResult->cigarString_1, mdStr1, &bestMismatchNum1, &bestGapOpen1, &bestGapExtend1, &avg_mismatch_qual1, boundTrim1 );
        }

        if ( bestResult->algnmt_2 != NOT_ALIGNED )
        {
            boundTrim2 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen2, bestResult->algnmt_2, bestResult->cigarString_2, &tp_2, &chr_2, &newCigar2 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen2 = boundTrim2 ? convertToCigarStr ( newCigar2, cigarStr2 ) : convertToCigarStr ( bestResult->cigarString_2, cigarStr2 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen2 = getMisInfoForDP ( hsp, query2, qualities2, readlen2, bestResult->algnmt_2, best_strand2,
                                          boundTrim2 ? newCigar2 : bestResult->cigarString_2, mdStr2, &bestMismatchNum2, &bestGapOpen2, &bestGapExtend2, &avg_mismatch_qual2, boundTrim2 );
        }

        if ( bestResult->algnmt_1 != NOT_ALIGNED )
        {
            r1 = readLengthWithCigar ( bestResult->cigarString_1 );
        }

        if ( bestResult->algnmt_2 != NOT_ALIGNED )
        {
            r2 = readLengthWithCigar ( bestResult->cigarString_2 );
        }

        // check for read-through reads
        // not well-defined, commented
        /*
        if ( bestResult->algnmt_1 != NOT_ALIGNED && bestResult->algnmt_2 != NOT_ALIGNED )
        {
            unsigned long long adjust1 = bestResult->algnmt_1 + ( boundTrim1 > 0 ? boundTrim1 : 0 );
            unsigned long long adjust2 = bestResult->algnmt_2 + ( boundTrim2 > 0 ? boundTrim2 : 0 );

            if ( ( bestResult->strand_1 == QUERY_POS_STRAND &&
                    ( adjust1 > adjust2 || adjust1 + r1 > adjust2 + r2 ) ) ||
                    ( bestResult->strand_1 == QUERY_NEG_STRAND &&
                      ( adjust2 > adjust1 || adjust2 + r2 > adjust1 + r1 ) ) )
            {
                if ( bestMismatchNum1 <= bestMismatchNum2 )
                {
                    bestResult->algnmt_2 = NOT_ALIGNED;
                    tp_2 = chr_2 = 0;
                }
                else
                {
                    bestResult->algnmt_1 = NOT_ALIGNED;
                    tp_1 = chr_1 = 0;
                }

                best_insert = 0;
            }
        }*/

        int bestPairNum = 0;
        int bestPairScore = 0;
        int secBestPairScore = 0;

        // obtain the number of pairs with maximum sum of scores
        if ( bestResult->algnmt_1 != NOT_ALIGNED && bestResult->algnmt_2 != NOT_ALIGNED )
        {
            bestPairNum = 1;
            bestPairScore = bestResult->score_1 + bestResult->score_2;

            if ( num > 1 )
            {
                for ( i = start; i < start + num; i++ )
                {
                    if ( & ( algnResult[i] ) == bestResult )
                    { continue; }

                    if ( algnResult[i].score_1 + algnResult[i].score_2 == bestPairScore )
                    { bestPairNum++; }
                    else if ( algnResult[i].score_1 + algnResult[i].score_2 > secBestPairScore )
                    { secBestPairScore = algnResult[i].score_1 + algnResult[i].score_2; }
                }
            }
        }

        // obtain the number of pairs with similar score to the best score (i.e. difference <= 1 mismatch score)
        int numSimilarBestPairs = 0;

        if ( bestResult->algnmt_1 != NOT_ALIGNED && bestResult->algnmt_2 != NOT_ALIGNED )
        {
            numSimilarBestPairs = 1; // includes the best result
            int bestPairScore1 = bestResult->score_1;
            int bestPairScore2 = bestResult->score_2;

            if ( num > 1 )
            {
                for ( i = start; i < start + num; i++ )
                {
                    if ( & ( algnResult[i] ) == bestResult )
                    { continue; }

                    if ( ( algnResult[i].score_1 >= ( bestPairScore1 + hspaux->dpMisMatchScore ) ) &&
                            ( algnResult[i].score_2 >= ( bestPairScore2 + hspaux->dpMisMatchScore ) ) )
                    {
                        numSimilarBestPairs++;
                    }
                }
            }
        }

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the first read
        int bestHitNum1 = 0;
        int secBestHitNum1 = 0;
        int bestScore1 = 0;
        int secBestScore1 = 0;
        char isBestHit1 = 1;
        unsigned long long bestPos1 = NOT_ALIGNED;
        unsigned long long secBestPos1 = NOT_ALIGNED;

        if ( bestResult->algnmt_1 != NOT_ALIGNED )
        {
            bestScore1 = bestResult->score_1;
            bestPos1 = bestResult->algnmt_1; // alignment position
            bestHitNum1 = 1;
            // bestHitNum1 = bestResult->num_sameScore_1;
        }

        if ( bestResult->algnmt_1 != NOT_ALIGNED && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].score_1 >= bestScore1 )
                {
                    if ( algnResult[i].score_1 == bestScore1 )
                    {
                        if ( algnResult[i].algnmt_1 != bestPos1 )
                            // bestHitNum1++;
                        { bestHitNum1 += algnResult[i].num_sameScore_1; }
                    }
                    else
                    {
                        secBestScore1 = bestScore1;
                        secBestHitNum1 = bestHitNum1;
                        secBestPos1 = bestPos1;
                        bestScore1 = algnResult[i].score_1;
                        // bestHitNum1 = 1;
                        bestHitNum1 = algnResult[i].num_sameScore_1;
                        bestPos1 = algnResult[i].algnmt_1;
                        isBestHit1 = 0;
                    }
                }
                else if ( algnResult[i].score_1 >= secBestScore1 )
                {
                    if ( algnResult[i].score_1 == secBestScore1 )
                    {
                        if ( algnResult[i].algnmt_1 != secBestPos1 )
                            // secBestHitNum1++;
                        { secBestHitNum1 += algnResult[i].num_sameScore_1; }
                    }
                    else
                    {
                        secBestScore1 = algnResult[i].score_1;
                        secBestPos1 = algnResult[i].algnmt_1;
                        // secBestHitNum1 = 1;
                        secBestHitNum1 = algnResult[i].num_sameScore_1;
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID] > 0 )
        {
            int x0_score = hspaux->mismatch_array[bestResult->readID] * hspaux->dpMisMatchScore
                           + ( readlen1 - hspaux->mismatch_array[bestResult->readID] ) * hspaux->dpMatchScore;

            if ( x0_score >= bestScore1 )
            {
                bestHitNum1 = ( bestHitNum1 > hspaux->x0_array[bestResult->readID] ) ? bestHitNum1 : hspaux->x0_array[bestResult->readID];
                secBestHitNum1 = ( secBestHitNum1 > hspaux->x1_array[bestResult->readID] ) ? secBestHitNum1 : hspaux->x1_array[bestResult->readID];

                if ( x0_score > bestScore1 )
                { isBestHit1 = 0; }
            }
        }

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the second read
        int bestHitNum2 = 0;
        int secBestHitNum2 = 0;
        int bestScore2 = 0;
        int secBestScore2 = 0;
        char isBestHit2 = 1;
        unsigned long long bestPos2 = NOT_ALIGNED;
        unsigned long long secBestPos2 = NOT_ALIGNED;

        if ( bestResult->algnmt_2 != NOT_ALIGNED )
        {
            bestScore2 = bestResult->score_2;
            bestPos2 = bestResult->algnmt_2; // alignment position
            bestHitNum2 = 1;
            // bestHitNum2 = bestResult->num_sameScore_2;
        }

        if ( bestResult->algnmt_2 != NOT_ALIGNED && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].score_2 >= bestScore2 )
                {
                    if ( algnResult[i].score_2 == bestScore2 )
                    {
                        if ( algnResult[i].algnmt_2 != bestPos2 )
                            // bestHitNum2++;
                        { bestHitNum2 += algnResult[i].num_sameScore_2; }
                    }
                    else
                    {
                        secBestScore2 = bestScore2;
                        secBestHitNum2 = bestHitNum2;
                        secBestPos2 = bestPos2;
                        bestScore2 = algnResult[i].score_2;
                        // bestHitNum2 = 1;
                        bestHitNum2 = algnResult[i].num_sameScore_2;
                        bestPos2 = algnResult[i].algnmt_2;
                        isBestHit2 = 0;
                    }
                }
                else if ( algnResult[i].score_2 >= secBestScore2 )
                {
                    if ( algnResult[i].score_2 == secBestScore2 )
                    {
                        if ( algnResult[i].algnmt_2 != secBestPos2 )
                            // secBestHitNum2++;
                        { secBestHitNum2 += algnResult[i].num_sameScore_2; }
                    }
                    else
                    {
                        secBestScore2 = algnResult[i].score_2;
                        secBestPos2 = algnResult[i].algnmt_2;
                        // secBestHitNum2 = 1;
                        secBestHitNum2 = algnResult[i].num_sameScore_2;
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID + 1] > 1 )
        {
            int x0_score = hspaux->mismatch_array[bestResult->readID + 1] * hspaux->dpMisMatchScore
                           + ( readlen2 - hspaux->mismatch_array[bestResult->readID + 1] ) * hspaux->dpMatchScore;

            if ( x0_score >= bestScore2 )
            {
                bestHitNum2 = ( bestHitNum2 > hspaux->x0_array[bestResult->readID + 1] ) ? bestHitNum2 : hspaux->x0_array[bestResult->readID + 1];
                secBestHitNum2 = ( secBestHitNum2 > hspaux->x1_array[bestResult->readID + 1] ) ? secBestHitNum2 : hspaux->x1_array[bestResult->readID + 1];

                if ( x0_score > bestScore2 )
                { isBestHit2 = 0; }
            }
        }

        if ( hspaux->alignmentType == OUTPUT_ALL_VALID || hspaux->alignmentType == OUTPUT_ALL_BEST )
        {
            if ( hspaux->bwaLikeScore )
            {
                bwaLikePairQualScore ( bestHitNum1, secBestHitNum1, bestHitNum2, secBestHitNum2, hspaux->g_log_n, bestPairScore, bestPairNum, secBestPairScore, num - bestPairNum, readlen1, readlen2, &mapq1, &mapq2 );
            }
            else
            {
                int mapping_qual_score1 = getMapQualScoreForDP2 ( bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestHitNum1, secBestHitNum1, bestScore1, secBestScore1, isBestHit1, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                int mapping_qual_score2 = getMapQualScoreForDP2 ( bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestHitNum2, secBestHitNum2, bestScore2, secBestScore2, isBestHit2, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                mapq1 = mapq2 = getMapQualScoreForPair ( mapping_qual_score1, mapping_qual_score2 );
            }

            if ( boundTrim1 ) { mapq1 = 0; }

            if ( boundTrim2 ) { mapq2 = 0; }

        }
        else
        {
            mapq1 = mapq2 = SAM_MAPQ_UNAVAILABLE;
        }

        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_1 != NOT_ALIGNED && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 + algnResult[i].score_2 < bestPairScore ) )
                {
                    continue;
                }
                if ( algnResult[i].algnmt_1 == NOT_ALIGNED) continue;
                getChrAndPos ( qInput, algnResult[i].algnmt_1, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_1;
                char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                curr_occStr[pos++] = ',';
                curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                curr_occStr[pos++] = ',';
                pos += convertToCigarStr ( algnResult[i].cigarString_1, & ( curr_occStr[pos] ) ); // cigar string
                curr_occStr[pos++] = ',';
                pos += writeNumToStr ( algnResult[i].editdist_1, & ( curr_occStr[pos] ) ); // edit distance
                curr_occStr[pos++] = ';';
                curr_occStr[pos] = '\0';
                curr_len = pos;
                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum1 = -1;
            secBestHitNum1 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum1 = -1;
        }

        if ( bestResult->algnmt_1 != NOT_ALIGNED )
        {
            if ( bestResult->algnmt_2 != NOT_ALIGNED )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length, cigarStr1, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1,
                                       bestGapExtend1, mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM, hspaux->dpStageId );
            }
            else
            {
                mapq1 = getMapQualScoreForDP ( bestHitNum1, bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, hspaux->maxMAPQ, hspaux->minMAPQ );

                if ( boundTrim1 ) { mapq1 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length, cigarStr1, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1,
                                       bestGapExtend1, mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM, hspaux->dpStageId );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                                  qualities1, best_strand1, NULL, 0, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != NOT_ALIGNED ) && ( bestResult->algnmt_2 != NOT_ALIGNED ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_FIRST_IN_PAIR;

        if ( bestResult->algnmt_1 == NOT_ALIGNED )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_2 == NOT_ALIGNED )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_1 != NOT_ALIGNED ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_2 != NOT_ALIGNED ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
        samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_1 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_2 - 1;
        //      samAlgnmt->core.mpos = tp_2 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_1 > tp_2 )
            { samAlgnmt->core.isize = - ( tp_1 + r1 - tp_2 ); }
            else
            { samAlgnmt->core.isize = tp_2 + r2 - tp_1; }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        if ( twoSamStrings != NULL )
            { twoSamStrings[0] = samString ( samFilePtr, samAlgnmt ); }
        else
            { samwrite ( samFilePtr, samAlgnmt ); }
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_2 != NOT_ALIGNED && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 + algnResult[i].score_2 < bestPairScore ) )
                {
                    continue;
                }

                if (algnResult[i].algnmt_2 == NOT_ALIGNED) continue;
                getChrAndPos ( qInput, algnResult[i].algnmt_2, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_2;
                char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                curr_occStr[pos++] = ',';
                curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                curr_occStr[pos++] = ',';
                pos += convertToCigarStr ( algnResult[i].cigarString_2, & ( curr_occStr[pos] ) ); // cigar string
                curr_occStr[pos++] = ',';
                pos += writeNumToStr ( algnResult[i].editdist_2, & ( curr_occStr[pos] ) ); // edit distance
                curr_occStr[pos++] = ';';
                curr_occStr[pos] = '\0';
                curr_len = pos;
                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum2 = -1;
            secBestHitNum2 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum2 = -1;
        }

        if ( bestResult->algnmt_2 != NOT_ALIGNED )
        {
            if ( bestResult->algnmt_1 != NOT_ALIGNED )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length, cigarStr2, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2,
                                       bestGapExtend2, mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM, hspaux->dpStageId );
            }
            else
            {
                mapq2 = getMapQualScoreForDP ( bestHitNum2, bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, hspaux->maxMAPQ, hspaux->minMAPQ );

                if ( boundTrim2 ) { mapq2 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length, cigarStr2, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2,
                                       bestGapExtend2, mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM, hspaux->dpStageId );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                                  qualities2, best_strand2, NULL, 0, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != NOT_ALIGNED ) && ( bestResult->algnmt_2 != NOT_ALIGNED ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_SECOND_IN_PAIR;

        if ( bestResult->algnmt_2 == NOT_ALIGNED )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_1 == NOT_ALIGNED )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_2 != NOT_ALIGNED ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_1 != NOT_ALIGNED ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
        samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_2 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_2 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_1 - 1;
        //      samAlgnmt->core.mpos = tp_1 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_2 > tp_1 )
            { samAlgnmt->core.isize = - ( tp_2 + r2 - tp_1 ); }
            else
            { samAlgnmt->core.isize = tp_1 + r1 - tp_2; }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        if ( twoSamStrings != NULL )
            { twoSamStrings[1] = samString ( samFilePtr, samAlgnmt ); }
        else
        {
            samwrite ( samFilePtr, samAlgnmt );
        }
    }
    else
    {
        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1, queryComment1,
                              qualities1, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        if ( twoSamStrings != NULL )
            { twoSamStrings[0] = samString ( samFilePtr, samAlgnmt ); }
        else
            { samwrite ( samFilePtr, samAlgnmt ); }
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2, queryComment2,
                              qualities2, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        if ( twoSamStrings != NULL )
            { twoSamStrings[1] = samString ( samFilePtr, samAlgnmt ); }
        else
        {
            samwrite ( samFilePtr, samAlgnmt );
        }
    }

    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }
    // printf("BGS-IO: %s\n%s\n", twoSamStrings[0] , twoSamStrings[1]);
}

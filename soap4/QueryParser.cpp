/*
 *
 *    QueryParser.c
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

#include "QueryParser.h"

// This function is to load reads for BAM format
int loadBAMReads ( bamFile bamQueryFile, bam_header_t * bamHeader, bam1_t * b, unsigned char * charMap,
                   uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char ** queryNames, char ** queryComments, uint maxReadLength, uint maxNumQueries,
                   uint wordPerQuery, int qualityConstant, int maxLenReadName )
{
    uint * queryPtr = queries;
    ullint queriesRead = 0;
    uint currentWord = 0;
    uint offset = 0;
    uint bits = 0;

    while ( queriesRead < maxNumQueries && bam_read1 ( bamQueryFile, b ) >= 0 )
    {
        // get the information of one read from BAM file
        char * rec_seq = bam_format1 ( bamHeader, b );
        // format: id >> flag >> chr >> pos >> mapQ >> cigar >> mateChr >> matePos >> insertSize >> seq >> qual;
        int rec_len = strlen ( rec_seq );
        int pos = 0;
        // store the read description
        int i = 0;

        while ( pos < rec_len && rec_seq[pos] != '\t' )
        {
            if ( i < maxLenReadName - 1 )
            {
                // TODO
                i++;
            }

            pos++;
        }

        // upkdQueryNames[queriesRead * maxLenReadName + i] = '\0';

        // skip the spaces
        while ( pos < rec_len && rec_seq[pos] == '\t' )
        { pos++; }

        // get the flag
        unsigned int flag = 0;

        while ( pos < rec_len && rec_seq[pos] != '\t' )
        {
            flag = flag * 10 + ( rec_seq[pos] - '0' );
            pos++;
        }

        // skip the spaces
        while ( pos < rec_len && rec_seq[pos] == '\t' )
        { pos++; }

        // skip seven columns
        // i.e. chr >> pos >> mapQ >> cigar >> mateChr >> matePos >> insertSize
        for ( i = 0; i < 7; i++ )
        {
            // skip the column
            while ( pos < rec_len && rec_seq[pos] != '\t' )
            { pos++; }

            // skip the spaces
            while ( pos < rec_len && rec_seq[pos] == '\t' )
            { pos++; }
        }

        // get the length of read sequence
        int seq_len = 0;

        while ( pos + seq_len < rec_len && rec_seq[pos + seq_len] != '\t' )
        { seq_len++; }

        // store the read sequence
        currentWord = 0;
        offset = 0;
        int isFWD = 1;

        if ( ( flag & SAM_FLAG_READ_ALGNMT_STRAND ) != 0 )
        {
            isFWD = 0; // reverse strand
        }

        if ( isFWD )
        {
            i = 0;

            while ( i < seq_len && i < maxReadLength )
            {
                bits = charMap[rec_seq[pos + i]];
                currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                offset++;
                i++;

                if ( offset == CHAR_PER_WORD )
                {
                    *queryPtr = currentWord;
                    queryPtr += 32;
                    offset = 0;
                    currentWord = 0;
                }
            }
        }
        else
        {
            i = seq_len - 1;

            while ( i >= 0 && seq_len - i <= maxReadLength )
            {
                bits = soap3DnaComplement[rec_seq[pos + i]];
                currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                offset++;
                i--;

                if ( offset == CHAR_PER_WORD )
                {
                    *queryPtr = currentWord;
                    queryPtr += 32;
                    offset = 0;
                    currentWord = 0;
                }
            }
        }

        pos += seq_len;

        if ( offset > 0 )
        { *queryPtr = currentWord; }

        readLengths[queriesRead] = ( seq_len > maxReadLength ) ? maxReadLength : seq_len;
        readIDs[queriesRead] = queriesRead + 1;

        // skip the spaces
        while ( pos < rec_len && rec_seq[pos] == '\t' )
        { pos++; }

        // store the quality values
        for ( i = 0; i < readLengths[queriesRead]; i++ )
        {
            if ( isFWD )
            { upkdQualities[queriesRead * maxReadLength + i] = rec_seq[pos + i] - qualityConstant; }
            else
            { upkdQualities[queriesRead * maxReadLength + i] = rec_seq[pos + seq_len - i - 1] - qualityConstant; }
        }

        queriesRead++;
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );
        free ( rec_seq );
    }

    return queriesRead;
}

// helper function trim /1 /2
static inline void trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
        s->l -= 2, s->s[s->l] = 0;
}

// helper function for loadPairReadsKseq
static inline void appendToQueryArrays(kseq_t *ks, unsigned char * charMap, uint * readLengths, uint * readIDs, uint * queries, char * upkdQualities, char ** queryNames, char ** queryComments,
                                uint maxReadLength, uint wordPerQuery, int qualityConstant, int queryID, char & isFastq)
{
    readIDs[queryID] = queryID + 1;
    readLengths[queryID] = ks->seq.l > maxReadLength - 1 ? maxReadLength - 1 : ks->seq.l;

    // sequence
    uint *queryPtr = queries + (queryID / 32 * 32 * wordPerQuery + queryID % 32);
    int offset = 0;
    uint word = 0;
    for (int i = 0; i < readLengths[queryID]; ++i) {
        word |= charMap[ks->seq.s[i]] << (offset * BIT_PER_CHAR);
        ++offset;
        if (offset == CHAR_PER_WORD) {
            offset = 0;
            *queryPtr = word;
            word = 0;
            queryPtr += 32;
        }
    }

    if (offset > 0) {
        *queryPtr = word;
    }

    // quality
    isFastq = ks->qual.l > 0;
    for (int i = 0; isFastq && i < readLengths[queryID]; ++i) {
        upkdQualities[queryID * maxReadLength + i] = ks->qual.s[i] - qualityConstant;
    }
    upkdQualities[queryID * maxReadLength + readLengths[queryID]] = 0;

    // name
    if (queryNames[queryID] != NULL) {
        free(queryNames[queryID]);
        queryNames[queryID] = NULL;
    }

    trim_readno(&ks->name);
    queryNames[queryID] = strdup(ks->name.s);

    // comments
    if (queryComments[queryID] != NULL) {
        free(queryComments[queryID]);
        queryComments[queryID] = NULL;
    }

    if (ks->comment.l > 0) queryComments[queryID] = strdup(ks->comment.s);
    // fprintf(stderr, "%d %s %s\n", queryID, queryNames[queryID], ks->seq.s);
}

// This function is to load the pair-end reads for at most "maxNumQueries/2" # of pairs of reads
// Here we store read names and comments without using fixed size arrays
int loadPairReadsKseq(kseq_t *ks1, kseq_t *ks2, unsigned char * charMap, uint * queries, uint * readLengths, uint * readIDs,
                      char * upkdQualities, char ** queryNames, char ** queryComments,
                      uint maxReadLength, uint maxNumQueries, uint wordPerQuery, int qualityConstant, char & isFastq)
{
    int queriesRead = 0;

    while (queriesRead < maxNumQueries) {
        int ksl1 = kseq_read(ks1);
        int ksl2 = kseq_read(ks2);
        if ((ksl1 >= 0 && ksl2 < 0) || (ksl1 < 0 && ksl2 >= 0)) {
            fprintf(stderr, "Error: number of sequences of pair-end files not matched.\n");
            exit(1);
        }

        if (ksl1 < 0) {
            break;
        }

        appendToQueryArrays(ks1, charMap, readLengths, readIDs, queries, upkdQualities, queryNames, queryComments,
                            maxReadLength, wordPerQuery, qualityConstant, queriesRead++, isFastq);
        appendToQueryArrays(ks2, charMap, readLengths, readIDs, queries, upkdQualities, queryNames, queryComments,
                            maxReadLength, wordPerQuery, qualityConstant, queriesRead++, isFastq);
    }

    return queriesRead;
}

// This function is to load the pair-end reads for at most "maxNumQueries/2" # of pairs of reads
int loadPairReadsGz2 ( gzFile queryFile, char * queryFileBuffer, gzFile queryFile2, char * queryFileBuffer2, unsigned char * charMap,
                       uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char ** queryNames, char ** queryComments, uint maxReadLength, uint maxNumQueries,
                       size_t & bufferSize, char & queryChar, uint & bufferIndex, size_t & bufferSize2, char & queryChar2, uint & bufferIndex2,
                       uint accumReadNum, uint wordPerQuery, int qualityConstant, char & isFastq, int maxLenReadName )
{
    // for pair-ended reads
    // load reads from the read file
    // return how many reads are loaded
    ullint queriesRead = 0;
    uint currentWord = 0;
    uint offset = 0;
    uint bits = 0;
    uint * queryPtr = queries;
    int i;
    uint nameCharId;
    int cropped;
    uint nucleoId;

    if ( bufferIndex >= bufferSize )
    {
        bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

        if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
        {
            const char * error_string;
            int err;
            error_string = gzerror ( queryFile, & err );

            if ( err )
            {
                fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                exit ( EXIT_FAILURE );
            }
        }

        bufferIndex = 0;
    }

    // FOR THE FIRST READ
    while ( bufferSize != 0 )
    {
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize != 0 && queryChar != '>' && queryChar != '@' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        isFastq = 0;

        if ( queryChar == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        nameCharId = 0;
        cropped = 0;

        while ( bufferSize != 0  && queryChar != '\n' )
        {
            if ( queryChar == ' ' || queryChar == '\t' ) {cropped = 1;}

            if ( !cropped && nameCharId < maxLenReadName - 1 )
            {
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        //Read the pattern body of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        nucleoId = 0;

        while ( bufferSize != 0  && queryChar != '>' && queryChar != '@' && queryChar != '+' )
        {
            if ( queryChar != '\n' )
            {
                bits = charMap[queryChar];

                if ( nucleoId < maxReadLength )
                {
                    // upkdQueries[queriesRead*maxReadLength+nucleoId]=bits;
                    currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                    offset++;

                    if ( offset == CHAR_PER_WORD )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        if ( offset > 0 )
        { *queryPtr = currentWord; }

        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;

        if ( nucleoId == 0 )
        {
            fprintf ( stderr, "Error! There exists a read with length 0 inside the first read file.\n" );
            exit ( EXIT_FAILURE );
        }

        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Left read #%u is longer than %u! Read truncated.\n", (queriesRead+accumReadNum)/2+1, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize != 0 && queryChar != '+' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize == 0 ) { break; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }

            while ( bufferSize != 0 && queryChar != '\n' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the first file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize == 0 ) { break; }

            for ( i = 0; i < nucleoId; i++ )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( queryChar == '\n' )
                { break; }

                upkdQualities[queriesRead * maxReadLength + i] = queryChar - qualityConstant;

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            if ( i != nucleoId )
            {
                fprintf ( stderr, "Error! Inside the first read file, there exists a read of which the number of qualities does not match with the sequence length.\n" );
                exit ( EXIT_FAILURE );
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        queriesRead += 2;
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

        if ( queriesRead >= maxNumQueries )
        { break; }
    }

    if ( queriesRead == 0 )
    {
        // no read is loaded
        return 0;
    }

    if ( maxNumQueries > queriesRead )
    { maxNumQueries = queriesRead; }

    // FOR THE SECOND READ
    if ( bufferIndex2 >= bufferSize2 )
    {
        bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

        if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
        {
            const char * error_string;
            int err;
            error_string = gzerror ( queryFile2, & err );

            if ( err )
            {
                fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                exit ( EXIT_FAILURE );
            }
        }

        bufferIndex2 = 0;
    }

    queriesRead = 1;
    queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

    while ( bufferSize2 != 0 )
    {
        // FOR THE SECOND READ
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize2 != 0 && queryChar2 != '>' && queryChar2 != '@' )
        {
            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        isFastq = 0;

        if ( queryChar2 == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize2 == 0 ) { break; }

        queryChar2 = queryFileBuffer2[bufferIndex2++];

        if ( bufferIndex2 >= bufferSize2 )
        {
            bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

            if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile2, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex2 = 0;
        }

        nameCharId = 0;
        cropped = 0;

        while ( bufferSize2 != 0  && queryChar2 != '\n' )
        {
            if ( queryChar2 == ' ' || queryChar2 == '\t' ) {cropped = 1;}

            if ( !cropped && nameCharId < maxLenReadName - 1 )
            {
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        // remove "/1" or "/2" if exists

        //Read the pattern body of a read
        if ( bufferSize2 == 0 ) { break; }

        queryChar2 = queryFileBuffer2[bufferIndex2++];

        if ( bufferIndex2 >= bufferSize2 )
        {
            bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

            if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile2, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex2 = 0;
        }

        nucleoId = 0;

        while ( bufferSize2 != 0  && queryChar2 != '>' && queryChar2 != '@' && queryChar2 != '+' )
        {
            if ( queryChar2 != '\n' )
            {
                bits = charMap[queryChar2];

                if ( nucleoId < maxReadLength )
                {
                    // upkdQueries[queriesRead*maxReadLength+nucleoId]=bits;
                    currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                    offset++;

                    if ( offset == CHAR_PER_WORD )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        if ( offset > 0 )
        { *queryPtr = currentWord; }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
        // reverse the second read
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );
        uint readLength = nucleoId;
        uint leftWordIndex = 0;
        uint rightWordIndex = ( readLength - 1 ) / CHAR_PER_WORD;
        uint leftWord = queryPtr[leftWordIndex * 32];
        uint rightWord = queryPtr[rightWordIndex * 32];

        for ( i = 0, j = readLength - 1; i <= j; ++i, --j )
        {
            // check if need to move to next word
            if ( i / CHAR_PER_WORD != leftWordIndex )
            {
                // write back leftword
                queryPtr[leftWordIndex * 32] = leftWord;
                // load next leftword
                leftWordIndex++;
                leftWord = queryPtr[leftWordIndex * 32];
            }

            if ( j / CHAR_PER_WORD != rightWordIndex )
            {
                // write back rightword
                queryPtr[rightWordIndex * 32] = rightWord;
                // load next rightword
                rightWordIndex--;
                rightWord = queryPtr[rightWordIndex * 32];
            }

            // swap left and right characters
            unsigned char leftChar = ( leftWord >> ( i % CHAR_PER_WORD * 2 ) ) & CHAR_MASK;
            unsigned char rightChar = ( rightWord >> ( j % CHAR_PER_WORD * 2 ) ) & CHAR_MASK;
            leftWord ^= ( leftChar ^ ( soap3DnaComplement[rightChar] ) ) << ( i % CHAR_PER_WORD * 2 );
            rightWord ^= ( ( soap3DnaComplement[leftChar] ) ^ rightChar ) << ( j % CHAR_PER_WORD * 2 );
        }

        // write back
        if ( leftWordIndex == rightWordIndex )
        {
            uint numLeftBits = ( ( readLength - 1 ) / 2 % CHAR_PER_WORD + 1 ) * BIT_PER_CHAR;
            uint numRightBits = 32 - numLeftBits;
            queryPtr[leftWordIndex * 32] = ( ( leftWord << numRightBits ) >> numRightBits ) |
                                           ( ( rightWord >> numLeftBits ) << numLeftBits );
        }
        else
        {
            queryPtr[leftWordIndex * 32] = leftWord;
            queryPtr[rightWordIndex * 32] = rightWord;
        }

#endif
        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;

        if ( nucleoId == 0 )
        {
            fprintf ( stderr, "Error! There exists a read with length 0 inside the second read file.\n" );
            exit ( EXIT_FAILURE );
        }

        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Right read #%u is longer than %u! Read truncated.\n", (queriesRead+accumReadNum)/2+1, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize2 != 0 && queryChar2 != '+' )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize2 == 0 ) { break; }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }

            while ( bufferSize2 != 0 && queryChar2 != '\n' )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize2 == 0 ) { break; }

            for ( i = 0; i < nucleoId; i++ )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( queryChar2 == '\n' )
                { break; }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
                // reverse the quality of the read
                upkdQualities[queriesRead * maxReadLength + nucleoId - 1 - i] = queryChar2 - qualityConstant;
#else
                upkdQualities[queriesRead * maxReadLength + i] = queryChar2 - qualityConstant;
#endif

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            if ( i != nucleoId )
            {
                fprintf ( stderr, "Error! Inside the second read file, there exists a read of which the number of qualities does not match with the sequence length.\n" );
                exit ( EXIT_FAILURE );
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        queriesRead += 2;
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

        if ( queriesRead >= maxNumQueries )
        { break; }
    }

    queriesRead--; // number of pairs of reads loaded
    return queriesRead;
}

// This function is used to load the reads with random skipping scheme
// for Score Recalibration in SNP-calling module.
// Initially, first "pairsToSkip" # of read pairs will be skipped.
// When "pairsToSkip" becomes 0, the read pairs will be included,
// a random number [0, 3] will be chosen as the next "pairsToSkip".
int loadPairReadsGzForScoreRecalibration ( gzFile queryFile, char * queryFileBuffer, gzFile queryFile2, char * queryFileBuffer2,
                       uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char ** queryNames, char ** queryComments, uint maxReadLength, uint maxNumQueries,
                       size_t & bufferSize, char & queryChar, uint & bufferIndex, size_t & bufferSize2, char & queryChar2, uint & bufferIndex2,
                       uint accumReadNum, uint wordPerQuery, int qualityConstant, char & isFastq, int maxLenReadName,
                       uint & pairsToSkip )
{
    // for pair-ended reads
    // load reads from the read file
    // return how many reads are loaded
    // Convert ACGT to 0123
    uint charMap[256];

    for ( int i = 0; i < 256; i++ )
        { charMap[i] = 0; }

    charMap['A'] = 0;
    charMap['C'] = 1;
    charMap['G'] = 2;
    charMap['T'] = 3;
    charMap['U'] = 3;
    charMap['N'] = 2; // N -> G
    charMap['a'] = 0;
    charMap['c'] = 1;
    charMap['g'] = 2;
    charMap['t'] = 3;
    charMap['u'] = 3;
    charMap['n'] = 2; // N -> G
    ullint queriesRead = 0;
    uint currentWord = 0;
    uint offset = 0;
    uint bits = 0;
    uint * queryPtr = queries;
    int i;
    uint nameCharId;
    int cropped;
    uint nucleoId;

    unsigned int toSkip = pairsToSkip;

    if ( bufferIndex >= bufferSize )
    {
        bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

        if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
        {
            const char * error_string;
            int err;
            error_string = gzerror ( queryFile, & err );

            if ( err )
            {
                fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                exit ( EXIT_FAILURE );
            }
        }

        bufferIndex = 0;
    }
    
    if ( bufferIndex2 >= bufferSize2 )
    {
        bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

        if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
        {
            const char * error_string;
            int err;
            error_string = gzerror ( queryFile2, & err );

            if ( err )
            {
                fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                exit ( EXIT_FAILURE );
            }
        }

        bufferIndex2 = 0;
    }

    // FOR THE FIRST READ
    while ( bufferSize != 0 && bufferSize2 != 0 )
    {
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize != 0 && queryChar != '>' && queryChar != '@' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        isFastq = 0;

        if ( queryChar == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        nameCharId = 0;
        cropped = 0;

        while ( bufferSize != 0  && queryChar != '\n' )
        {
            if ( queryChar == ' ' || queryChar == '\t' ) {cropped = 1;}

            if ( !cropped && nameCharId < maxLenReadName - 1 )
            {
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        // remove "/1" or "/2" if exists

        //Read the pattern body of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        nucleoId = 0;

        while ( bufferSize != 0  && queryChar != '>' && queryChar != '@' && queryChar != '+' )
        {
            if ( queryChar != '\n' )
            {
                bits = charMap[queryChar];

                if ( nucleoId < maxReadLength )
                {
                    // upkdQueries[queriesRead*maxReadLength+nucleoId]=bits;
                    currentWord |= ( bits << ( offset * 2 ) );
                    offset++;

                    if ( offset == 16 )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        if ( offset > 0 )
            { *queryPtr = currentWord; }

        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;

        if ( nucleoId == 0 )
        {
            fprintf ( stderr, "Error! There exists a read with length 0 inside the first read file.\n" );
            exit ( EXIT_FAILURE );
        }

        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Left read #%u is longer than %u! Read truncated.\n", (queriesRead+accumReadNum)/2+1, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize != 0 && queryChar != '+' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize == 0 ) { break; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }

            while ( bufferSize != 0 && queryChar != '\n' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the first file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize == 0 ) { break; }

            for ( i = 0; i < nucleoId; i++ )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( queryChar == '\n' )
                    { break; }

                upkdQualities[queriesRead * maxReadLength + i] = queryChar - qualityConstant;

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            if ( i != nucleoId )
            {
                fprintf ( stderr, "Error! Inside the first read file, there exists a read of which the number of qualities does not match with the sequence length.\n" );
                exit ( EXIT_FAILURE );
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        if ( toSkip == 0 )
        {
            queriesRead++;
        }

        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );
        
        while ( bufferSize2 != 0 && queryChar2 != '>' && queryChar2 != '@' )
        {
            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        isFastq = 0;

        if ( queryChar2 == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize2 == 0 ) { break; }

        queryChar2 = queryFileBuffer2[bufferIndex2++];

        if ( bufferIndex2 >= bufferSize2 )
        {
            bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

            if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile2, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex2 = 0;
        }

        nameCharId = 0;
        cropped = 0;

        while ( bufferSize2 != 0  && queryChar2 != '\n' )
        {
            if ( queryChar2 == ' ' || queryChar2 == '\t' ) {cropped = 1;}

            if ( !cropped && nameCharId < maxLenReadName - 1 )
            {
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        // remove "/1" or "/2" if exists

        //Read the pattern body of a read
        if ( bufferSize2 == 0 ) { break; }

        queryChar2 = queryFileBuffer2[bufferIndex2++];

        if ( bufferIndex2 >= bufferSize2 )
        {
            bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

            if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile2, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex2 = 0;
        }

        nucleoId = 0;

        while ( bufferSize2 != 0  && queryChar2 != '>' && queryChar2 != '@' && queryChar2 != '+' )
        {
            if ( queryChar2 != '\n' )
            {
                bits = charMap[queryChar2];

                if ( nucleoId < maxReadLength )
                {
                    // upkdQueries[queriesRead*maxReadLength+nucleoId]=bits;
                    currentWord |= ( bits << ( offset * 2 ) );
                    offset++;

                    if ( offset == 16 )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        if ( offset > 0 )
            { *queryPtr = currentWord; }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
        // reverse the second read
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );
        uint readLength = nucleoId;
        uint leftWordIndex = 0;
        uint rightWordIndex = ( readLength - 1 ) / CHAR_PER_WORD;
        uint leftWord = queryPtr[leftWordIndex * 32];
        uint rightWord = queryPtr[rightWordIndex * 32];

        for ( i = 0, j = readLength - 1; i <= j; ++i, --j )
        {
            // check if need to move to next word
            if ( i / CHAR_PER_WORD != leftWordIndex )
            {
                // write back leftword
                queryPtr[leftWordIndex * 32] = leftWord;
                // load next leftword
                leftWordIndex++;
                leftWord = queryPtr[leftWordIndex * 32];
            }

            if ( j / CHAR_PER_WORD != rightWordIndex )
            {
                // write back rightword
                queryPtr[rightWordIndex * 32] = rightWord;
                // load next rightword
                rightWordIndex--;
                rightWord = queryPtr[rightWordIndex * 32];
            }

            // swap left and right characters
            unsigned char leftChar = ( leftWord >> ( i % CHAR_PER_WORD * 2 ) ) & 3;
            unsigned char rightChar = ( rightWord >> ( j % CHAR_PER_WORD * 2 ) ) & 3;
            leftWord ^= ( leftChar ^ ( 3 - rightChar ) ) << ( i % CHAR_PER_WORD * 2 );
            rightWord ^= ( ( 3 - leftChar ) ^ rightChar ) << ( j % CHAR_PER_WORD * 2 );
        }

        // write back
        if ( leftWordIndex == rightWordIndex )
        {
            uint numLeftBits = ( ( readLength - 1 ) / 2 % CHAR_PER_WORD + 1 ) * 2;
            uint numRightBits = 32 - numLeftBits;
            queryPtr[leftWordIndex * 32] = ( ( leftWord << numRightBits ) >> numRightBits ) |
                                           ( ( rightWord >> numLeftBits ) << numLeftBits );
        }
        else
        {
            queryPtr[leftWordIndex * 32] = leftWord;
            queryPtr[rightWordIndex * 32] = rightWord;
        }

#endif
        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;

        if ( nucleoId == 0 )
        {
            fprintf ( stderr, "Error! There exists a read with length 0 inside the second read file.\n" );
            exit ( EXIT_FAILURE );
        }

        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Right read #%u is longer than %u! Read truncated.\n", (queriesRead+accumReadNum)/2+1, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize2 != 0 && queryChar2 != '+' )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize2 == 0 ) { break; }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }

            while ( bufferSize2 != 0 && queryChar2 != '\n' )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize2 == 0 ) { break; }

            for ( i = 0; i < nucleoId; i++ )
            {
                queryChar2 = queryFileBuffer2[bufferIndex2++];

                if ( queryChar2 == '\n' )
                    { break; }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
                // reverse the quality of the read
                upkdQualities[queriesRead * maxReadLength + nucleoId - 1 - i] = queryChar2 - qualityConstant;
#else
                upkdQualities[queriesRead * maxReadLength + i] = queryChar2 - qualityConstant;
#endif

                if ( bufferIndex2 >= bufferSize2 )
                {
                    bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                    if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile2, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex2 = 0;
                }
            }

            if ( i != nucleoId )
            {
                fprintf ( stderr, "Error! Inside the second read file, there exists a read of which the number of qualities does not match with the sequence length.\n" );
                exit ( EXIT_FAILURE );
            }

            queryChar2 = queryFileBuffer2[bufferIndex2++];

            if ( bufferIndex2 >= bufferSize2 )
            {
                bufferSize2 = gzread ( queryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

                if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile2 ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile2, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex2 = 0;
            }
        }

        if ( toSkip == 0 )
        {
            queriesRead++;
            toSkip = rand() % 4;
        }
        else
        {
            toSkip--;
        }

        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

        pairsToSkip = toSkip;

        if ( queriesRead >= maxNumQueries )
            { break; }
    }

    return queriesRead;
}

// This function is to load the single reads for at most "maxNumQueries" # of reads
int loadSingleReadsGz ( gzFile queryFile, char * queryFileBuffer, unsigned char * charMap,
                        uint * queries, uint * readLengths, uint * readIDs,
                        char * upkdQualities, char ** queryNames, char ** queryComments, uint maxReadLength, uint maxNumQueries,
                        size_t & bufferSize, char & queryChar, uint & bufferIndex, uint accumReadNum, uint wordPerQuery,
                        int qualityConstant, char & isFastq, int maxLenReadName )
{
    // for single reads
    // load reads from the read file
    // return how many reads are loaded
    ullint queriesRead = 0;
    uint currentWord = 0;
    uint offset = 0;
    uint bits = 0;
    uint * queryPtr = queries;
    int i;

    while ( bufferSize != 0 )
    {
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize != 0 && queryChar != '>' && queryChar != '@' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        isFastq = 0;

        if ( queryChar == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        uint nameCharId = 0;
        int cropped = 0;

        while ( bufferSize != 0  && queryChar != '\n' )
        {
            if ( queryChar == ' ' || queryChar == '\t' ) {cropped = 1;}

            if ( !cropped && nameCharId < maxLenReadName - 1 )
            {
                // TODO
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        // remove "/1" or "/2" if exists
        // TODO

        //Read the pattern body of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

            if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( queryFile, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex = 0;
        }

        uint nucleoId = 0;

        while ( bufferSize != 0  && queryChar != '>' && queryChar != '@' && queryChar != '+' )
        {
            if ( queryChar != '\n' )
            {
                bits = charMap[queryChar];

                if ( nucleoId < maxReadLength )
                {
                    // upkdQueries[queriesRead*maxReadLength+nucleoId]=bits;
                    currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                    offset++;

                    if ( offset == CHAR_PER_WORD )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        if ( offset > 0 )
        { *queryPtr = currentWord; }

        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;

        if ( nucleoId == 0 )
        {
            fprintf ( stderr, "Error! There exists a read with length 0 inside the read file.\n" );
            exit ( EXIT_FAILURE );
        }

        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Read #%u is longer than %u! Read truncated.\n", queriesRead+1+accumReadNum, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize != 0 && queryChar != '+' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize == 0 ) { break; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }

            while ( bufferSize != 0 && queryChar != '\n' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize == 0 ) { break; }

            for ( i = 0; i < nucleoId; i++ )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( queryChar == '\n' )
                { break; }

                upkdQualities[queriesRead * maxReadLength + i] = queryChar - qualityConstant;

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                    if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                    {
                        const char * error_string;
                        int err;
                        error_string = gzerror ( queryFile, & err );

                        if ( err )
                        {
                            fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                            exit ( EXIT_FAILURE );
                        }
                    }

                    bufferIndex = 0;
                }
            }

            if ( i != nucleoId )
            {
                fprintf ( stderr, "Error! Inside the read file, there exists a read of which the number of qualities does not match with the sequence length.\n" );
                exit ( EXIT_FAILURE );
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = gzread ( queryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

                if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( queryFile ) ) )
                {
                    const char * error_string;
                    int err;
                    error_string = gzerror ( queryFile, & err );

                    if ( err )
                    {
                        fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                        exit ( EXIT_FAILURE );
                    }
                }

                bufferIndex = 0;
            }
        }

        queriesRead++;
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

        if ( queriesRead >= maxNumQueries )
        { break; }
    }

    return queriesRead;
}


uint GetReadLength ( uint * readLengths, uint numQueries, int sample )
{
    if ( sample <= 0 ) {
        sample = 1;
    }
        
    // scan the first ONE MILLION reads
    // and get the max read length among them
    if ( numQueries <= 0 )
    { return 100; }

    int i = 0;
    int j = 1;
    uint maxReadLength = readLengths[i];
    while (i<numQueries && j<1000000)
    {
        if ( maxReadLength < readLengths[i] )
        { maxReadLength = readLengths[i]; }
        
        j++;
        i+=sample;
    }

    return maxReadLength;
}


/*
 *
 *    dependencies.h
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

#ifndef __DEPENDENCIES_H__
#define __DEPENDENCIES_H__

// #define BGS_OUTPUT_KERNEL_DEBUG_MESSAGE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "2bwt-flex/OCC.h"
#include "2bwt-flex/SRAArguments.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"
#include "definitions.h"

static void getSeedPositions ( int stage, int readLength, int * seedLength, int * seedPositions, int * seedNum )
{
    if ( stage == STAGE_SINGLE_DP ||
            stage == STAGE_NEW_DEFAULT_DP )
    {
        if ( readLength > 300 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_V_LONG_READ; }
        else if ( readLength > 80 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_LONG_READ; }
        else if ( readLength > 60 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_MEDIAN_READ; }
        else if ( readLength > 40 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_SHORT_READ; }
        else
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_V_SHORT_READ; }

        // update: for reads longer than 120, one more seed for every extra 100 bases
        if ( readLength > 120 )
        {
            ( *seedNum ) = SEED_NUM_SINGLE_DP + readLength / 100;
        }
        else
        {
            ( *seedNum ) = SEED_NUM_SINGLE_DP;
        }

        int H = 0;

        if ( readLength > 300 )
        { H = readLength * 0.15; }

        int X;

        if ( readLength > 300 )
        { X = readLength * 0.15; }
        else if ( readLength > 80 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_LONG_READ; }
        else if ( readLength > 60 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_MEDIAN_READ; }
        else if ( readLength > 40 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_SHORT_READ; }
        else
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_V_SHORT_READ; }

        int seed_apart = ( readLength - X - H ) / ( *seedNum );
        int i;

        for ( i = 0; i < ( *seedNum ); i++ )
        {
            seedPositions[i] = H + i * seed_apart;
        }

        if ( seedPositions[ ( *seedNum ) - 1] > readLength - ( *seedLength ) - X )
        { seedPositions[ ( *seedNum ) - 1] = readLength - ( *seedLength ) - X; }
    }
    else if ( stage == STAGE_DEEP_DP_ROUND1 )
    {
        (*seedNum) = NUMBER_OF_SEED;
        (*seedLength) = HALF_SEED_SIZE_FOR_V_LONG_READ_UB * 2;
        seedPositions[0] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_1;
        seedPositions[1] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_2;
        seedPositions[2] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_3;
        seedPositions[3] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_4;
        seedPositions[4] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_5;
        seedPositions[5] = SEED_POS_DEEP_DP_FOR_V_LONG_READ_6;
       // seedPositions[1]=40;
    }
    else if ( stage == STAGE_DEEP_DP_ROUND2A )
    {
        /*
        int seed_size_lb, seed_size_ub;
        if ( readLength > 150 )
        {
             seed_size_lb = HALF_SEED_SIZE_FOR_V_LONG_READ_LB;
             seed_size_ub = HALF_SEED_SIZE_FOR_V_LONG_READ_UB;
             ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_LONG_READ;
        }
        else if ( readLength > 80 )
        {
             seed_size_lb = HALF_SEED_SIZE_FOR_LONG_READ_LB;
             seed_size_ub = HALF_SEED_SIZE_FOR_LONG_READ_UB;
             ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_LONG_READ;
        }

        int idx = static_cast <int> (floor( (readLength-1) / seed_size_lb )) - 1;
        int * pos = (int *)malloc( idx * sizeof( idx ) );
        pos[idx-1] = readLength - ( *seedLength ) + ( seed_size_ub - seed_size_lb );
        int i;
        for ( i=idx-1 ; i>0 ;--i)
            pos[i-1] = pos[i]-seed_size_lb;
		
        ( *seedNum ) = 0;
        // seedPositions[ (*seedNum)++ ] = pos[1];
        for ( i=0; i<idx; i+=2, ++ (*seedNum) )
            seedPositions[( *seedNum ) ] = pos[i];
        
        if ( seedPositions[0] < 0 )
            seedPositions[0] = 0;
        */
        /*
        if ( readLength > 150 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_LONG_READ_2; }
        else if ( readLength > 80 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_LONG_READ_2; }
        else if ( readLength > 60 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ_2; }
        else if ( readLength > 40 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_SHORT_READ_2; }
        else
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_SHORT_READ_2; }

        int H = 0;
        int T = 0;

        if ( readLength > 150 )
        {
            H = readLength * 0.1;
            T = readLength * 0.2;
        }

        ( *seedNum ) = 0;

        for ( int i = readLength - ( *seedLength ) - T; i >= H; i -= ( *seedLength ) / 2 )
        {
            seedPositions[ ( *seedNum ) ++] = i;
        }

        if ( seedPositions[ ( *seedNum ) - 1] > H )
        { seedPositions[ ( *seedNum ) ++] = H; }
        */
    }
    else if ( stage == STAGE_DEEP_DP_ROUND2B )
    {
        int seed_size_lb, seed_size_ub;
        if ( readLength > 150 )
        {
             seed_size_lb = HALF_SEED_SIZE_FOR_V_LONG_READ_LB;
             seed_size_ub = HALF_SEED_SIZE_FOR_V_LONG_READ_UB;
             ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_LONG_READ;
        }
        else if ( readLength > 80 )
        {
             seed_size_lb = HALF_SEED_SIZE_FOR_LONG_READ_LB;
             seed_size_ub = HALF_SEED_SIZE_FOR_LONG_READ_UB;
             ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_LONG_READ;
        }

        int idx = static_cast <int> (floor( (readLength-1) / seed_size_lb )) - 1;
        int * pos = (int *)malloc( idx * sizeof( idx ) );
        pos[idx-1] = readLength - ( seed_size_lb * 2 ) - ( seed_size_ub - seed_size_lb );
        int i;
        for ( i=idx-1 ; i>0 ;--i)
            pos[i-1] = pos[i]-seed_size_lb;
		
        ( *seedNum ) = 0;
        // seedPositions[ (*seedNum)++ ] = pos[1];
        for ( i=0; i<idx; i+=2, ++ (*seedNum) )
            seedPositions[( *seedNum ) ] = pos[i];
        
        if ( seedPositions[0] < 0 )
            seedPositions[0] = 0;
    }

}

static unsigned int getWordPerQuery ( unsigned int maxReadLength )
{
    return (maxReadLength + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
}

#endif

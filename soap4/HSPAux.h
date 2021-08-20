/*
 *
 *    HSPAux.h
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
 
#ifndef __HSP_AUXILLIARY_H__
#define __HSP_AUXILLIARY_H__

#define MAX_NUMBER_OF_OVERFLOW_FILE_PTR 16 // should be same as MAX_NUM_CPU_THREADS

typedef struct HSPAux
{
	char isFastq; // input file format - 0: No; 1: Yes
	int dpMatchScore; // match sccore for DP
	int dpMisMatchScore; // mismatch score forDP
	int alignmentType; // alignment type
	int readType; // read type: single end or paired end
	int peMaxOutputPerRead; // maximum number of hits reported per read for paired-end alignment

//add these back for calculating MAPQ
    int * x0_array; // the x0 array of the soap3 result
    int * x1_array; // the x1 array of the soap3 result
    int * mismatch_array; // the min number of mismatches of the soap3 result
//

	int minMAPQ; // minimum value of MAPQ
	int maxMAPQ; // maximum value of MAPQ
	int bwaLikeScore; // whether bwa-like MAPQ score is reported
	int g_log_n[256]; // store the value of 4.343 * log(i), for i = 0...255;
	char * readGroup; // read group
	char * sampleName; // sample name
	char * readGrpOption; // read group option
	int maxLenReadName; // max length for read name
	int ProceedDPForTooManyHits; // whether the seed will proceed to perform DP if there are too many hits
	double singleDPcutoffRatio;
	double singleDPcutoffLB;
    int isPrintMDNM;

	// the following arrays for storing the unaligned read IDs
	// for proceeding single deep-dp on them
	void * readsIDForSingleDP; // type: UnalignedSinglesArrays
	void * allHits; // store the corresponding alignment result

    int dpStageId;

    int megapathMode;
    int outputBAM;
    int quality_constant;
    double top_percentage;
} HSPAux;

#endif

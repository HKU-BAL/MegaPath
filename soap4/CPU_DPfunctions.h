#ifndef _CPU_DP_FUNCTIONS_H
#define _CPU_DP_FUNCTIONS_H

#define USE_CPU_DP
#ifdef USE_CPU_DP

#include "AlgnResult.h"
#include "alignment.h"

#ifdef __AVX2__
#include <immintrin.h>
#else
#include <xmmintrin.h>
#endif

#include <sys/time.h>
#include <pthread.h>
#include <semaphore.h>

#include <string>
#include <stack>
#include <map>
#include <vector>
using namespace std;

#define DP_THREADS_PER_BLOCK 128

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

//////////////////////////////////////////////////////////////////////////////////////
// Functions for the DP module //
//////////////////////////////////////////////////////////////////////////////////////

#define MC_MemberCopy(type,prefix,x) type x = prefix x
#define MC_MemberCopy2(type,prefix,x,y) MC_MemberCopy(type,prefix,x); MC_MemberCopy(type,prefix,y)
#define MC_MemberCopy3(type,prefix,x,y,z) MC_MemberCopy2(type,prefix,x,y); MC_MemberCopy(type,prefix,z)
#define MC_MemberCopy4(type,prefix,x,y,z,u) MC_MemberCopy3(type,prefix,x,y,z); MC_MemberCopy(type,prefix,u)
#define MC_MemberCopy5(type,prefix,x,y,z,u,v) MC_MemberCopy4(type,prefix,x,y,z,u); MC_MemberCopy(type,prefix,v)

#define MC_Max(x,y) (x > y ? x : y)
#define MC_CeilDivide16(x) ((x+15)>>4)

static void * addr;
#define MC_CheckMalloc(var,type,size) \
	addr = malloc((size) * sizeof(type)); \
	if (addr == NULL) { \
		fprintf(stderr, "[DPfunc] error: main memory allocation failed\n"); \
        fprintf ( stderr, "error at %s: %u\n", __FILE__, __LINE__ ); \
        fprintf ( stderr, "hihi %u %u\n", (size),sizeof (type) ); \
		exit(-1); \
	} \
	var = (type *) addr; \

struct DPThreadsWrapperObj;

class SemiGlobalAligner
{
		int n_conf, blockConf[16];
		double coefConf[16];

		int batchSize, maxReadLength, maxDNALength, maxDPTableLength;
		DPParameters dpPara;
		int alignmentScheme;

		void * _DPTable;
#ifdef __AVX2__
		__m256i* _DPTempRow;
#else
		__m128i* _DPTempRow;
#endif
		uint * _packedDNASequence, *_DNALengths;
		uint * _packedReadSequence, *_readLengths;
		uint * _startLocs, *_startOffsets, *_hitLocs;
		int * _scores, *_cutoffThresholds;
		uint * _clipLtSizes, *_clipRtSizes;
		uint * _anchorLeftLocs, *_anchorRightLocs;
		uchar * _pattern;
		uint * _maxScoreCounts;

        
		uint estimateThreadSize ( int maxReadLength, int maxDNALength );
		int tryAlloc ( size_t estimatedThreadSize, size_t numOfBlocks );

	public:
		SemiGlobalAligner();
		void decideConfiguration (
		    int maxReadLength, int maxDNALength,
		    int & maxDPTableLength, int & numOfBlocks,
		    int & patternLength, DPParameters & dpPara
		);
        void decideConfigurationForIndelRA (
		    int maxReadLength, int maxDNALength,
		    int & maxDPTableLength, int & numOfBlocks,
		    int & patternLength, DPParameters & dpPara
		);
		void init (
		    int batchSize,
		    int maxReadLength, int maxDNALength, int maxDPTableLength,
		    DPParameters & dpPara
		);
		void performAlignment (
		    uint * packedDNASequence, uint * DNALengths,
		    uint * packedReadSequence, uint * readLengths,
		    int * cutoffThresholds, int * scores, uint * hitLocs,
		    uint * maxScoreCounts,
		    uchar * pattern, int numOfThreads,
		    uint * clipLtSizes = NULL, uint * clipRtSizes = NULL,
		    uint * anchorLeftLocs = NULL, uint * anchorRightLocs = NULL
		);

		void freeMemory();

};

#endif

#endif

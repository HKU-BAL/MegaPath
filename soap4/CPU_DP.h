#include <cstdio>
#include <cstring>

#ifdef __AVX2__
#include <immintrin.h>
#else
#include <smmintrin.h>
#endif

void callDP ( uint * DNASequences, uint * DNALengths, uint maxDNALength,
              uint * readSequences, uint * readLengths, uint maxReadLength,
              int clipLtSizes, int clipRtSizes,
              uint * anchorLeftLocs, uint * anchorRightLocs,
              uint MatchScore, uint MismatchScore,
              uint GapOpenScore, uint GapExtendScore,
              int * cutOffThresholds,
#ifdef __AVX2__
              void * DPTable, __m256i* DPTempRow, uint numDPInstances,
#else
              void * DPTable, __m128i* DPTempRow, uint numDPInstances,
#endif
              // output
              int * scores, uint * hitLocs,
              uint * maxScoreCounts, uchar * pattern );

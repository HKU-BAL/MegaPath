#define INVALID_CIGAR_STRING_MARK 12345
#include <assert.h>
#include "definitions.h"
#include <stdint.h>
#include <cstdio>
#include <cstring>
#include <algorithm>
#define MC_CeilDivide16(x) ((x+15)>>4)
#define MC_DnaUnpack(X,i) ((X[dnaTPARA + (((i)>>4)<<5)] >> ((15-((i)&0xF))<<1)) & 3)
#define MC_ReadUnpack(X,i) ((readSequences[readTPARA + (((i)>>4)<<5)] >> ((15-((i)&0xF))<<1)) & 3)

#define MC_ScoreAddr(X,j,i) *(X + ( ((j)) * ( maxReadLength + 1 ) + (i) ) * 2)
#define MC_getOffset(j,i) ( ( ((j)) * ( maxReadLength + 1 ) + (i) ) * 2 )

#define DP_SCORE_NEG_INFINITY -32000
typedef unsigned int uint;
typedef unsigned char uchar;



#ifdef __AVX2__
// see manual at https://software.intel.com/sites/landingpage/IntrinsicsGuide/
#include <immintrin.h>
#define MC_SetEpi8(Y,base) _mm256_set_epi8( Y[(base)], Y[(base)+1], Y[(base)+2], Y[(base)+3], Y[(base)+4], Y[(base)+5], Y[(base)+6], Y[(base)+7], Y[(base)+8], Y[(base)+9], Y[(base)+10], Y[(base)+11], Y[(base)+12], Y[(base)+13], Y[(base)+14], Y[(base)+15], Y[(base)+16], Y[(base)+17], Y[(base)+18], Y[(base)+19], Y[(base)+20], Y[(base)+21], Y[(base)+22], Y[(base)+23], Y[(base)+24], Y[(base)+25], Y[(base)+26], Y[(base)+27], Y[(base)+28], Y[(base)+29], Y[(base)+30], Y[(base)+31] )
#define N_DP_PACKED 32
typedef __m256i SIMDtype;
#define SIMD_SET1(x) _mm256_set1_epi8((x))
#define SIMD_AND(x,y) _mm256_and_si256((x),(y))
#define SIMD_OR(x,y) _mm256_or_si256((x),(y))
#define SIMD_XOR(x,y) _mm256_xor_si256((x),(y))
#define SIMD_NAND(x,y) _mm256_andnot_si256((x),(y))
#define SIMD_SRLI(x,y) _mm256_srli_epi16((x),(y))
#define SIMD_SLLI(x,y) _mm256_slli_epi16((x),(y))
#define SIMD_ADDS_EPI(x,y) _mm256_adds_epi8((x),(y))
#define SIMD_SUB_EPI(x,y) _mm256_sub_epi8((x),(y))
#define SIMD_CMPEQ(x,y) _mm256_cmpeq_epi8((x),(y))
#define SIMD_ADDS_EPU(x,y) _mm256_adds_epu8((x),(y))
#define SIMD_ADDS16_EPU(x,y) _mm256_adds_epu16((x),(y))
#define SIMD_SUBS_EPU(x,y) _mm256_subs_epu8((x),(y))
#define SIMD_SUBS16_EPU(x,y) _mm256_subs_epu16((x),(y))
#define SIMD_STORE(x,y) _mm256_store_si256((x),(y))
#define SIMD_LOAD(x) _mm256_load_si256((x))
#define SIMD_MIN_EPU(x,y) _mm256_min_epu8((x),(y))
#define SIMD_MAX_EPU(x,y) _mm256_max_epu8((x),(y))
#define SIMD_MIN16_EPU(x,y) _mm256_min_epu16((x),(y))
#define SIMD_MAX16_EPU(x,y) _mm256_max_epu16((x),(y))
#define SIMD_BLENDV(x,y,z) _mm256_blendv_epi8((x),(y),(z)) // z?x:y
#define SIMD_SET_SHUFFLE_MASK(x,y,z,w) _mm256_set_epi32((x),(y),(z),(w),(x),(y),(z),(w)) 
#define SIMD_MAX128(x) _mm_max_epu16(_mm256_extracti128_si256((x), 1), _mm256_extracti128_si256((x), 0))
#define SIMD_MIN128(x) _mm_min_epu16(_mm256_extracti128_si256((x), 1), _mm256_extracti128_si256((x), 0))
#define SIMD_TESTZ(x,y) _mm256_testz_si256((x),(y))
#define SIMD_UNPACKLO(x,y) _mm256_unpacklo_epi8((x),(y))
#define SIMD_UNPACKHI(x,y) _mm256_unpackhi_epi8((x),(y))
#define SIMD_SHUFFLE(x,y) _mm256_shuffle_epi8((x),(y))


#else

#include <smmintrin.h>
#define MC_SetEpi8(Y,base) _mm_set_epi8( Y[(base)], Y[(base)+1], Y[(base)+2], Y[(base)+3], Y[(base)+4], Y[(base)+5], Y[(base)+6], Y[(base)+7], Y[(base)+8], Y[(base)+9], Y[(base)+10], Y[(base)+11], Y[(base)+12], Y[(base)+13], Y[(base)+14], Y[(base)+15] )
#define N_DP_PACKED 16
typedef __m128i SIMDtype;
#define SIMD_SET1(x) _mm_set1_epi8((x))
#define SIMD_AND(x,y) _mm_and_si128((x),(y))
#define SIMD_OR(x,y) _mm_or_si128((x),(y))
#define SIMD_XOR(x,y) _mm_xor_si128((x),(y))
#define SIMD_NAND(x,y) _mm_andnot_si128((x),(y))
#define SIMD_SRLI(x,y) _mm_srli_epi16((x),(y))
#define SIMD_SLLI(x,y) _mm_slli_epi16((x),(y))
#define SIMD_ADDS_EPI(x,y) _mm_adds_epi8((x),(y))
#define SIMD_SUB_EPI(x,y) _mm_sub_epi8((x),(y))
#define SIMD_CMPEQ(x,y) _mm_cmpeq_epi8((x),(y))
#define SIMD_ADDS_EPU(x,y) _mm_adds_epu8((x),(y))
#define SIMD_ADDS16_EPU(x,y) _mm_adds_epu16((x),(y))
#define SIMD_SUBS_EPU(x,y) _mm_subs_epu8((x),(y))
#define SIMD_SUBS16_EPU(x,y) _mm_subs_epu16((x),(y))
#define SIMD_STORE(x,y) _mm_store_si128((x),(y))
#define SIMD_LOAD(x) _mm_load_si128((x))
#define SIMD_MIN_EPU(x,y) _mm_min_epu8((x),(y))
#define SIMD_MAX_EPU(x,y) _mm_max_epu8((x),(y))
#define SIMD_MIN16_EPU(x,y) _mm_min_epu16((x),(y))
#define SIMD_MAX16_EPU(x,y) _mm_max_epu16((x),(y))
#define SIMD_BLENDV(x,y,z) _mm_blendv_epi8((x),(y),(z))
#define SIMD_SET_SHUFFLE_MASK(x,y,z,w) _mm_set_epi32((x),(y),(z),(w))
#define SIMD_MAX128(x) (x)
#define SIMD_MIN128(x) (x)
#define SIMD_TESTZ(x,y) _mm_testz_si128((x),(y))
#define SIMD_UNPACKLO(x,y) _mm_unpacklo_epi8((x),(y))
#define SIMD_UNPACKHI(x,y) _mm_unpackhi_epi8((x),(y))
#define SIMD_SHUFFLE(x,y) _mm_shuffle_epi8((x),(y))


#endif

inline uint8_t extract_epi8 (const SIMDtype &vec, const short idx) {
    return ((uint8_t*)&vec)[N_DP_PACKED -1 -idx];
}

#define PACK_DPTEMP(MhOrMv, DmOrIm) SIMD_OR ( MhOrMv, SIMD_SLLI ( DmOrIm, 4 ) )
#define UNPACK_DPTEMP_MhOrMv(DP) SIMD_AND ( DP, BITS_MASK_4_OF_8 )
#define UNPACK_DPTEMP_DmOrIm(DP) SIMD_AND ( SIMD_SRLI ( DP, 4 ) , BITS_MASK_4_OF_8 )

//#define PACK_DPTABLE(MhOrMv, Md, DmOrImBit) SIMD_OR ( Md, SIMD_OR( SIMD_SLLI ( MhOrMv, 3 ), SIMD_SLLI ( DmOrImBit, 7 ) ) )
#define PACK_DPTABLE(MhOrMv, Md, DmOrImBit) SIMD_ADDS_EPU( SIMD_SHUFFLE(Md_Shuffle, Md), SIMD_ADDS_EPU( SIMD_SHUFFLE( MhOrMv_Shuffle, MhOrMv), DmOrImBit ) )



const SIMDtype SIMD_ALL_ZERO = SIMD_SET1(0);
const SIMDtype SIMD_ALL_ONE = SIMD_SET1(1);
const SIMDtype SIMD_ALL_255 = SIMD_SET1(0xFFu);
const SIMDtype BITS_MASK_4_OF_8 = SIMD_SET1((1 << 4) - 1);
const SIMDtype BITS_MASK_3_OF_8 = SIMD_SET1((1 << 3) - 1);
const SIMDtype BITS_MASK_1_OF_8 = SIMD_SET1((1 << 1) - 1);

#define OUT_OF_RANGE_BIT_MASK_REF 4   // Add this to reference character if out of reference's range
#define OUT_OF_RANGE_BIT_MASK_READ 8  // Add this to read character if out of read's range - Reserve only at the moment

// All SIMDtype are 32 x unsigned 8 bit integers
// matchScore = 1, gapExtendScore = -1, 2 * gapOpenScore <= mismatchScore <= gapExtendScore
// minMismatch <= mismatchScore, minSGO <= gapOpenScore < gapExtendScore; where minMismatch and minSGO are constants set in the program
// refSequence and readSequence can only accept 0, 1, 2, 3
int GenerateDPTable(SIMDtype * refSequence, const SIMDtype & refLengthHi, const SIMDtype & refLengthLo, const int maxRefLength,
    SIMDtype * readSequence, const SIMDtype & readLengthHi, const SIMDtype & readLengthLo, const int maxReadLength,
    const int mismatchScore, const int gapOpenScore,
    const SIMDtype & scoreThreshold,
    const int softClipLengthLeft, const SIMDtype & softClipLengthRight,
    SIMDtype & maxScoreHi, SIMDtype & maxScoreLo, SIMDtype & hitPosHi, SIMDtype & hitPosLo, SIMDtype & softClipRight, SIMDtype & maxScoreCount,
    SIMDtype * __restrict DPTable, SIMDtype * __restrict DPTempRow)
{
    //
    // readLength < 255 + gapOpenScore - 1 + scoreThreshold
    // scoreThreshold must <= readLength, and 1 <= scoreThreshold <= 255
    //
    // DPTable Dimension: 1 reference character vs reads as 1 row; 1 read character vs reference as 1 column
    //
    // DP scores are defined as (they are always >= 0):
    //  Im[i, j] = Max(I[i, j] - M[i, j] - SGO + SGE, 0)  =  Max(I[i, j] + SGE, M[i, j] + SGO) - (M[i, j] + SGO)
    //  Dm[i, j] = Max(D[i, j] - M[i, j] - SGO + SGE, 0)  =  Max(D[i, j] + SGE, M[i, j] + SGO) - (M[i, j] + SGO)
    //  Mv[i, j] = M[i, j] - M[i - 1, j] - SGO
    //  Mh[i, j] = M[i, j] - M[i, j - 1] - SGO
    //  Md[i, j] = M[i, j] - M[i - 1, j - 1] - SGO
    //
    // Recurive formula as follows:
    //  D_Delta = Mh[i - 1, j] + Dm[i - 1, j] + 2SGO    -> D_Delta is the score delta from M[i-1, j-1] to [i, j] through "Deletion" path
    //  I_Delta = Mv[i, j - 1] + Im[i, j - 1] + 2SGO    -> I_Delta is the score delta from M[i-1, j-1] to [i, j] through "Insertion" path
    //  M_Delta = MatchScore/MismatchScore              -> M_Delta is the score delta from M[i-1, j-1] to [i, j] through "Match/Mismatch" path
    //  Md[i, j] = Max(I_Delta, D_Delta, S(i,j) ) - SGO
    //  Dm[i, j] = Max(D_Delta - Md[i, j] + SGE, 0)
    //  Im[i, j] = Max(I_Delta - Md[i, j] + SGE, 0)
    //  Mv[i, j] = Md[i, j] - Mh[i - 1, j] - SGO
    //  Mh[i, j] = Md[i, j] - Mv[i, j - 1] - SGO
    //
    // All variables above are always positive
    //
    // During backtracking, Dm == SGE - SGO and Im == SGE - SGO imply that alignment can go back on "Deletion" path and "Insertion" path respectively
    // If backtracking cannot go back on "Match/Mismatch" path, either Dm or Im must be equal to SGE - SGO
    // Also, during backtracking and when a gap has been initiated (in a cell with Dm == SGE - SGO or Im == SGE - SGO), 
    // Dm > 0 and Im > 0 imply that the gap should be extended
    //
    //  M = optimal score
    //  M is stored as (M + readLength - scoreThreshold + 1 - col)
    //    (M + readLength - scoreThreshold + 1 - col) <= 255 (with matchScore == 1 and read length <= 255 + scoreThreshold - 1) 
    //    Any cells on column col can add at most readLength - col to the alignment score in the alignments that follow (with matchScore == 1) 
    //    Therefore, any cells with (M + readLength - scoreThreshold + 1 - col) <= 0 will not be part of any alignment with score >= scoreThreshold
    //    and (M + readLength - scoreThreshold + 1 - col) can be trimmed at zero without affecting alignment result
    //
    //  In general, the score trimmed will be handled similarly as left softClip (i.e. setting score to 0 in local alignment).
    //  Delta scores keep track of score delta after the score trimming.
    //  If the delta scores are used to recover the DP scores, the DP scores after score trimming will be recovered.
    //
    //  Let ST stand for scoreThreshold
    //  Let STM1 as scoreThreshold - 1
    //  Let RL as readLength
    //  M is stored as (M + RL - STM1 - col)
    //
    // Cells in DPTempRow (DPTempRow is used for DP):
    //    0th to 3nd bit = Mh
    //    4rd to 7th bit = Dm
    //
    // DPTempRow must be allocated with SIMDtype * (maxReadLength + 1)  - DPTempRow is used internally as temporary variables only
    //
    // Cells in DPTable (DPTable is used for backtracking):
    //    (Md - SGO + mismatchScore) x 42 + Mh x 3 + (0 -> softclipped, 1 -> Dm == SGE - SGO, 2 -> Dm != SGE - SGO)
    //
    // DPTable must be allocated with SIMDtype * (maxReadLength + 1) * (maxRefLength + 1)
    //
    // True readLength = readLengthHi * 256 + readLengthLo
    // True refLength = refLengthHi * 256 + refLengthLo
    // True maxScore = maxScoreHi * 256 + maxScoreLo
    // True hitPos = hitPosHi * 256 + hitPosLo
    //
    // All SIMDtype variables are 32 x 8 bit unsigned
    // Saturation arithmetics is used extensively in this module and the order of arithmetics operations are important
    //

    const char moduleName[] = "GenerateDPTable_SIMD_256_DPBit8";

    // Fixed scores
    const int gapExtendScore = -1;
    const int matchScore = 1;
    const int minSGO = -6;
    const int minMismatch = -4;

    if (mismatchScore > gapExtendScore || mismatchScore < gapOpenScore * 2 || mismatchScore < minMismatch || gapOpenScore < minSGO || gapOpenScore >= gapExtendScore)
    {
        fprintf(stderr, "%s : mismatchScore > gapExtendScore || mismatchScore < gapOpenScore * 2 || mismatchScore < minMismatch || gapOpenScore < minSGO || gapOpenScore >= gapExtendScore\n", moduleName);
        return EXIT_FAILURE;
    }

    // Derived scores
    const SIMDtype negGapOpenScore = SIMD_SET1(0 - gapOpenScore);
    const SIMDtype negGapExtendScore = SIMD_SET1(0 - gapExtendScore);
    const SIMDtype negGapOpenScorePlus1 = SIMD_ADDS_EPI(negGapOpenScore, SIMD_ALL_ONE);
    const SIMDtype negGapOpenScoreX2 = SIMD_ADDS_EPI(negGapOpenScore, negGapOpenScore);
    const SIMDtype negGapOpenScoreX2Plus1 = SIMD_ADDS_EPI(negGapOpenScoreX2, SIMD_ALL_ONE);
    const SIMDtype negGapInitScore = SIMD_SUB_EPI(negGapOpenScore, negGapExtendScore);
    const SIMDtype negGapOpenScoreAddMismatchScore = SIMD_SUB_EPI(negGapOpenScore, SIMD_SET1(0 - mismatchScore));

    // SGO = GapOpenScore; SGE = GapExtendScore

    SIMDtype matchMismatchScoreMinus2SGO_Shuffle;  // Shuffle mask = ref char XOR read char
    SIMDtype DmBitMaskShuffle;                     // Shuffle mask = Dm; if Dm == NegGapInitScore, set to 1, 2 otherwise
    SIMDtype Md_Shuffle;                           // Shuffle mask = Md_MinusSGO, set to (Md + SGO - mismatchScore) * 42
    SIMDtype MhOrMv_Shuffle;                       // Shuffle mask = MhOrMv, set to MhOrMv * 3
    {
        //for-z-loop is for hanlding sse/avx
        //as avx instructions can only shuffle within 128-bit lanes 
        matchMismatchScoreMinus2SGO_Shuffle = SIMD_SET1(mismatchScore - gapOpenScore * 2);
        unsigned char * matchMismatchScoreMinus2SGO_ShuffleArray = (unsigned char *)(&matchMismatchScoreMinus2SGO_Shuffle);
        for (int z = 0; z<N_DP_PACKED; z += 16) matchMismatchScoreMinus2SGO_ShuffleArray[z] = matchScore - gapOpenScore * 2;

        DmBitMaskShuffle = SIMD_SET1(2);
        unsigned char * DmBitMaskShuffleArray = (unsigned char *)(&DmBitMaskShuffle);
        for (int z = 0; z<N_DP_PACKED; z += 16) DmBitMaskShuffleArray[z + gapExtendScore - gapOpenScore] = 1;

        Md_Shuffle = SIMD_ALL_ZERO;
        unsigned char * Md_ShuffleArray = (unsigned char *)(&Md_Shuffle);
        for (int i = mismatchScore - gapOpenScore * 2; i < 16; ++i)
            for (int z = 0; z<N_DP_PACKED; z += 16)
                Md_ShuffleArray[z + i] = (i + gapOpenScore * 2 - mismatchScore) * 42;

        MhOrMv_Shuffle = SIMD_ALL_ZERO;
        unsigned char * MhOrMv_ShuffleArray = (unsigned char *)(&MhOrMv_Shuffle);
        for (int i = 0; i < 16; ++i)
            for (int z = 0; z<N_DP_PACKED; z += 16)
                MhOrMv_ShuffleArray[z + i] = i * 3;
    }

    // Constant values below are precalculated for the corresponding swap operations using _mm256_shuffle_epi8
    // for calculating horizontal max/min
    const SIMDtype swap8BitMask = SIMD_SET_SHUFFLE_MASK(0x0E0F0C0D, 0x0A0B0809, 0x06070405, 0x02030001);
    const SIMDtype swap16BitMask = SIMD_SET_SHUFFLE_MASK(0x0D0C0F0E, 0x09080B0A, 0x05040706, 0x01000302);
    const SIMDtype swap32BitMask = SIMD_SET_SHUFFLE_MASK(0x0B0A0908, 0x0F0E0D0C, 0x03020100, 0x07060504);
    const SIMDtype swap64BitMask = SIMD_SET_SHUFFLE_MASK(0x07060504, 0x03020100, 0x0F0E0D0C, 0x0B0A0908);


    // Input validation and determine readPosToTakeMax
    int minReadPosToTakeMax;
    int minReadPosAllCanTakeMax;
    {
        // Check refLength against maximum refLength
        SIMDtype refLength16Bit;
        refLength16Bit = SIMD_MAX16_EPU(SIMD_UNPACKLO(refLengthLo, refLengthHi), SIMD_UNPACKHI(refLengthLo, refLengthHi));
        refLength16Bit = SIMD_MAX16_EPU(refLength16Bit, SIMD_SHUFFLE(refLength16Bit, swap16BitMask));
        refLength16Bit = SIMD_MAX16_EPU(refLength16Bit, SIMD_SHUFFLE(refLength16Bit, swap32BitMask));
        refLength16Bit = SIMD_MAX16_EPU(refLength16Bit, SIMD_SHUFFLE(refLength16Bit, swap64BitMask));
        __m128i refLength16Bit128 = SIMD_MAX128(refLength16Bit);
        if ((uint16_t)_mm_extract_epi16(refLength16Bit128,0) > maxRefLength)
        {
            fprintf(stderr, "%s : refLength > maxRefLength\n", moduleName);
            return EXIT_FAILURE;
        }

        // Check readLength against maximum readLength
        SIMDtype readLength16Bit_0 = SIMD_UNPACKLO(readLengthLo, readLengthHi);
        SIMDtype readLength16Bit_1 = SIMD_UNPACKHI(readLengthLo, readLengthHi);
        SIMDtype readLength16Bit;
        readLength16Bit = SIMD_MAX16_EPU(readLength16Bit_0, readLength16Bit_1);
        readLength16Bit = SIMD_MAX16_EPU(readLength16Bit, SIMD_SHUFFLE(readLength16Bit, swap16BitMask));
        readLength16Bit = SIMD_MAX16_EPU(readLength16Bit, SIMD_SHUFFLE(readLength16Bit, swap32BitMask));
        readLength16Bit = SIMD_MAX16_EPU(readLength16Bit, SIMD_SHUFFLE(readLength16Bit, swap64BitMask));
        __m128i readLength16Bit128 = SIMD_MAX128(readLength16Bit);
        int *maxReadLength32Bit = (int*)(&readLength16Bit128);
        if ((uint16_t)_mm_extract_epi16(readLength16Bit128,0) > maxReadLength)
        {
            fprintf(stderr, "%s : readLength > maxReadLength\n", moduleName);
            return EXIT_FAILURE;
        }

        SIMDtype scoreThreshold_0 = SIMD_UNPACKLO(scoreThreshold, SIMD_ALL_ZERO);
        SIMDtype scoreThreshold_1 = SIMD_UNPACKHI(scoreThreshold, SIMD_ALL_ZERO);

        // Check scoreThreshold against readLength
        if (SIMD_TESTZ(SIMD_SUBS16_EPU(scoreThreshold_0, readLength16Bit_0), SIMD_ALL_255) != 1 ||
            SIMD_TESTZ(SIMD_SUBS16_EPU(scoreThreshold_1, readLength16Bit_1), SIMD_ALL_255) != 1)
        {
            fprintf(stderr, "%s : scoreThreshold > readLength\n", moduleName);
            return EXIT_FAILURE;
        }

        // Check that scoreThreshold must > 0
        if (SIMD_TESTZ(SIMD_SUBS_EPU(SIMD_ALL_ONE, scoreThreshold), SIMD_ALL_255) != 1)
        {
            fprintf(stderr, "%s : scoreThreshold == 0\n", moduleName);
            return EXIT_FAILURE;
        }

        // Check that scoreThreshold must < 255 - (gapExtendScore - gapOpenScore) / (gapExtendScore - mismatchScore) + 1;
        int scoreTakingBuffer = (gapExtendScore - gapOpenScore) / (gapExtendScore - mismatchScore) + 1;
        if (SIMD_TESTZ(SIMD_SUBS_EPU(scoreThreshold, SIMD_SET1(255 - scoreTakingBuffer)), SIMD_ALL_255) != 1)
        {
            fprintf(stderr, "%s : scoreThreshold > 255 - (gapExtendScore - gapOpenScore) / (gapExtendScore - mismatchScore) + 1\n", moduleName);
            return EXIT_FAILURE;
        }

        // Check maximum supported read length against readLength
        SIMDtype maxSupportedReadLength_0 = SIMD_ADDS16_EPU(SIMD_UNPACKLO(SIMD_SET1(255 + gapOpenScore - 2), SIMD_ALL_ZERO), scoreThreshold_0);
        SIMDtype maxSupportedReadLength_1 = SIMD_ADDS16_EPU(SIMD_UNPACKLO(SIMD_SET1(255 + gapOpenScore - 2), SIMD_ALL_ZERO), scoreThreshold_1);
        if (SIMD_TESTZ(SIMD_SUBS16_EPU(readLength16Bit_0, maxSupportedReadLength_0), SIMD_ALL_255) != 1 ||
            SIMD_TESTZ(SIMD_SUBS16_EPU(readLength16Bit_1, maxSupportedReadLength_1), SIMD_ALL_255) != 1)
        {
            fprintf(stderr, "%s : readLength >= 255 + gapOpenScore - 1 + scoreThreshold\n", moduleName);
            return EXIT_FAILURE;
        }

        // Determine minReadPosToTakeMax and minReadPosAllCanTakeMax
        SIMDtype readLengthMinusSCRight_0 = SIMD_SUBS16_EPU(readLength16Bit_0, SIMD_UNPACKLO(softClipLengthRight, SIMD_ALL_ZERO));
        SIMDtype readLengthMinusSCRight_1 = SIMD_SUBS16_EPU(readLength16Bit_1, SIMD_UNPACKHI(softClipLengthRight, SIMD_ALL_ZERO));
        SIMDtype readPosToTakeMax_0 = SIMD_MAX16_EPU(readLengthMinusSCRight_0, scoreThreshold_0);
        SIMDtype readPosToTakeMax_1 = SIMD_MAX16_EPU(readLengthMinusSCRight_1, scoreThreshold_1);

        SIMDtype minReadPosToTakeMax16Bit;
        minReadPosToTakeMax16Bit = SIMD_MIN16_EPU(readPosToTakeMax_0, readPosToTakeMax_1);
        minReadPosToTakeMax16Bit = SIMD_MIN16_EPU(minReadPosToTakeMax16Bit, SIMD_SHUFFLE(minReadPosToTakeMax16Bit, swap16BitMask));
        minReadPosToTakeMax16Bit = SIMD_MIN16_EPU(minReadPosToTakeMax16Bit, SIMD_SHUFFLE(minReadPosToTakeMax16Bit, swap32BitMask));
        minReadPosToTakeMax16Bit = SIMD_MIN16_EPU(minReadPosToTakeMax16Bit, SIMD_SHUFFLE(minReadPosToTakeMax16Bit, swap64BitMask));
        __m128i minReadPosToTakeMax16Bit128 = SIMD_MIN128(minReadPosToTakeMax16Bit);
        minReadPosToTakeMax = (uint16_t)_mm_extract_epi16(minReadPosToTakeMax16Bit128,0);

        SIMDtype maxReadPosToTakeMax16Bit;
        maxReadPosToTakeMax16Bit = SIMD_MAX16_EPU(readPosToTakeMax_0, readPosToTakeMax_1);
        maxReadPosToTakeMax16Bit = SIMD_MAX16_EPU(maxReadPosToTakeMax16Bit, SIMD_SHUFFLE(maxReadPosToTakeMax16Bit, swap16BitMask));
        maxReadPosToTakeMax16Bit = SIMD_MAX16_EPU(maxReadPosToTakeMax16Bit, SIMD_SHUFFLE(maxReadPosToTakeMax16Bit, swap32BitMask));
        maxReadPosToTakeMax16Bit = SIMD_MAX16_EPU(maxReadPosToTakeMax16Bit, SIMD_SHUFFLE(maxReadPosToTakeMax16Bit, swap64BitMask));
        __m128i maxReadPosToTakeMax16Bit128 = SIMD_MAX128(maxReadPosToTakeMax16Bit);
        // (gapExtendScore - gapOpenScore) / (gapExtendScore - mismatchScore) + 1 is added to make sure no score at out-of-ref-range positions will exceed or having the same score as in-ref-range positions
        minReadPosAllCanTakeMax = (uint16_t)_mm_extract_epi16(maxReadPosToTakeMax16Bit128,0) + (gapExtendScore - gapOpenScore) / (gapExtendScore - mismatchScore) + 1;
    }

    // Pre-calculate frequently used numbers

    SIMDtype RL_MinusSTM1;
    {
        SIMDtype tempWithHiBit = SIMD_NAND(SIMD_CMPEQ(readLengthHi, SIMD_ALL_ZERO),
            SIMD_ADDS_EPU(SIMD_SET1(2),
                SIMD_ADDS_EPU(readLengthLo, SIMD_SUBS_EPU(SIMD_ALL_255, scoreThreshold))));
        SIMDtype tempNoHiBit = SIMD_AND(SIMD_CMPEQ(readLengthHi, SIMD_ALL_ZERO),
            SIMD_ADDS_EPU(SIMD_ALL_ONE, SIMD_SUBS_EPU(readLengthLo, scoreThreshold)));
        RL_MinusSTM1 = SIMD_OR(tempWithHiBit, tempNoHiBit);
    }

    const SIMDtype STM1 = SIMD_SUBS_EPU(scoreThreshold, SIMD_ALL_ONE);
    const SIMDtype _255MinusSTM1 = SIMD_SUBS_EPU(SIMD_ALL_255, STM1);

    // Initially, softClipRightPosMinusSTM1 is set as below; when DP goes out of ref range, it will be set to 0xFF to stop taking maxScore
    SIMDtype softClipRightPosMinusSTM1 = SIMD_SUBS_EPU(RL_MinusSTM1, softClipLengthRight);

    // Initialize maxScore and read position
    SIMDtype maxScore_MinusSTM1 = SIMD_ALL_ZERO;           // Plus 1 so 0 implies no max score has been taken
    SIMDtype maxScore_MinusSTM1_LastRow = SIMD_ALL_ZERO;
    SIMDtype maxScorePos_MinusSTM1 = SIMD_ALL_ZERO;     // Read position giving maxScore
                                                          
    SIMDtype * DPTableCurr = DPTable;    // Set current pointer to DPTable

    // Filling the first row

    // Filling the first cell of the first row

    {
        // DPTable is undefined at the 0th column of the 0th row and should never be used
        SIMD_STORE(DPTableCurr++, SIMD_ALL_ZERO);
    }

    // Set current pointer to DPTemp; DPTemp values for ith column is written to the (i-1)th cell in DPTemp
    SIMDtype * DPTempRowCurr = DPTempRow;

    // Filling other cells in the first row, i.e. i = 0

    {
        SIMDtype M_AddRL_MinusSTM1_MinusCol = RL_MinusSTM1;
        SIMDtype M_AddRL_MinusSTM1_MinusCol_NegDelta;
        SIMDtype M_Increment;

        for (int j = 1; j <= maxReadLength; ++j)
        {
            // Need to keep track of M_AddRL_MinusSTM1_MinusCol to track for score trimming; trimmed score will be used to adjust delta scores
            if (j <= softClipLengthLeft)
            {
                M_AddRL_MinusSTM1_MinusCol_NegDelta = SIMD_ALL_ONE;
                M_Increment = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol_NegDelta, M_AddRL_MinusSTM1_MinusCol);
                M_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol_NegDelta);
                // DPTemp: Mh = negGapOpenScore; Dm = 0
                SIMD_STORE(DPTempRowCurr, PACK_DPTEMP(negGapOpenScore, SIMD_ALL_ZERO));
                // DPTable: Mh = negGapOpenScore; Md is undefined; Dm != SGO - SGE; softclipped
                SIMD_STORE(DPTableCurr, PACK_DPTABLE(negGapOpenScore, SIMD_ALL_ZERO, SIMD_ALL_ZERO));
            }
            else
            {
                if (j == softClipLengthLeft + 1)
                {
                    M_AddRL_MinusSTM1_MinusCol_NegDelta = negGapOpenScorePlus1;
                    M_Increment = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol_NegDelta, M_AddRL_MinusSTM1_MinusCol);
                    M_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol_NegDelta);
                    // DPTemp: Mh = 0; Dm = 0
                    SIMD_STORE(DPTempRowCurr, PACK_DPTEMP(SIMD_ALL_ZERO, SIMD_ALL_ZERO));
                    // DPTable: Mh = 0; Md is undefined; Dm != SGO - SGE; not softclipped
                    SIMD_STORE(DPTableCurr, PACK_DPTABLE(SIMD_ALL_ZERO, SIMD_ALL_ZERO, SIMD_ALL_ONE));
                }
                else
                {
                    M_AddRL_MinusSTM1_MinusCol_NegDelta = SIMD_ADDS_EPU(SIMD_ALL_ONE, SIMD_ALL_ONE);
                    M_Increment = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol_NegDelta, M_AddRL_MinusSTM1_MinusCol);
                    M_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol_NegDelta);
                    // DPTemp: Mh = negGapInitScore; Dm = 0
                    SIMD_STORE(DPTempRowCurr, PACK_DPTEMP(negGapInitScore, SIMD_ALL_ZERO));
                    // DPTable: Mh = negGapInitScore; Md is undefined; Dm != SGO - SGE; not softclipped
                    SIMD_STORE(DPTableCurr, PACK_DPTABLE(negGapInitScore, SIMD_ALL_ZERO, SIMD_ALL_ONE));
                }
            }
            DPTableCurr++;
            DPTempRowCurr++;
        }
    }

    // Filling other rows

    for (int i = 1; i <= maxRefLength; ++i)
    {
        SIMDtype refChar = SIMD_LOAD(refSequence + i);          // refSequence is 1-based
        {
            // Mask reference character if out of reference's range so the character will always mismatch with read characters
            SIMDtype vectorI_Hi = SIMD_SET1(i / 256);
            SIMDtype vectorI_Lo = SIMD_SET1(i % 256);
            SIMDtype temp_HiLessEqual = SIMD_CMPEQ(SIMD_MIN_EPU(refLengthHi, vectorI_Hi), vectorI_Hi);
            SIMDtype temp_HiNotEqual = SIMD_XOR(SIMD_ALL_255, SIMD_CMPEQ(refLengthHi, vectorI_Hi));
            SIMDtype temp_LoLessEqual = SIMD_CMPEQ(SIMD_MIN_EPU(refLengthLo, vectorI_Lo), vectorI_Lo);
            SIMDtype notInRefRange = SIMD_XOR(SIMD_ALL_255, SIMD_AND(temp_HiLessEqual, SIMD_OR(temp_HiNotEqual, temp_LoLessEqual)));
            // Set softClipRightPosMinusMinSTM1 when out of ref range to stop maxScore taking
            softClipRightPosMinusSTM1 = SIMD_OR(notInRefRange, softClipRightPosMinusSTM1);
            // Mask reference character if out of reference's range so the character will always mismatch with read characters
            // This is necessary as softClipRightPosMinusSTM1 is not checked after read position reaches minReadPosAllCanTakeMax
            refChar = SIMD_OR(refChar, SIMD_AND(notInRefRange, SIMD_SET1(OUT_OF_RANGE_BIT_MASK_REF)));
        }

        // Semi-global alignment
        SIMDtype Im = SIMD_ALL_ZERO;
        SIMDtype Mv = negGapOpenScore;
        SIMDtype zero_AddRL_MinusSTM1_MinusCol = RL_MinusSTM1;
        SIMDtype M_PrevRow_AddRL_MinusSTM1_MinusCol = RL_MinusSTM1;

        // Filling the first cell of the first row
        {
            // DPTable: Mh is undefined; Md is undefined;  Dm != SGO - SGE; softclipped
            SIMD_STORE(DPTableCurr++, SIMD_ALL_ZERO);
        }

        DPTempRowCurr = DPTempRow;     // Go back to the first cells in DPTempRow
        SIMDtype TempDPCell = SIMD_LOAD(DPTempRowCurr);

        for (int j = 1; j <= maxReadLength; j++)
        {
            SIMDtype readChar = SIMD_LOAD(readSequence + j);   // readSequence is 1-based
            SIMDtype refCharXORReadChar = SIMD_XOR(refChar, readChar);
            SIMDtype M_Delta_Minus2SGO = SIMD_SHUFFLE(matchMismatchScoreMinus2SGO_Shuffle, refCharXORReadChar);

            SIMDtype M_PrevRowPrevCol_MinusSTM1_MinusCol = M_PrevRow_AddRL_MinusSTM1_MinusCol;

            zero_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(zero_AddRL_MinusSTM1_MinusCol, SIMD_ALL_ONE);

            SIMDtype Mh = UNPACK_DPTEMP_MhOrMv(TempDPCell);
            SIMDtype Dm = UNPACK_DPTEMP_DmOrIm(TempDPCell);

            // Read length < 255 + gapOpenScore - 1 + ST so the following statement will not saturate above
            // This statement will not saturate below as Mh stored in DPTemp is the delta between saturated scores
            M_PrevRow_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(SIMD_ADDS_EPU(M_PrevRow_AddRL_MinusSTM1_MinusCol, Mh), negGapOpenScorePlus1);

            SIMDtype D_Delta_Minus2SGO, I_Delta_Minus2SGO, ID_Delta_Minus2SGO;
            D_Delta_Minus2SGO = SIMD_ADDS_EPU(Mh, Dm);
            I_Delta_Minus2SGO = SIMD_ADDS_EPU(Mv, Im);
            ID_Delta_Minus2SGO = SIMD_MAX_EPU(D_Delta_Minus2SGO, I_Delta_Minus2SGO);
            SIMDtype Md_MinusSGO = SIMD_MAX_EPU(ID_Delta_Minus2SGO, M_Delta_Minus2SGO);

            SIMDtype M_AddRL_MinusSTM1_MinusCol;
            SIMDtype M_Increment;
            {
                // Calculating M_AddRL_MinusSTM1_MinusCo = M_PrevRowPrevCol_MinusSTM1_MinusCol + Md_MinusSGO + 2SGO - 1 (note that Md = M[i-1] - M[i-1,j-1] - SGO)

                // Need to keep the score being trimmed in M_Increment to make sure delta scores are deltas between saturated scores
                // Note that Md_MinusSGO + 2SGO - 1 <= 0 so negGapOpenScoreX2Plus1 - Md_MinusSGO >= 0
                SIMDtype M_AddRL_MinusSTM1_MinusCol_NegDelta = SIMD_SUBS_EPU(negGapOpenScoreX2Plus1, Md_MinusSGO);
                M_Increment = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol_NegDelta, M_PrevRowPrevCol_MinusSTM1_MinusCol);
                M_AddRL_MinusSTM1_MinusCol = SIMD_SUBS_EPU(M_PrevRowPrevCol_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol_NegDelta);
            }

            if (j <= softClipLengthLeft)
            {
                M_Increment = SIMD_ADDS_EPU(M_Increment, SIMD_SUBS_EPU(zero_AddRL_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol));
                // There is no need to add SIMD_SUBS_EPU(zero_AddRL_MinusSTM1_MinusCol, M_AddRL_MinusSTM1_MinusCol) to M_AddRL_MinusSTM1_MinusCol
                // as M_AddRL_MinusSTM1_MinusCol will only be used to keep track of maxScore and if M_Increment > 0, M_AddRL_MinusSTM1_MinusCol + M_Increment cannot be a maxscore 
            }

            Md_MinusSGO = SIMD_ADDS_EPU(Md_MinusSGO, M_Increment);
            {
                SIMDtype Md_MinusSGE = SIMD_SUBS_EPU(Md_MinusSGO, negGapInitScore);
                Im = SIMD_SUBS_EPU(I_Delta_Minus2SGO, Md_MinusSGE);
                Dm = SIMD_SUBS_EPU(D_Delta_Minus2SGO, Md_MinusSGE);
            }
            {
                SIMDtype temp_Mv = Mv;
                Mv = SIMD_SUBS_EPU(Md_MinusSGO, Mh);
                Mh = SIMD_SUBS_EPU(Md_MinusSGO, temp_Mv);
            }

            SIMD_STORE(DPTempRowCurr, PACK_DPTEMP(Mh, Dm));
            DPTempRowCurr++;
            TempDPCell = SIMD_LOAD(DPTempRowCurr);

            {
                SIMDtype DmBit = SIMD_SHUFFLE(DmBitMaskShuffle, Dm);
                // softClipDmMask also considered score trimming but it's ok as trimmed score will never be used
                SIMDtype softClipDmMask = SIMD_CMPEQ(M_Increment, SIMD_ALL_ZERO);
                DmBit = SIMD_AND(softClipDmMask, DmBit);  // 0 -> softclipped, 1 -> Dm == SGE - SGO, 2 -> Dm != SGE - SGO
                SIMD_STORE(DPTableCurr, PACK_DPTABLE(Mh, Md_MinusSGO, DmBit));
            }

            DPTableCurr++;

            if (j >= minReadPosToTakeMax)
            {
                SIMDtype vectorJ_MinusSTM1;
                
                // Recover read position
                if (j <= 255)
                {
                    vectorJ_MinusSTM1 = SIMD_SUBS_EPU(SIMD_SET1(j), STM1);
                }
                else
                {
                    vectorJ_MinusSTM1 = SIMD_ADDS_EPU(SIMD_SET1(j - 255), _255MinusSTM1);
                    // j - STM1 >= 255 will make canTakeMax == true but it is ok as readLength < 255 + gapOpenScore - 1 + scoreThreshold
                    // so, j >= readLength + 1 - gapOpenScore and read positions after readLength always generate mismatches that score will be less than maxScore
                }

                // Recover score
                SIMDtype M_MinusSTM1 = SIMD_SUBS_EPU(M_AddRL_MinusSTM1_MinusCol, SIMD_SUBS_EPU(RL_MinusSTM1, vectorJ_MinusSTM1));

                if (j < minReadPosAllCanTakeMax)
                {
                    // As out-of-range readSequence is loaded with unmatchable characters, there is no need to check read length
                    SIMDtype softClipRightPosMinusVectorJ = SIMD_SUBS_EPU(softClipRightPosMinusSTM1, vectorJ_MinusSTM1);
                    // softClipRightPosMinusMinSTM1 is updated to 0xFF when out of ref range so there is no need to check ref length

                    SIMDtype canTakeMax = SIMD_CMPEQ(softClipRightPosMinusVectorJ, SIMD_ALL_ZERO);

                    // Set score to zero if cannot take max
                    M_MinusSTM1 = SIMD_AND(canTakeMax, M_MinusSTM1);
                }

                // Update maxScore; maxIncrement is used to check if maxScore is exceeded
                SIMDtype maxIncrement = SIMD_SUBS_EPU(M_MinusSTM1, maxScore_MinusSTM1);
                maxScore_MinusSTM1 = SIMD_ADDS_EPU(maxScore_MinusSTM1, maxIncrement);

                // Update maxScorePos to current position and maxScoreCount to zero when score > maxScore
                SIMDtype M_LessThanOrEqualMaxScore = SIMD_CMPEQ(maxIncrement, SIMD_ALL_ZERO);
                maxScorePos_MinusSTM1 = SIMD_BLENDV(vectorJ_MinusSTM1, maxScorePos_MinusSTM1, M_LessThanOrEqualMaxScore);
                maxScoreCount = SIMD_AND(M_LessThanOrEqualMaxScore, maxScoreCount);

                // Add 1 to maxScoreCount when score >= maxScore
                SIMDtype M_GreaterThanOrEqualMaxScore = SIMD_CMPEQ(M_MinusSTM1, maxScore_MinusSTM1);
                maxScoreCount = SIMD_ADDS_EPU(maxScoreCount, SIMD_AND(M_GreaterThanOrEqualMaxScore, SIMD_ALL_ONE));
            }
        }
        {
            // Update hitPos if maxScore was updated
            SIMDtype maxScoreNotUpdated = SIMD_CMPEQ(maxScore_MinusSTM1, maxScore_MinusSTM1_LastRow);
            SIMDtype vectorI_Hi = SIMD_SET1(i / 256);
            SIMDtype vectorI_Lo = SIMD_SET1(i % 256);
            hitPosHi = SIMD_BLENDV(vectorI_Hi, hitPosHi, maxScoreNotUpdated);  // if maxScoreNotUpdated is true, get hitPostHi; 2 cycles
            hitPosLo = SIMD_BLENDV(vectorI_Lo, hitPosLo, maxScoreNotUpdated);  // if maxScoreNotUpdated is true, get hitPostLo; 2 cycles
            maxScore_MinusSTM1_LastRow = maxScore_MinusSTM1;
        }
    }
    {
        // Format output parameters
        softClipRight = SIMD_SUBS_EPU(RL_MinusSTM1, maxScorePos_MinusSTM1);

        SIMDtype maxScore_MinusST = SIMD_SUBS_EPU(maxScore_MinusSTM1, SIMD_ALL_ONE);
        SIMDtype STM1 = SIMD_SUBS_EPU(scoreThreshold, SIMD_ALL_ONE);
        SIMDtype maxScoreGT255 = SIMD_CMPEQ(SIMD_MAX_EPU(maxScore_MinusST, SIMD_SUBS_EPU(SIMD_ALL_255, STM1)), maxScore_MinusST);

        maxScoreHi = SIMD_AND(SIMD_ALL_ONE, maxScoreGT255);

        // maxScoreLoIfMaxScoreGT255 = (maxScore - ST) - 255 + ST - 1 = maxScore - 256
        SIMDtype maxScoreLoIfMaxScoreGT255 = SIMD_SUBS_EPU(maxScore_MinusST, SIMD_SUBS_EPU(SIMD_ALL_255, STM1));
        SIMDtype maxScoreLoIfMaxScoreLE255 = SIMD_ADDS_EPU(maxScore_MinusST, scoreThreshold);

        maxScoreLo = SIMD_BLENDV(maxScoreLoIfMaxScoreLE255, maxScoreLoIfMaxScoreGT255, maxScoreGT255);

        // Set output to zero if maxScore_MinusSTM1 == 0, i.e. no maxScore was set
        SIMDtype maxScoreNotTaken = SIMD_CMPEQ(maxScore_MinusSTM1, SIMD_ALL_ZERO);
        maxScoreHi = SIMD_NAND(maxScoreNotTaken, maxScoreHi);
        maxScoreLo = SIMD_NAND(maxScoreNotTaken, maxScoreLo);
        maxScoreCount = SIMD_NAND(maxScoreNotTaken, maxScoreCount);
        hitPosHi = SIMD_NAND(maxScoreNotTaken, hitPosHi);
        hitPosLo = SIMD_NAND(maxScoreNotTaken, hitPosLo);
        softClipRight = SIMD_NAND(maxScoreNotTaken, softClipRight);
    }
    return EXIT_SUCCESS;
}


void GPUBacktrack ( SIMDtype * DNASequence, short DNALength, uint maxDNALength,
                    SIMDtype * readSequence, short readLength, uint maxReadLength,
                    short & hitPos,
                    short currScore, // this is useless, for debug only
                    uint8_t clipLtCheckLoc, uint8_t clipRtLength,
                    int8_t MatchScore, int8_t MismatchScore,
                    int8_t GapOpenScore, int8_t GapExtendScore,
                    void * DPTable, uchar * pattern, const short iterIdx )
{
    SIMDtype * dptable = ( SIMDtype * ) ( DPTable );
    uchar * curPattern = pattern;
    uint pIndex = 0;
#define MC_PatternAppend(x) { curPattern[pIndex] = (x); ++pIndex; }

/*
//skip backtracking 
MC_PatternAppend ('M');
MC_PatternAppend ('V');
MC_PatternAppend (readLength);
MC_PatternAppend (0);
return;
*/

    if ( clipRtLength > 0 )
    {
        MC_PatternAppend ( 'S' );
        MC_PatternAppend ( 'V' );
        MC_PatternAppend ( clipRtLength );
    }
    short i = readLength - clipRtLength;
    short j = hitPos;
    // start backtracking at (readPos = i, refPos = j)
#define GET_CELL(j,i) extract_epi8(*(dptable+(j)*(maxReadLength+1)+(i)), iterIdx)

    uint8_t readChar = extract_epi8 ( readSequence[i], iterIdx );
    uint8_t refChar = extract_epi8 ( DNASequence[j], iterIdx );
#define MC_NextRefChar() { --j; refChar = extract_epi8 ( DNASequence[j], iterIdx );}
#define MC_NextReadChar() { --i; readChar = extract_epi8 ( readSequence[i], iterIdx );; }


    enum DP_BacktrackState { DP_BT_NORMAL, DP_BT_I_EXT, DP_BT_D_EXT, DP_BT_SM_EXIT, DP_BT_SI_EXIT, DP_BT_SD_EXIT };
    //normal, insertion extension, deletion extension, softclip match/mismatch exit, softclip insertion exit
    DP_BacktrackState state = DP_BT_NORMAL;

    int8_t accum = 0;
    while ( i > 0 && j > 0 )
    {
// solve equation: 42x+3y+z = RHS, (x<6, y<14, z<3, x,y,z nonneg) 
// x = RHS/42
// y = RHS/3%14 
// z = RHS%3

        uint8_t this_score = GET_CELL(j,i); // this is unsigned!
        int8_t flag = this_score % 3;
        int8_t upper_left_flag = GET_CELL(j-1,i-1) % 3;
#define FLAG_SOFTCLIP 0
#define FLAG_DELETION 1
        int8_t hori_diff = GapOpenScore + (int8_t)(this_score / 3 % 14);
        int8_t diag_diff = MismatchScore + (int8_t)(this_score / 42);
        int8_t vert_diff = diag_diff - (int8_t)(GET_CELL(j-1,i) / 3 % 14) - GapOpenScore;
        int8_t match_score = refChar==readChar?MatchScore:MismatchScore;
        if (state == DP_BT_NORMAL) {
            if (upper_left_flag == FLAG_SOFTCLIP && diag_diff == match_score && i!=1) {
                state = DP_BT_SM_EXIT;
                break;
            } else if (diag_diff == match_score) {
                currScore -= diag_diff;
                MC_PatternAppend ( 'm' - ( ( refChar == readChar ) << 5 ) );
                MC_NextRefChar( );
                MC_NextReadChar( );
            } else if (flag == FLAG_DELETION) {
                currScore -= vert_diff;
                if (vert_diff == GapOpenScore ) {
                    MC_PatternAppend ( 'D' );
                    MC_NextRefChar( );
                } else {
                    MC_PatternAppend ( 'D' );
                    MC_NextRefChar ( );
                    accum = vert_diff - GapExtendScore; // accum = -1 if vert_diff == -2
                    state = DP_BT_D_EXT;
                }
            } else {
                currScore -= hori_diff;
                if (hori_diff == GapOpenScore) {
                    MC_PatternAppend ( 'I' );
                    MC_NextReadChar ( );
                } else {
                    MC_PatternAppend ( 'I' );
                    MC_NextReadChar ( );
                    accum = hori_diff - GapExtendScore;
                    state = DP_BT_I_EXT;
                }
            }
        } else if (state == DP_BT_D_EXT) {
            if (GET_CELL(j-1,i)%3==FLAG_SOFTCLIP && vert_diff+accum == GapOpenScore) {
                state = DP_BT_SD_EXIT; 
                break;
            }
            MC_PatternAppend ( 'D' );
            MC_NextRefChar ( );
            currScore -= GapExtendScore;
            if (vert_diff+accum == GapOpenScore ) {
                state = DP_BT_NORMAL;
                currScore -= (GapOpenScore - GapExtendScore);
            } else accum += vert_diff - GapExtendScore; 
        } else {
            if (GET_CELL(j,i-1)%3==FLAG_SOFTCLIP && hori_diff+accum == GapOpenScore) {
                state = DP_BT_SI_EXIT;
                break;
            }
            MC_PatternAppend ( 'I' );
            MC_NextReadChar ( );
            currScore -= GapExtendScore;
            if (hori_diff+accum == GapOpenScore) {
                state = DP_BT_NORMAL;
                currScore -= GapOpenScore - GapExtendScore;
            } else accum += hori_diff - GapExtendScore;
        }
    }

    //last proc
    if ( j == 0 )
    {
        uint scNum = ( clipLtCheckLoc > i ? i : clipLtCheckLoc );

        if ( scNum < i )
        {
            MC_PatternAppend ( 'I' );
            MC_PatternAppend ( 'V' );
            MC_PatternAppend ( i - scNum );
        }

        MC_PatternAppend ( 'S' );
        MC_PatternAppend ( 'V' );
        MC_PatternAppend ( scNum );
    }
    else if ( state == DP_BT_SI_EXIT )
    {
        MC_PatternAppend ( 'I' );
        MC_PatternAppend ( 'S' );
        MC_PatternAppend ( 'V' );
        MC_PatternAppend ( i - 1 );
    }
    else if ( state == DP_BT_SD_EXIT) {
        MC_PatternAppend ( 'D' );
        MC_PatternAppend ( 'S' );
        MC_PatternAppend ( 'V' );
        MC_PatternAppend ( i - 1 );
        MC_PatternAppend ( 0 );
    hitPos = INVALID_CIGAR_STRING_MARK;
    return;
    // ignore this alignment, will set dpScore <- 0
    }
    else if ( state == DP_BT_SM_EXIT )
    {
        MC_PatternAppend ( 'm' - ( ( refChar == readChar ) << 5 ) );
        MC_PatternAppend ( 'S' );
        MC_PatternAppend ( 'V' );
        MC_PatternAppend ( i - 1 );
        j -= 1;
    }

    MC_PatternAppend ( 0 );
    hitPos = j;
}

void SemiGlobalAlignment ( SIMDtype * DNASequence, SIMDtype & DNALengthsHi, SIMDtype & DNALengthsLo, uint maxDNALength,
                           SIMDtype * readSequence, SIMDtype & readLengthsHi, SIMDtype & readLengthsLo, uint maxReadLength,
                           SIMDtype & maxScoresHi, SIMDtype & maxScoresLo, SIMDtype & hitLocsHi, SIMDtype & hitLocsLo,
                           int clipLtSizes, int clipRtSize,
                           int8_t MismatchScore,
                           int8_t GapOpenScore,
                           void * DPTable, SIMDtype* DPTempRow, SIMDtype & maxScoreCounts,
                           SIMDtype & cutOffThreshold, uchar * pattern, uchar numDPInstances)
{
    hitLocsLo = SIMD_ALL_ZERO;
    hitLocsHi = SIMD_ALL_ZERO;
    maxScoreCounts = SIMD_ALL_ZERO;
    SIMDtype clipRtSizes = SIMD_SET1(clipRtSize);
    SIMDtype scRight = SIMD_ALL_ZERO;

    GenerateDPTable ( DNASequence, DNALengthsHi, DNALengthsLo, maxDNALength,
                        readSequence, readLengthsHi, readLengthsLo, maxReadLength,
                        MismatchScore,GapOpenScore,
                        cutOffThreshold, 
                        clipLtSizes, clipRtSizes,
                        maxScoresHi, maxScoresLo, hitLocsHi, hitLocsLo, scRight, maxScoreCounts,
                        (SIMDtype*) DPTable, DPTempRow );

    short hitLocations[N_DP_PACKED]={0};

    int8_t matchScore = 1;
    int8_t mismatchScore = MismatchScore;
    int8_t gapOpenScore = GapOpenScore;
    int8_t gapExtendScore = -1;


    for ( int j=0;j<numDPInstances;++j )
    {
        short maxScore = extract_epi8 ( maxScoresHi, j ) * 256 + extract_epi8(maxScoresLo,j);
        uint8_t cutOffScore = extract_epi8 ( cutOffThreshold, j );
        if ( maxScore >= cutOffScore )
        {
            short DNALength = extract_epi8 ( DNALengthsHi, j )*256 + extract_epi8(DNALengthsLo,j);
            short readLength = extract_epi8 ( readLengthsHi, j ) * 256 + extract_epi8(readLengthsLo,j);
            hitLocations[j] = extract_epi8 ( hitLocsHi, j )* 256 + extract_epi8(hitLocsLo,j);
            uint8_t rightClip = extract_epi8 ( scRight, j );
            
            if (readLength < cutOffScore)
                hitLocations[j] = INVALID_CIGAR_STRING_MARK;
            else 
                GPUBacktrack ( DNASequence, DNALength, maxDNALength,
                                readSequence, readLength, maxReadLength,
                                hitLocations[j], // as a return value
                                maxScore,
                                clipLtSizes, rightClip,
                                matchScore, mismatchScore,
                                gapOpenScore, gapExtendScore,
                                DPTable, pattern + j * ( maxReadLength + maxDNALength ), j );

           if (hitLocations[j] == INVALID_CIGAR_STRING_MARK) {
/*
3 cases:
1. DP start with Deletion
2. DP start with Softclip and then Deletion
3. readlen < cutOffScore

Case 1.2. should disappear if No Anchor is set
*/
//puts("Ignore 1 Invalid Alignment");
                uint8_t* ptr = (uint8_t*)&maxScoresLo;
                ptr[N_DP_PACKED-1-j] = 0;
                ptr = (uint8_t*)&maxScoresHi;
                ptr[N_DP_PACKED-1-j] = 0;
                hitLocations[j] = 0;
            }
        }
    }

    uint8_t bufhi[N_DP_PACKED];
    uint8_t buflo[N_DP_PACKED];

    for (int i=0;i<N_DP_PACKED;i++) {
        bufhi[i] = hitLocations[i] / 256;
        buflo[i] = hitLocations[i] % 256;
    }
    hitLocsHi = MC_SetEpi8 ( bufhi , 0 );
    hitLocsLo = MC_SetEpi8 ( buflo , 0 );
    return;
}

//this is very dummy!
void dummy_pack(uint8_t* bufhi, uint8_t* buflo, uint* source, int offset) {
    for (int i=0;i<N_DP_PACKED;i++) {
        bufhi[i] = source[offset+i] / 256;
        buflo[i] = source[offset+i] % 256;
    }
}

void callDP ( uint * DNASequences, uint * DNALengths, uint maxDNALength,
              uint * readSequences, uint * readLengths, uint maxReadLength,
              int clipLtSizes, int clipRtSizes,
              uint * anchorLeftLocs, uint * anchorRightLocs,
              uint MatchScore, uint MismatchScore,
              uint GapOpenScore, uint GapExtendScore,
              int * cutOffThresholds,
              void * DPTable, SIMDtype* DPTempRow,  uint numDPInstances,
              // output
              int * scores, uint * hitLocs,
              uint * maxScoreCounts, uchar * pattern )
{
//puts(N_DP_PACKED==16?"This is SSE!!!":"This is AVX!!!");
//printf("I have %d dp tables!\n", numDPInstances);
    for (int i=0;i<numDPInstances;i++) {
        if (cutOffThresholds[i] != (int)std::max(readLengths[i] * DP_SCORE_THRESHOLD_RATIO,DP_SCORE_THRESHOLD_LOWER_BOUND)) {fprintf(stderr, "%d %d\n", readLengths[i], cutOffThresholds[i]);exit(1) ;}
    }
    SIMDtype _DNASequence[maxDNALength+1];
    SIMDtype _readSequence[maxReadLength+1];
    memset(_DNASequence, 0, sizeof(SIMDtype) * (maxDNALength+1));
    memset(_readSequence, 0, sizeof(SIMDtype) * (maxReadLength+1));
    
    uint8_t bufferHi[N_DP_PACKED];
    uint8_t bufferLo[N_DP_PACKED];

    SIMDtype _readLengthsHi;
    SIMDtype _DNALengthsHi;
    SIMDtype _cutOffThresholds;
    SIMDtype _DNALengthsLo;
    SIMDtype _readLengthsLo;
    
    SIMDtype _hitLocsHi;
    SIMDtype _hitLocsLo;
    SIMDtype _maxScoresHi;
    SIMDtype _maxScoresLo;
    SIMDtype _maxScoreCounts;

    uint round = ( numDPInstances + N_DP_PACKED -1 ) / N_DP_PACKED;
//!! sizeof DNALengths, readLengths, cutOffThresholds  are 32768
    for (int i=round*N_DP_PACKED-1;i>=numDPInstances;--i) {
        DNALengths[i] = DNALengths[0];
        readLengths[i] = readLengths[0];
        cutOffThresholds[i] = cutOffThresholds[0];
    }

    for ( int i=0; i < round; ++i )
    {
    //printf("round %d\n", i);
        int caseBaseNo = N_DP_PACKED * i;
        int numCase = ( numDPInstances - caseBaseNo > N_DP_PACKED ? N_DP_PACKED : ( numDPInstances - caseBaseNo ) );

        for ( int idx = 0; idx < numCase; ++idx )
        {
            uint8_t* dna_ptr = (uint8_t*) (_DNASequence+1) + N_DP_PACKED -1 - idx;
            uint8_t* read_ptr = (uint8_t*) (_readSequence+1) + N_DP_PACKED -1  - idx;
            int caseId = idx + caseBaseNo < numDPInstances ? idx + caseBaseNo : 0;
            uint dnaTPARA = ( caseId >> 5 ) * ( MC_CeilDivide16 ( maxDNALength ) << 5 ) + ( caseId & 0x1F );
            uint readTPARA = ( caseId >> 5 ) * ( MC_CeilDivide16 ( maxReadLength ) << 5 ) + ( caseId & 0x1F );
            for ( int j=1;j<=maxDNALength;++j )
                { *dna_ptr = j>DNALengths[caseId]?4:MC_DnaUnpack ( DNASequences, j ); dna_ptr += N_DP_PACKED; }

            for ( int j=1;j<=maxReadLength;++j )
                { *read_ptr = j>readLengths[caseId]?8:MC_ReadUnpack ( ReadSequences, j ); read_ptr += N_DP_PACKED; }
        }

        // need optimize?
        dummy_pack(bufferHi, bufferLo, DNALengths, caseBaseNo);
        _DNALengthsHi = MC_SetEpi8 ( bufferHi, 0 );
        _DNALengthsLo = MC_SetEpi8 ( bufferLo, 0 );
    

        dummy_pack(bufferHi, bufferLo, readLengths, caseBaseNo);
        _readLengthsHi = MC_SetEpi8 ( bufferHi, 0 );
        _readLengthsLo = MC_SetEpi8 ( bufferLo, 0 );

        _cutOffThresholds = MC_SetEpi8 ( cutOffThresholds, caseBaseNo );

        SemiGlobalAlignment ( _DNASequence, _DNALengthsHi, _DNALengthsLo, maxDNALength,
                              _readSequence, _readLengthsHi, _readLengthsLo,  maxReadLength,
                              _maxScoresHi,_maxScoresLo,  _hitLocsHi, _hitLocsLo,
                              clipLtSizes, clipRtSizes,
                              (int8_t)MismatchScore,
                              (int8_t)GapOpenScore,
                              DPTable, DPTempRow, _maxScoreCounts,
                              _cutOffThresholds, pattern + caseBaseNo * ( maxDNALength + maxReadLength ), 
                              numCase);

        for ( int j=0;j<N_DP_PACKED;++j )
        {
//int z = caseBaseNo+j;
            scores[caseBaseNo+j] = extract_epi8 ( _maxScoresHi, j ) * 256 + extract_epi8(_maxScoresLo, j);
            hitLocs[caseBaseNo+j] = extract_epi8 ( _hitLocsHi, j ) * 256 + extract_epi8 ( _hitLocsLo, j ) ;
            maxScoreCounts[caseBaseNo+j] = extract_epi8 ( _maxScoreCounts, j );
//if (READS[j].find(FOO) != -1 || READS[j].find(FOO_REV) != -1) if (z < numDPInstances ) fprintf(stderr, "READS %s\nREF%s\n%d %d %d %d\n",READS[j].c_str(), REFS[j].c_str(), z, scores[z], hitLocs[z], maxScoreCounts[z]);
//if (READS[j].find(BAR) != -1 || READS[j].find(BAR_REV) != -1) if (z < numDPInstances ) fprintf(stderr, "READS %s\nREF%s\n%d %d %d %d\n",READS[j].c_str(), REFS[j].c_str(), z, scores[z], hitLocs[z], maxScoreCounts[z]);
        }
    }
}

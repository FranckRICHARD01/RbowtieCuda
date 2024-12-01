/*
 * nvbio
 * Copyright (c) 2011-2014, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the NVIDIA CORPORATION nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// alignment_test.cu
//

#pragma once

#include <nvbio-test/alignment_test_utils.h>
#include <nvbio/basic/timer.h>
#include <nvbio/basic/console.h>
#include <nvbio/basic/cuda/ldg.h>
#include <nvbio/basic/cached_iterator.h>
#include <nvbio/basic/packedstream.h>
#include <nvbio/basic/packedstream_loader.h>
#include <nvbio/basic/vector_view.h>
#include <nvbio/basic/shared_pointer.h>
#include <nvbio/basic/dna.h>
#include <nvbio/alignment/alignment.h>
#include <nvbio/alignment/batched.h>
#include <nvbio/alignment/sink.h>
#include <thrust/device_vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

namespace nvbio
{
    namespace aln
    {

        enum AlignmentTest
        {
            ALL = 0xFFFFFFFFu,
            ED = 1u,
            SW = 2u,
            GOTOH = 4u,
            ED_BANDED = 8u,
            SW_BANDED = 16u,
            GOTOH_BANDED = 32u,
            SW_WARP = 64u,
            SW_STRIPED = 128u,
            FUNCTIONAL = 256u,
            WFA = 512u,
            WFA_BANDED = 1024u
        };

        // make a light-weight string from an ASCII char string
        vector_view<const char *> make_string(const char *str)
        {
            return vector_view<const char *>(uint32(strlen(str)), str);
        }

        // run-length encode a string
        std::string rle(const char *input)
        {
            char buffer[1024];
            buffer[0] = '\0';

            char prev = '\0';
            uint32 cnt = 0;

            for (const char *p = input; *p != '\0'; ++p)
            {
                if (*p == prev)
                    ++cnt;
                else
                {
                    if (p != input)
                        sprintf(buffer + strlen(buffer), "%u%c", cnt, prev);
                    prev = *p;
                    cnt = 1u;
                }
            }
            if (cnt)
                sprintf(buffer + strlen(buffer), "%u%c", cnt, prev);
            return std::string(buffer);
        }

        // fill an array with random symbols
        template <uint32 BITS, typename rand_type>
        void fill_packed_stream(rand_type &rand, const uint32 value_range, const uint32 n, uint32 *storage)
        {
            typedef PackedStream<uint32 *, uint8, BITS, false> packed_stream_type;

            packed_stream_type stream(storage);

            for (uint32 i = 0; i < n; ++i)
                stream[i] = rand.next() % value_range;
        }

        // abstract container for the score matrices
        template <uint32 N, uint32 M, typename aligner_tag>
        struct ScoreMatrices
        {
        };

        // ScoreMatrices EditDistanceTag-specialization
        template <uint32 N, uint32 M>
        struct ScoreMatrices<N, M, aln::EditDistanceTag>
        {
            int32 H[N + 1][M + 1];
            char H_flow[N + 1][M + 1];

            void print()
            {
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "H[%2u]  : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, " %3d", H[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "Hd[%2u] : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, "   %c", H_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
            }
        };

        // ScoreMatrices SmithWatermanTag-specialization
        template <uint32 N, uint32 M>
        struct ScoreMatrices<N, M, aln::SmithWatermanTag>
        {
            int32 H[N + 1][M + 1];
            char H_flow[N + 1][M + 1];

            void print()
            {
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "H[%2u]  : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, " %3d", H[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "Hd[%2u] : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, "   %c", H_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
            }
        };

        // ScoreMatrices GotohTag-specialization
        template <uint32 N, uint32 M>
        struct ScoreMatrices<N, M, aln::GotohTag>
        {
            int32 H[N + 1][M + 1];
            int32 E[N + 1][M + 1];
            int32 F[N + 1][M + 1];
            char H_flow[N + 1][M + 1];
            char E_flow[N + 1][M + 1];
            char F_flow[N + 1][M + 1];

            void print()
            {
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "H[%2u]  : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, " %3d", H[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "E[%2u]  : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, " %3d", E[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "F[%2u]  : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, " %3d", F[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "Hd[%2u] : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, "   %c", H_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "Ed[%2u] : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, "   %c", E_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= N; ++i)
                {
                    fprintf(stderr, "Fd[%2u] : ", i);
                    for (uint32 j = 0; j <= M; ++j)
                        fprintf(stderr, "   %c", F_flow[i][j]);
                    fprintf(stderr, "\n");
                }
            }
        };

        // ScoreMatrices WfahTag-specialization
        template <uint32 N, uint32 M>
        struct ScoreMatrices<N, M, aln::WfahTag>
        {
            int32 Score_H_band[N][M];
            int32 Score_E_band[N][M];
            int32 Score_F_band[N][M];
            int32 Point_H_BAND[N];
            int32 HH[N][M];
            int32 EE[N][M];
            int32 FF[N][M];
            int32 Max_H_band[N];
            int32 Min_H_band[N];
            int32 Max_F_band[N];
            int32 Min_F_band[N];
            int32 Max_E_band[N];
            int32 Min_E_band[N];
            char H_flow[N][M];
            char E_flow[N][M];
            char F_flow[N][M];
            char HH_flow[N][M];
            char EE_flow[N][M];
            char FF_flow[N][M];

            void print()
            {
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "H[%2u]  : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, " %3d", (HH[i][j] < -1000) ? -1 : HH[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "E[%2u]  : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, " %3d", (EE[i][j] < -1000) ? -1 : EE[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "F[%2u]  : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, " %3d", (FF[i][j] < -1000) ? -1 : FF[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "Hd[%2u] : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, "   %c", HH_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "Ed[%2u] : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, "   %c", EE_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                for (uint32 i = 0; i <= M / 20; i++)
                {
                    fprintf(stderr, "Fd[%2u] : ", i);
                    for (uint32 j = 0; j <= N / 20; j++)
                        fprintf(stderr, "   %c", FF_flow[i][j]);
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
            }
        };

        template <uint32 M, uint32 N, uint32 BAND_LEN, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_banded_sw(const uint8 *str, const uint8 *ref, const uint32 pos, const aln::SmithWatermanAligner<TYPE, scheme_type> aligner)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G = scoring.deletion();
            const int32 I = scoring.insertion();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            const uint32 start = pos;

            int32 best_score = Field_traits<int32>::min();

            int32 band[BAND_LEN];
            for (int32 j = 0; j < BAND_LEN; ++j)
                band[j] = 0;

            for (int32 i = 0; i < M; ++i)
            {
                const uint8 q = str[i];

                int32 top, left, diagonal, hi;

                // j == 0 case
                {
                    const int32 S_ij = (ref[start + i] == q) ? V : S;
                    diagonal = band[0] + S_ij;
                    top = band[1] + G;
                    hi = nvbio::max(top, diagonal);
                    if (TYPE == aln::LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0));           // clamp to zero
                        best_score = nvbio::max(best_score, hi); // save the highest score in the matrix
                    }
                    band[0] = hi;
                }

                for (uint32 j = 1; j < BAND_LEN - 1; ++j)
                {
                    const int32 S_ij = (ref[start + i + j] == q) ? V : S;
                    diagonal = band[j] + S_ij;
                    top = band[j + 1] + G;
                    left = band[j - 1] + I;
                    hi = nvbio::max3(top, left, diagonal);
                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0));           // clamp to zero
                        best_score = nvbio::max(best_score, hi); // save the highest score in the matrix
                    }
                    band[j] = hi;
                }

                // j == BAND_LEN-1 case
                const int32 S_ij = (ref[start + i + BAND_LEN - 1] == q) ? V : S;
                diagonal = band[BAND_LEN - 1] + S_ij;
                left = band[BAND_LEN - 2] + I;
                hi = nvbio::max(left, diagonal);
                if (TYPE == LOCAL)
                {
                    hi = nvbio::max(hi, int32(0));           // clamp to zero
                    best_score = nvbio::max(best_score, hi); // save the highest score in the matrix
                }
                band[BAND_LEN - 1] = hi;
            }

            if (TYPE == GLOBAL)
                best_score = band[BAND_LEN - 1];
            else if (TYPE == SEMI_GLOBAL)
            {
                // get the highest score along the long edge of the path graph
                best_score = band[0];
                for (uint32 j = 1; j < BAND_LEN; ++j)
                    best_score = nvbio::max(best_score, band[j]);
            }
            return best_score;
        }

        template <uint32 M, uint32 N, uint32 BAND_LEN, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_banded_sw(const uint8 *pattern, const uint8 *text, const uint32 pos, const aln::GotohAligner<TYPE, scheme_type> aligner)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G_o = scoring.pattern_gap_open();
            const int32 G_e = scoring.pattern_gap_extension();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            const uint32 start = pos;

            int32 best_score = Field_traits<int32>::min();
            int32 infimum = Field_traits<int32>::min() - G_e;

            int32 H_band[BAND_LEN];
            int32 F_band[BAND_LEN];

            H_band[0] = 0;
            for (uint32 j = 1; j < BAND_LEN; ++j)
                H_band[j] = TYPE == GLOBAL ? scoring.text_gap_open() + (j - 1) * scoring.text_gap_extension() : 0;

            for (uint32 j = 0; j < BAND_LEN; ++j)
                F_band[j] = infimum;

#define DEBUG_THIS
#if defined(DEBUG_THIS)
            FILE *file = fopen("debug.txt", "w");
#define DEBUG_THIS_STATEMENT(x) x
#else
#define DEBUG_THIS_STATEMENT(x)
#endif

            for (uint32 i = 0; i < M; ++i)
            {
#if defined(DEBUG_THIS)
                {
                    fprintf(file, "F[%2u] = ", i);
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        fprintf(file, " %3d", F_band[j]);
                    fprintf(file, "\n");

                    fprintf(file, "H[%2u] = ", i);
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        fprintf(file, " %3d", H_band[j]);
                    fprintf(file, "\n\n");
                }
#endif

                const uint8 q = pattern[i];

                DEBUG_THIS_STATEMENT(fprintf(file, "U[%2u] = ", i));

                // update F for the next row
                for (uint32 j = 0; j < BAND_LEN - 1; ++j)
                {
                    const int32 ftop = F_band[j + 1] + G_e;
                    const int32 htop = H_band[j + 1] + G_o;
                    F_band[j] = nvbio::max(ftop, htop);
                    if (i && (j == BAND_LEN - 2))
                        assert(F_band[j] == htop);

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", ftop > htop ? 'F' : 'H'));
                }
                F_band[BAND_LEN - 1] = infimum;

                // if (i < 8) fprintf(stderr,"\nU[%2u] = ", i);
                //  j == 0 case
                {
                    const int32 S_ij = (text[start + i] == q) ? V : S;
                    const int32 diagonal = H_band[0] + S_ij;
                    const int32 top = F_band[0];
                    int32 hi = nvbio::max(top, diagonal);

                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[0] = hi;

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", top > diagonal ? 'I' : 'S'));
                }

                // compute E_1
                int32 E_j = H_band[0] + G_o;

                for (uint32 j = 1; j < BAND_LEN - 1; ++j)
                {
                    const uint32 g = text[start + i + j];

                    const int32 S_ij = (g == q) ? V : S;
                    const int32 diagonal = H_band[j] + S_ij;
                    const int32 top = F_band[j];
                    const int32 left = E_j;
                    int32 hi = nvbio::max3(top, left, diagonal);
                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[j] = hi;
                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", (top > left ? (top > diagonal ? 'I' : 'S') : (left > diagonal ? 'D' : 'S'))));

                    // update E for the next round, i.e. j+1
                    const int32 eleft = E_j + G_e;
                    const int32 ediagonal = hi + G_o;
                    E_j = nvbio::max(ediagonal, eleft);
                }

                // j == BAND_LEN-1 case
                {
                    const uint8 g = text[start + i + BAND_LEN - 1];

                    const int32 S_ij = (g == q) ? V : S;
                    const int32 diagonal = H_band[BAND_LEN - 1] + S_ij;
                    const int32 left = E_j;
                    int32 hi = nvbio::max(left, diagonal);
                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[BAND_LEN - 1] = hi;

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", left > diagonal ? 'D' : 'S'));
                }
                DEBUG_THIS_STATEMENT(fprintf(file, "\n"));
            }

            DEBUG_THIS_STATEMENT(fclose(file));

            if (TYPE == GLOBAL)
                best_score = H_band[BAND_LEN - 1];
            else if (TYPE == SEMI_GLOBAL)
            {
                // get the highest score along the long edge of the path graph
                best_score = H_band[0];
                for (uint32 j = 1; j < BAND_LEN; ++j)
                    best_score = nvbio::max(best_score, H_band[j]);
            }
            return best_score;

#if defined(DEBUG_THIS)
#undef DEBUG_THIS_STATEMENT
#undef DEBUG_THIS
#endif
        }

        template <uint32 M, uint32 N, uint32 BAND_LEN, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_banded_sw(const uint8 *pattern, const uint8 *text, const uint32 pos, const aln::WfahAligner<TYPE, scheme_type> aligner)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G_o = scoring.pattern_gap_open();
            const int32 G_e = scoring.pattern_gap_extension();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            const uint32 start = pos;

            int32 best_score = Field_traits<int32>::min();
            int32 infimum = Field_traits<int32>::min() - G_e;

            int32 H_band[BAND_LEN];
            int32 F_band[BAND_LEN];

            H_band[0] = 0;
            for (uint32 j = 1; j < BAND_LEN; ++j)
                H_band[j] = TYPE == GLOBAL ? scoring.text_gap_open() + (j - 1) * scoring.text_gap_extension() : 0;

            for (uint32 j = 0; j < BAND_LEN; ++j)
                F_band[j] = infimum;

#define DEBUG_THIS
#if defined(DEBUG_THIS)
            FILE *file = fopen("debug.txt", "w");
#define DEBUG_THIS_STATEMENT(x) x
#else
#define DEBUG_THIS_STATEMENT(x)
#endif

            for (uint32 i = 0; i < M; ++i)
            {
#if defined(DEBUG_THIS)
                {
                    fprintf(file, "F[%2u] = ", i);
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        fprintf(file, " %3d", F_band[j]);
                    fprintf(file, "\n");

                    fprintf(file, "H[%2u] = ", i);
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        fprintf(file, " %3d", H_band[j]);
                    fprintf(file, "\n\n");
                }
#endif

                const uint8 q = pattern[i];

                DEBUG_THIS_STATEMENT(fprintf(file, "U[%2u] = ", i));

                // update F for the next row
                for (uint32 j = 0; j < BAND_LEN - 1; ++j)
                {
                    const int32 ftop = F_band[j + 1] + G_e;
                    const int32 htop = H_band[j + 1] + G_o;
                    F_band[j] = nvbio::max(ftop, htop);
                    if (i && (j == BAND_LEN - 2))
                        assert(F_band[j] == htop);

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", ftop > htop ? 'F' : 'H'));
                }
                F_band[BAND_LEN - 1] = infimum;

                // if (i < 8) fprintf(stderr,"\nU[%2u] = ", i);
                //  j == 0 case
                {
                    const int32 S_ij = (text[start + i] == q) ? V : S;
                    const int32 diagonal = H_band[0] + S_ij;
                    const int32 top = F_band[0];
                    int32 hi = nvbio::max(top, diagonal);

                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[0] = hi;

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", top > diagonal ? 'I' : 'S'));
                }

                // compute E_1
                int32 E_j = H_band[0] + G_o;

                for (uint32 j = 1; j < BAND_LEN - 1; ++j)
                {
                    const uint32 g = text[start + i + j];

                    const int32 S_ij = (g == q) ? V : S;
                    const int32 diagonal = H_band[j] + S_ij;
                    const int32 top = F_band[j];
                    const int32 left = E_j;
                    int32 hi = nvbio::max3(top, left, diagonal);
                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[j] = hi;
                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", (top > left ? (top > diagonal ? 'I' : 'S') : (left > diagonal ? 'D' : 'S'))));

                    // update E for the next round, i.e. j+1
                    const int32 eleft = E_j + G_e;
                    const int32 ediagonal = hi + G_o;
                    E_j = nvbio::max(ediagonal, eleft);
                }

                // j == BAND_LEN-1 case
                {
                    const uint8 g = text[start + i + BAND_LEN - 1];

                    const int32 S_ij = (g == q) ? V : S;
                    const int32 diagonal = H_band[BAND_LEN - 1] + S_ij;
                    const int32 left = E_j;
                    int32 hi = nvbio::max(left, diagonal);
                    if (TYPE == LOCAL)
                    {
                        hi = nvbio::max(hi, int32(0)); // clamp to zero
                        best_score = nvbio::max(best_score, hi);
                    }
                    H_band[BAND_LEN - 1] = hi;

                    DEBUG_THIS_STATEMENT(fprintf(file, "   %c", left > diagonal ? 'D' : 'S'));
                }
                DEBUG_THIS_STATEMENT(fprintf(file, "\n"));
            }

            DEBUG_THIS_STATEMENT(fclose(file));

            if (TYPE == GLOBAL)
                best_score = H_band[BAND_LEN - 1];
            else if (TYPE == SEMI_GLOBAL)
            {
                // get the highest score along the long edge of the path graph
                best_score = H_band[0];
                for (uint32 j = 1; j < BAND_LEN; ++j)
                    best_score = nvbio::max(best_score, H_band[j]);
            }
            return best_score;

#if defined(DEBUG_THIS)
#undef DEBUG_THIS_STATEMENT
#undef DEBUG_THIS
#endif
        }

        template <uint32 M, uint32 N, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_sw(
            const uint8 *str,
            const uint8 *ref,
            const uint32 length_str,
            const uint32 length_ref,
            const aln::SmithWatermanAligner<TYPE, scheme_type> aligner,
            ScoreMatrices<M, N, aln::SmithWatermanTag> *mat)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G = scoring.deletion();
            const int32 I = scoring.insertion();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            mat->H[0][0] = 0;
            mat->H_flow[0][0] = '*';
            for (uint32 j = 1; j <= M; ++j)
            {
                mat->H[0][j] = TYPE != LOCAL ? G * j : 0;
                mat->H_flow[0][j] = '*';
            }
            for (uint32 i = 1; i <= N; ++i)
            {
                mat->H[i][0] = TYPE == GLOBAL ? G * i : 0;
                mat->H_flow[i][0] = '*';
            }

            int32 best_score = Field_traits<int32>::min();

            for (uint32 i = 1; i <= N; ++i)
            {
                for (uint32 j = 1; j <= M; ++j)
                {
                    const uint8 r_i = ref[i - 1];
                    const uint8 s_j = str[j - 1];
                    const int32 S_ij = (r_i == s_j) ? V : S;

                    const int32 top = mat->H[i - 1][j] + G;
                    const int32 left = mat->H[i][j - 1] + I;
                    const int32 diag = mat->H[i - 1][j - 1] + S_ij;

                    mat->H[i][j] = nvbio::max3(top, left, diag);

                    mat->H_flow[i][j] = top > left ? (top > diag ? '|' : '\\') : (left > diag ? '-' : '\\');

                    if (TYPE == aln::LOCAL)
                        best_score = nvbio::max(best_score, mat->H[i][j] = nvbio::max(mat->H[i][j], int32(0)));
                }
                if (TYPE == aln::SEMI_GLOBAL)
                    best_score = nvbio::max(best_score, mat->H[i][M]);
            }
            if (TYPE == aln::LOCAL || TYPE == aln::SEMI_GLOBAL)
                return best_score;
            else
                return mat->H[N][M];
        }

        template <uint32 M, uint32 N, aln::AlignmentType TYPE>
        int32 ref_sw(
            const uint8 *str,
            const uint8 *ref,
            const uint32 length_str,
            const uint32 length_ref,
            const aln::EditDistanceAligner<TYPE> aligner,
            ScoreMatrices<N, M, aln::EditDistanceTag> *mat)
        {
            return ref_sw<M, N>(
                str,
                ref,
                length_str,
                length_ref,
                aln::make_smith_waterman_aligner<TYPE>(priv::EditDistanceSWScheme()),
                (ScoreMatrices<N, M, aln::SmithWatermanTag> *)mat);
        }

        template <uint32 M, uint32 N, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_sw(
            const uint8 *str,
            const uint8 *ref,
            const uint32 length_str,
            const uint32 length_ref,
            const aln::GotohAligner<TYPE, scheme_type> aligner,
            ScoreMatrices<N, M, aln::GotohTag> *mat)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G_o = scoring.pattern_gap_open();
            const int32 G_e = scoring.pattern_gap_extension();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            mat->H[0][0] = 0;
            mat->F[0][0] = (TYPE != LOCAL) ? -100000 : 0;
            mat->E[0][0] = (TYPE != LOCAL) ? -100000 : 0;
            mat->H_flow[0][0] = '*';
            mat->E_flow[0][0] = '*';
            mat->F_flow[0][0] = '*';
            for (uint32 j = 1; j <= M; ++j)
            {
                mat->H[0][j] = (TYPE != LOCAL) ? G_o + G_e * (j - 1) : 0;
                mat->E[0][j] = mat->F[0][j] = (TYPE != LOCAL) ? -100000 : 0;
                mat->H_flow[0][j] = '*';
                mat->E_flow[0][j] = '*';
                mat->F_flow[0][j] = '*';
            }
            for (uint32 i = 1; i <= N; ++i)
            {
                mat->H[i][0] = (TYPE == aln::GLOBAL) ? G_o + G_e * (i - 1) : 0;
                mat->E[i][0] = mat->F[i][0] = (TYPE != LOCAL) ? -100000 : 0;
                mat->H_flow[i][0] = '*';
                mat->E_flow[i][0] = '*';
                mat->F_flow[i][0] = '*';
            }

            int32 best_score = Field_traits<int32>::min();
            uint2 best_sink = make_uint2(uint32(-1), uint32(-1));

            for (uint32 i = 1; i <= N; ++i)
            {
                for (uint32 j = 1; j <= M; ++j)
                {
                    const uint8 r_i = ref[i - 1];
                    const uint8 s_j = str[j - 1];
                    const int32 S_ij = (r_i == s_j) ? V : S;

                    mat->E_flow[i][j] = mat->E[i][j - 1] + G_e > mat->H[i][j - 1] + G_o ? '-' : 'H';
                    mat->F_flow[i][j] = mat->F[i - 1][j] + G_e > mat->H[i - 1][j] + G_o ? '|' : 'H';

                    mat->E[i][j] = nvbio::max(mat->E[i][j - 1] + G_e, mat->H[i][j - 1] + G_o);
                    mat->F[i][j] = nvbio::max(mat->F[i - 1][j] + G_e, mat->H[i - 1][j] + G_o);

                    mat->H_flow[i][j] = mat->F[i][j] > mat->E[i][j] ? (mat->F[i][j] > mat->H[i - 1][j - 1] + S_ij ? 'F' : '\\') : (mat->E[i][j] > mat->H[i - 1][j - 1] + S_ij ? 'E' : '\\');

                    mat->H[i][j] = nvbio::max3(mat->H[i - 1][j - 1] + S_ij, mat->E[i][j], mat->F[i][j]);

                    if (TYPE == aln::LOCAL)
                    {
                        mat->H[i][j] = nvbio::max(mat->H[i][j], int32(0));

                        if (mat->H[i][j] == 0)
                            mat->H_flow[i][j] = '0';

                        if (best_score < mat->H[i][j])
                        {
                            best_score = mat->H[i][j];
                            best_sink = make_uint2(i, j);
                        }
                    }
                }
                if (TYPE == aln::SEMI_GLOBAL)
                {
                    if (best_score < mat->H[i][M])
                    {
                        best_score = mat->H[i][M];
                        best_sink = make_uint2(i, M);
                    }
                }
            }
            if (TYPE == aln::LOCAL || TYPE == aln::SEMI_GLOBAL)
                return best_score;
            else
                return mat->H[N][M];
        }

        // wavefront_compute.c (576)
        template <int32 N, int32 M>
        void trim_ends(
            int32 mat[][M],
            int32 *max,
            int32 *min,
            int32 s,
            int32 deltaH,
            int32 deltaV,
            int32 length_ref,
            int32 length_str)
        {
            int32 i, k;

            for (k = max[deltaV + s]; k >= min[deltaV + s]; --k)
            {
                int32 offset = mat[deltaV + s][deltaH + k];
                uint32 v = (uint32)(offset - k);
                uint32 h = (uint32)offset;

                if (h <= (uint32)length_str && v <= (uint32)length_ref)
                    break; // ?
            }

            for (i = k + 1; i <= max[deltaV + s]; i++)
                mat[deltaV + s][deltaH + i] = WFA_MIN + 1;

            max[deltaV + s] = k;

            for (k = min[deltaV + s]; k <= max[deltaV + s]; ++k)
            {
                int32 offset = mat[deltaV + s][deltaH + k];
                uint32 v = (uint32)(offset - k);
                uint32 h = (uint32)offset;

                if (h <= (uint32)length_str && v <= (uint32)length_ref)
                    break; // ?
            }

            for (i = min[deltaV + s]; i < k - 1; i++)
                mat[deltaV + s][deltaH + i] = WFA_MIN + 1;

            min[deltaV + s] = k;

            if (min[deltaV + s] > max[deltaV + s])
                mat[deltaV + s][deltaH] = WFA_MIN;
        }

       int32 f(
            int32 i,
            int32 deltaV,
            int32 *mat)
        {
            if (i == 0)
                return 0;

            int32 k = 0;

            while (mat[deltaV + k] < i)
                k++;

            if (mat[deltaV + k] == i)
                return k + 1;

            return -1;
        }
        
        template <int32 N, int32 M>
        bool testEnd(
            int32 length_ref,
            int32 length_str,
            int32 mat[][M],
            int32 *max,
            int32 *min,
            int32 deltaVss,
            int32 deltaH,
            int32 *alignemnt_k,
            int32 *alignemnt_offset)
        {
            *alignemnt_k = length_str - length_ref;
            *alignemnt_offset = length_str;

            int32 min_lo;
            int32 max_hi;

            min_lo = min[deltaVss];
            max_hi = max[deltaVss];

            if (min_lo > max_hi)
                return false;

            if (min_lo > *alignemnt_k || *alignemnt_k > max_hi)
                return false;

            int32 moffset = mat[deltaVss][deltaH];
            if (moffset < *alignemnt_offset)
                return false;

            return true;
        }

        template <int32 N, int32 M, aln::AlignmentType TYPE, typename scheme_type>
        int32 ref_sw(
            const uint8 *str,
            const uint8 *ref,
            const int32 length_str,
            const int32 length_ref,
            const aln::WfahAligner<TYPE, scheme_type> aligner,
            ScoreMatrices<N, M, aln::WfahTag> *mat)
        {
            const scheme_type &scoring = aligner.scheme;

            const int32 G_o = scoring.pattern_gap_open();
            const int32 G_e = scoring.pattern_gap_extension();
            const int32 S = scoring.mismatch();
            const int32 V = scoring.match();

            int32 min_lo;
            int32 max_hi;
            int32 s = 0;
            int32 ss = 0;
            int32 last_s = 0;
            int32 indice = 0;
            int32 mismatch = 0;
            int32 gap_open1 = 0;
            int32 gap_extend1 = 0;

            int32 deltaH = M / 2;
            int32 deltaV = 10;

            int32 h_save[N][M];

            bool debug = false;

            for (uint32 i = 0; i < N; i++)
            {
                for (uint32 j = 0; j < M; j++)
                {
                    mat->Score_H_band[i][j] = mat->Score_E_band[i][j] = mat->Score_F_band[i][j] = WFA_MIN;
                    mat->H_flow[i][j] = '\\';
                    mat->E_flow[i][j] = ' ';
                    mat->F_flow[i][j] = ' ';
                    mat->HH_flow[i][j] = '\\';
                    mat->EE_flow[i][j] = '-';
                    mat->FF_flow[i][j] = '|';
                    h_save[i][j] = WFA_MIN;
                }
            }
            mat->Score_H_band[0][0] = 1;
            mat->Score_F_band[0][0] = 0;
            mat->Score_E_band[0][0] = 0;
            mat->H_flow[0][0] = '*';
            mat->E_flow[0][0] = '*';
            mat->F_flow[0][0] = '*';
            for (int32 j = 1; j < M; j++)
            {
                mat->Score_H_band[0][j] = 0;
                mat->Score_E_band[0][j] = 0;
                mat->Score_F_band[0][j] = 0;
                mat->H_flow[0][j] = '*';
                mat->E_flow[0][j] = '*';
                mat->F_flow[0][j] = '*';
                mat->HH_flow[0][j] = '*';
                mat->EE_flow[0][j] = '*';
                mat->FF_flow[0][j] = '*';
            }
            for (int32 i = 1; i < N; i++)
            {
                mat->Max_H_band[i] = (i < deltaV) ? -1 : 0;
                mat->Min_H_band[i] = (i < deltaV) ? 1 : 0;
                mat->Max_F_band[i] = -1;
                mat->Min_F_band[i] = 1;
                mat->Max_E_band[i] = -1;
                mat->Min_E_band[i] = 1;
                mat->H_flow[i][0] = '*';
                mat->E_flow[i][0] = '*';
                mat->F_flow[i][0] = '*';
                mat->HH_flow[i][0] = '*';
                mat->EE_flow[i][0] = '*';
                mat->FF_flow[i][0] = '*';
                mat->Point_H_BAND[i] = 0;
            }

            int32 Score_H_band_last = WFA_MIN;
            int32 Score_E_band_last = WFA_MIN;
            int32 Score_F_band_last = WFA_MIN;

            int32 alignement_k;
            int32 alignement_offset;

            // int32 best_score = Field_traits<int32>::min();
            // uint2 best_sink  = make_uint2(uint32(-1),uint32(-1));

            mat->Score_H_band[deltaV][deltaH] = 0;

            do
            {
                if (ss >= 0 && mat->Score_H_band[deltaV + ss][deltaH] != WFA_MIN)
                {
                    min_lo = mat->Min_H_band[deltaV + ss];
                    max_hi = mat->Max_H_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_H_band[deltaV + ss][deltaH + k];

                        int32 h = offset - k;
                        int32 v = offset;

                        while (h >= 0 && v >= 0 && h < length_ref && v < length_str && ref[h] == str[v])
                        {
                            mat->Score_H_band[deltaV + ss][deltaH + k]++;
                            h++;
                            v++;
                            if (debug)
                                fprintf(stderr, "* v=%d h=%d mat->H[%d][%d]=%d\n", v - 1, h - 1, ss, k, mat->Score_H_band[deltaV + ss][deltaH + k]);
                        }

                        char hdir = '\\';

                        // if (debug) fprintf(stderr, "Calcul Hdir[%d][%d]: ", ss, k);

                        if (mat->Score_F_band[deltaV + ss][deltaH + k] > mat->Score_E_band[deltaV + ss][deltaH + k])
                        {
                            // if (debug) fprintf(stderr, "F(%d) > E(%d), ", mat->Score_F_band[deltaV + ss][deltaH + k], mat->Score_E_band[deltaV + ss][deltaH + k]);
                            if (mat->Score_F_band[deltaV + ss][deltaH + k] != WFA_MIN && h_save[deltaV + ss][deltaH + k] != WFA_MIN)
                            {
                                if (mat->Score_F_band[deltaV + ss][deltaH + k] > h_save[deltaV + ss][deltaH + k])
                                {
                                    // if (debug) fprintf(stderr, "F(%d) > H(%d)\n", mat->Score_F_band[deltaV + ss][deltaH + k], h_save[deltaV + ss][deltaH + k]);
                                    hdir = 'F';
                                }
                                // else
                                // if (debug) fprintf(stderr, "\n");
                            }
                        }
                        else
                        {
                            // if (debug) fprintf(stderr, "E(%d) > F(%d), ", mat->Score_E_band[deltaV + ss][deltaH + k], mat->Score_F_band[deltaV + ss][deltaH + k]);
                            if (mat->Score_E_band[deltaV + ss][deltaH + k] != WFA_MIN && h_save[deltaV + ss][deltaH + k] != WFA_MIN)
                            {
                                if (mat->Score_E_band[deltaV + ss][deltaH + k] > h_save[deltaV + ss][deltaH + k])
                                {
                                    // if (debug) fprintf(stderr, "E(%d) > H(%d)\n", mat->Score_E_band[deltaV + ss][deltaH + k], h_save[deltaV + ss][deltaH + k]);
                                    hdir = 'E';
                                }
                                // else
                                // if (debug) fprintf(stderr, "\n");
                            }
                        }

                        mat->H_flow[deltaV + ss][deltaH + k] = hdir;
                    }

                    trim_ends<N, M>(mat->Score_H_band, mat->Max_H_band, mat->Min_H_band, ss, deltaH, deltaV, length_ref, length_str);
                    trim_ends<N, M>(mat->Score_F_band, mat->Max_F_band, mat->Min_F_band, ss, deltaH, deltaV, length_ref, length_str);
                    trim_ends<N, M>(mat->Score_E_band, mat->Max_E_band, mat->Min_E_band, ss, deltaH, deltaV, length_ref, length_str);

                    if (debug)
                    {
                        fprintf(stderr, "\n");
                        fprintf(stderr, "  H[%d] lo=%d hi=%d\n", ss, mat->Min_H_band[deltaV + ss], mat->Max_H_band[deltaV + ss]);
                        fprintf(stderr, "  F[%d] lo=%d hi=%d\n", ss, mat->Min_F_band[deltaV + ss], mat->Max_F_band[deltaV + ss]);
                        fprintf(stderr, "  E[%d] lo=%d hi=%d\n", ss, mat->Min_E_band[deltaV + ss], mat->Max_E_band[deltaV + ss]);

                        fprintf(stderr, "final  s=%d:\n", s);
                    }
                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        if (debug)
                            fprintf(stderr, "k = %3d : H = %5d\tE = %5d\tF = %5d\n",
                                    k,
                                    mat->Score_H_band[deltaV + ss][deltaH + k],
                                    mat->Score_E_band[deltaV + ss][deltaH + k],
                                    mat->Score_F_band[deltaV + ss][deltaH + k]);
                    }
                }

                s++;

                // wavefront_terminations.c (16)
                if (mat->Min_H_band[deltaV + ss] <= length_str - length_ref &&
                    length_str - length_ref <= mat->Max_H_band[deltaV + ss] &&
                    mat->Score_H_band[deltaV + ss][deltaH + length_str - length_ref] >= length_str)
                    break;

                if (s > 10000)
                    break;

                if (testEnd<N, M>(N, M, mat->Score_H_band, mat->Max_H_band, mat->Min_H_band, deltaV + s, deltaH, &alignement_k, &alignement_offset))
                    break;

                if (testEnd<N, M>(N, M, mat->Score_E_band, mat->Max_E_band, mat->Min_E_band, deltaV + s, deltaH, &alignement_k, &alignement_offset))
                    break;

                if (testEnd<N, M>(N, M, mat->Score_F_band, mat->Max_F_band, mat->Min_F_band, deltaV + s, deltaH, &alignement_k, &alignement_offset))
                    break;

                mat->Point_H_BAND[deltaV + indice++] = s;

                ss = f(s, deltaV, mat->Point_H_BAND);

                if (ss >= 0)
                {
                    mismatch = f(s - S, deltaV, mat->Point_H_BAND);
                    gap_open1 = f(s - G_o - G_e, deltaV, mat->Point_H_BAND);
                    gap_extend1 = f(s - G_e, deltaV, mat->Point_H_BAND);

                    bool testH = true;
                    bool testO = true;
                    bool testE = true;
                    bool testF = true;

                    if (mismatch < 0 || mat->Score_H_band[deltaV + mismatch][deltaH] == WFA_MIN)
                        testH = false;
                    if (gap_open1 < 0 || mat->Score_H_band[deltaV + gap_open1][deltaH] == WFA_MIN)
                        testO = false;
                    if (gap_extend1 < 0 || mat->Score_E_band[deltaV + gap_extend1][deltaH] == WFA_MIN)
                        testE = false;
                    if (gap_extend1 < 0 || mat->Score_F_band[deltaV + gap_extend1][deltaH] == WFA_MIN)
                        testF = false;

                    if (testH || testO || testE || testF)
                    {
                        //    limits_inputs

                        min_lo = mat->Min_H_band[deltaV + mismatch];
                        max_hi = mat->Max_H_band[deltaV + mismatch];

                        if (min_lo > mat->Min_H_band[deltaV + gap_open1] - 1)
                            min_lo = mat->Min_H_band[deltaV + gap_open1] - 1;
                        if (max_hi < mat->Max_H_band[deltaV + gap_open1] + 1)
                            max_hi = mat->Max_H_band[deltaV + gap_open1] + 1;

                        if (min_lo > mat->Min_E_band[deltaV + gap_extend1] + 1)
                            min_lo = mat->Min_E_band[deltaV + gap_extend1] + 1;
                        if (max_hi < mat->Max_E_band[deltaV + gap_extend1] + 1)
                            max_hi = mat->Max_E_band[deltaV + gap_extend1] + 1;
                        if (min_lo > mat->Min_F_band[deltaV + gap_extend1] - 1)
                            min_lo = mat->Min_F_band[deltaV + gap_extend1] - 1;
                        if (max_hi < mat->Max_F_band[deltaV + gap_extend1] - 1)
                            max_hi = mat->Max_F_band[deltaV + gap_extend1] - 1;

                        if (debug)
                            fprintf(stderr, "\n\nminmax: s=%d ss=%d mat_lo=%d mat_hi=%d\n", s, ss, min_lo, max_hi);

                        // wavefront_compute.c (445)
                        if (Score_H_band_last != WFA_MIN || Score_E_band_last != WFA_MIN)
                        {
                            mat->Max_E_band[deltaV + ss] = max_hi;
                            mat->Min_E_band[deltaV + ss] = min_lo;
                        }
                        else
                        {
                            mat->Max_E_band[deltaV + ss] = -1;
                            mat->Min_E_band[deltaV + ss] = 1;
                        }
                        if (Score_H_band_last != WFA_MIN || Score_F_band_last != WFA_MIN)
                        {
                            mat->Max_F_band[deltaV + ss] = max_hi;
                            mat->Min_F_band[deltaV + ss] = min_lo;
                        }
                        else
                        {
                            mat->Max_F_band[deltaV + ss] = -1;
                            mat->Min_F_band[deltaV + ss] = 1;
                        }
                        mat->Max_H_band[deltaV + ss] = max_hi;
                        mat->Min_H_band[deltaV + ss] = min_lo;

                        // affine dispatcher

                        for (int32 k = min_lo; k <= max_hi; k++)
                        {
                            if (mat->Score_E_band[deltaV + gap_extend1][deltaH + k] != WFA_MIN && mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1] != WFA_MIN)
                            {
                                mat->E_flow[deltaV + ss][deltaH + k] = mat->Score_E_band[deltaV + gap_extend1][deltaH + k - 1] > mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1] ? '-' : 'H';

                                // if (debug) fprintf(stderr, "E_flow[%d][%d]=%c   %d    %d\n", ss, k, mat->E_flow[deltaV + ss][deltaH + k], mat->Score_E_band[deltaV + gap_extend1][deltaH + k - 1], mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1]);
                            }

                            if (mat->Score_F_band[deltaV + gap_extend1][deltaH + k] != WFA_MIN && mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1] != WFA_MIN)
                            {
                                mat->F_flow[deltaV + ss][deltaH + k] = mat->Score_F_band[deltaV + gap_extend1][deltaH + k + 1] > mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1] ? '|' : 'H';

                                // if (debug) fprintf(stderr, "F_flow[%d][%d]=%c   %d    %d\n", ss, k, mat->F_flow[deltaV + ss][deltaH + k], mat->Score_F_band[deltaV + gap_extend1][deltaH + k + 1], mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1]);
                            }

                            int32 E = nvbio::max(mat->Score_E_band[deltaV + gap_extend1][deltaH + k - 1], mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1]) + 1u;
                            int32 F = nvbio::max(mat->Score_F_band[deltaV + gap_extend1][deltaH + k + 1], mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1]);
                            mat->Score_E_band[deltaV + ss][deltaH + k] = E; // I
                            mat->Score_F_band[deltaV + ss][deltaH + k] = F; // D

                            if (false && debug)
                            {
                                if (mat->Score_E_band[deltaV + gap_extend1][deltaH + k - 1] > mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1])
                                {
                                    fprintf(stderr, "E[%d][%d] = E[%d][%d] + 1 = %d\n", ss, k, gap_extend1, k - 1, mat->Score_E_band[deltaV + gap_extend1][deltaH + k - 1] + 1);
                                }
                                else
                                {
                                    fprintf(stderr, "E[%d][%d] = H[%d][%d] + 1 = %d\n", ss, k, gap_open1, k - 1, mat->Score_H_band[deltaV + gap_open1][deltaH + k - 1] + 1);
                                }

                                if (mat->Score_F_band[deltaV + gap_extend1][deltaH + k + 1] > mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1])
                                {
                                    fprintf(stderr, "F[%d][%d] = F[%d][%d] = %d\n", ss, k, gap_extend1, k + 1, mat->Score_F_band[deltaV + gap_extend1][deltaH + k + 1]);
                                }
                                else
                                {
                                    fprintf(stderr, "F[%d][%d] = H[%d][%d] = %d\n", ss, k, gap_open1, k + 1, mat->Score_H_band[deltaV + gap_open1][deltaH + k + 1]);
                                }
                            }

                            int32 H = nvbio::max3(mat->Score_H_band[deltaV + mismatch][deltaH + k] + 1, mat->Score_E_band[deltaV + ss][deltaH + k], mat->Score_F_band[deltaV + ss][deltaH + k]);

                            if (debug)
                                fprintf(stderr, "H[%d][%d] = max(H[%d][%d] + 1, F[%d][%d], EF[%d][%d] = %d\n", ss, k, mismatch, k, ss, k, ss, k, (int)H);

                            // wavefront_compute_affine.c (80)
                            char c = ' ';

                            if (H > length_str || H - k > length_ref)
                            {
                                H = WFA_MIN; // offset null
                                // fprintf(stderr, "  ******correc***** H[%d][%d]\n", ss, k);
                                c = 'c';
                            }

                            mat->Score_H_band[deltaV + ss][deltaH + k] = H;

                            h_save[deltaV + ss][deltaH + k] = mat->Score_H_band[deltaV + mismatch][deltaH + k] + 1;

                            if (mat->Score_H_band[deltaV + ss][deltaH + k] == WFA_MIN)
                                mat->Score_H_band[deltaV + ss][deltaH + k]++;
                            if (mat->Score_E_band[deltaV + ss][deltaH + k] == WFA_MIN)
                                mat->Score_E_band[deltaV + ss][deltaH + k]++;
                            if (mat->Score_F_band[deltaV + ss][deltaH + k] == WFA_MIN)
                                mat->Score_F_band[deltaV + ss][deltaH + k]++;

                            if (debug)
                                fprintf(stderr, "k = %3d : H = %5d%c\tE = %5d\tF = %5d\n",
                                        k,
                                        mat->Score_H_band[deltaV + ss][deltaH + k],
                                        c,
                                        mat->Score_E_band[deltaV + ss][deltaH + k],
                                        mat->Score_F_band[deltaV + ss][deltaH + k]);

                            NVBIO_CUDA_DEBUG_ASSERT(deltaV + ss < N, "indice problem deltaV + ss >= N, s=%d deltaH=%d, (%d >= %d)\n", s, deltaH, deltaV + ss , N);
                            NVBIO_CUDA_DEBUG_ASSERT(deltaH + k < M, "indice problem deltaH + k >= M, s=%d deltaH=%d, (%d >= %d)\n", s, deltaH, deltaH + k , M);
                            NVBIO_CUDA_DEBUG_ASSERT(deltaH + k >= 0, "indice problem deltaH + k < 0, s=%d deltaH=%d, (%d < 0)\n", s, deltaH, deltaH + k );
                        }

                        if (mat->Score_H_band[deltaV + ss][deltaH] == WFA_MIN)
                            mat->Score_H_band[deltaV + ss][deltaH]++;

                        Score_H_band_last = mat->Score_H_band[deltaV + ss][deltaH];
                        Score_E_band_last = mat->Score_E_band[deltaV + ss][deltaH];
                        Score_F_band_last = mat->Score_F_band[deltaV + ss][deltaH];

                        // trim ends

                        if (debug)
                            fprintf(stderr, "\n");

                        last_s = s;
                    }
                    else
                        indice--;
                }

                if (debug)
                    fprintf(stderr, "\n\n");

            } while (true);

            // mat->Score_H_band[deltaV + indice][0] = 0;

            for (uint32 ii = 0; ii < N; ii++)
            {
                for (uint32 jj = 0; jj < M; jj++)
                {
                    mat->HH[ii][jj] = mat->FF[ii][jj] = mat->EE[ii][jj] = WFA_MIN;
                }

                mat->HH[ii][0] = G_o + G_e * (ii - 1);
                if (ii < M)
                    mat->HH[0][ii] = mat->HH[ii][0];
            }

            mat->HH[0][0] = 0;

            for (int32 col = 1; col <= (int32)length_ref; col++)
            {
                for (int32 ii = 0; ii <= indice; ii++)
                {
                    ss = mat->Point_H_BAND[deltaV + ii - 1];

                    // H_band

                    min_lo = mat->Min_H_band[deltaV + ss];
                    max_hi = mat->Max_H_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_H_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            if (k == 0)
                            {
                                for (int32 l = 1; l <= offset; l++)
                                {
                                    if (l == col && mat->HH[l][col] < 0)
                                    {
                                        mat->HH[l][col] = ss;
                                    }
                                }
                            }
                            else
                            {
                                int32 v = offset - k;
                                int32 h = offset;

                                if (v == col)
                                {
                                    mat->HH[h][col] = ss;
                                }
                            }
                        }
                    }

                    // F_band

                    min_lo = mat->Min_F_band[deltaV + ss];
                    max_hi = mat->Max_F_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_F_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            if (k == 0)
                            {
                                for (int32 l = 1; l <= offset; l++)
                                {
                                    if (l == col && mat->FF[l][col] < 0)
                                    {
                                        mat->FF[l][col] = ss;
                                    }
                                }
                            }
                            else
                            {
                                int32 v = offset - k;
                                int32 h = offset;

                                if (v == col)
                                {
                                    mat->FF[h][col] = ss;
                                }
                            }
                        }
                    }

                    // E_band

                    min_lo = mat->Min_E_band[deltaV + ss];
                    max_hi = mat->Max_E_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_E_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            if (k == 0)
                            {
                                for (int32 l = 1; l <= offset; l++)
                                {
                                    if (l == col && mat->EE[l][col] < 0)
                                    {
                                        mat->EE[l][col] = ss;
                                    }
                                }
                            }
                            else
                            {
                                int32 v = offset - k;
                                int32 h = offset;

                                if (v == col && h >= 0 && h <= 8)
                                {
                                    mat->EE[h][col] = ss;
                                }
                            }
                        }
                    }

                    // HH_flow

                    min_lo = mat->Min_H_band[deltaV + ss];
                    max_hi = mat->Max_H_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_H_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            int32 v = offset - k;
                            int32 h = offset;

                            if (v == col && h >= 1)
                            {
                                mat->HH_flow[h][col] = mat->H_flow[deltaV + ss][deltaH + k];
                            }
                        }
                    }

                    // FF_flow

                    min_lo = mat->Min_F_band[deltaV + ss];
                    max_hi = mat->Max_F_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_F_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            int32 v = offset - k;
                            int32 h = offset;

                            if (v == col && h >= 1)
                            {
                                mat->FF_flow[h][col] = mat->F_flow[deltaV + ss][deltaH + k];
                            }
                        }
                    }

                    // EE_flow

                    min_lo = mat->Min_E_band[deltaV + ss];
                    max_hi = mat->Max_E_band[deltaV + ss];

                    for (int32 k = min_lo; k <= max_hi; k++)
                    {
                        int32 offset = mat->Score_E_band[deltaV + ss][deltaH + k];

                        if (offset > 0)
                        {
                            int32 v = offset - k;
                            int32 h = offset;

                            if (v == col && h >= 1)
                            {
                                mat->EE_flow[h][col] = mat->E_flow[deltaV + ss][deltaH + k];
                            }
                        }
                    }
                }
            }

            return last_s;
        }

        //
        // A test Backtracer class, \see Traceback.
        //
        struct TestBacktracker
        {
            void clear()
            {
                aln[0] = '\0';
                aln_len = 0;
            }

            NVBIO_HOST_DEVICE
            void clip(const uint32 len)
            {
                for (uint32 i = 0; i < len; ++i)
                    aln[aln_len + i] = 'S';
                aln[aln_len + len] = '\0';
            }
            NVBIO_HOST_DEVICE
            void push(const uint8 op)
            {
                aln[aln_len++] = "MID"[op];
            }

            // compute the score of the resulting alignment
            //
            template <AlignmentType TYPE>
            int32 score(
                EditDistanceAligner<TYPE> aligner,
                const uint32 offset,
                const uint8 *str,
                const uint8 *ref)
            {
                int32 score = 0;

                for (uint32 a = 0, j = 0, k = offset; a < aln_len; ++a)
                {
                    const char op = aln[aln_len - a - 1];

                    if (op == 'S')
                        ++j;
                    else if (op == 'I')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first insertion or not
                        --score;

                        // mover in the pattern
                        ++j;
                    }
                    else if (op == 'D')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first deletion or not
                        --score;

                        // move in the text
                        ++k;
                    }
                    else if (op == 'M')
                    {
                        // add a match or mismatch score
                        score -= (str[j] == ref[k]) ? 0 : 1;

                        // move diagonally
                        ++j;
                        ++k;
                    }
                }
                return score;
            }

            // compute the score of the resulting alignment
            //
            template <AlignmentType TYPE, typename scoring_type>
            int32 score(
                SmithWatermanAligner<TYPE, scoring_type> aligner,
                const uint32 offset,
                const uint8 *str,
                const uint8 *ref)
            {
                const scoring_type &scoring = aligner.scheme;

                int32 score = 0;

                for (uint32 a = 0, j = 0, k = offset; a < aln_len; ++a)
                {
                    const char op = aln[aln_len - a - 1];

                    if (op == 'S')
                        ++j;
                    else if (op == 'I')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first insertion or not
                        score += scoring.insertion();

                        // mover in the pattern
                        ++j;
                    }
                    else if (op == 'D')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first deletion or not
                        score += scoring.deletion();

                        // move in the text
                        ++k;
                    }
                    else if (op == 'M')
                    {
                        // add a match or mismatch score
                        score += (str[j] == ref[k]) ? scoring.match() : scoring.mismatch();

                        // move diagonally
                        ++j;
                        ++k;
                    }
                }
                return score;
            }

            // compute the score of the resulting alignment
            //
            template <AlignmentType TYPE, typename scoring_type>
            int32 score(
                GotohAligner<TYPE, scoring_type> aligner,
                const uint32 offset,
                const uint8 *str,
                const uint8 *ref)
            {
                const scoring_type &scoring = aligner.scheme;

                int32 score = 0;
                char state = 'S';

                for (uint32 a = 0, j = 0, k = offset; a < aln_len; ++a)
                {
                    const char op = aln[aln_len - a - 1];

                    if (op == 'S')
                        ++j;
                    else if (op == 'I')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first insertion or not
                        score += (state != op) ? scoring.pattern_gap_open() : scoring.pattern_gap_extension();

                        // mover in the pattern
                        ++j;
                    }
                    else if (op == 'D')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first deletion or not
                        score += (state != op) ? scoring.text_gap_open() : scoring.text_gap_extension();

                        // move in the text
                        ++k;
                    }
                    else if (op == 'M')
                    {
                        // add a match or mismatch score
                        score += (str[j] == ref[k]) ? scoring.match() : scoring.mismatch();

                        // move diagonally
                        ++j;
                        ++k;
                    }

                    // store new state
                    state = op;
                }
                return score;
            }

            // compute the score of the resulting alignment
            //
            template <AlignmentType TYPE, typename scoring_type>
            int32 score(
                WfahAligner<TYPE, scoring_type> aligner,
                const uint32 offset,
                const uint8 *str,
                const uint8 *ref)
            {
                const scoring_type &scoring = aligner.scheme;

                int32 score = 0;
                char state = 'S';

                for (uint32 a = 0, j = 0, k = offset; a < aln_len; ++a)
                {
                    const char op = aln[aln_len - a - 1];

                    if (op == 'S')
                        ++j;
                    else if (op == 'I')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first insertion or not
                        score += (state != op) ? scoring.pattern_gap_open() + scoring.pattern_gap_extension() : scoring.pattern_gap_extension();

                        // mover in the pattern
                        ++j;
                    }
                    else if (op == 'D')
                    {
                        // apply a gap-open or gap-extension penalty according to whether it's the first deletion or not
                        score += (state != op) ? scoring.text_gap_open() + scoring.text_gap_extension() : scoring.text_gap_extension();

                        // move in the text
                        ++k;
                    }
                    else if (op == 'M')
                    {
                        // add a match or mismatch score
                        score += (str[j] == ref[k]) ? scoring.match() : scoring.mismatch();

                        // move diagonally
                        ++j;
                        ++k;
                    }

                    // store new state
                    state = op;
                }
                return score;
            }

            char aln[1024];
            uint32 aln_len;
        };

    } // namespace sw
} // namespace nvbio

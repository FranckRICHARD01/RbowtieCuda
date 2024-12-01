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

#pragma once

#include <nvbio/basic/types.h>
#include <nvbio/basic/numbers.h>
#include <nvbio/alignment/wfa.h>

namespace nvbio
{
    namespace aln
    {
       namespace priv
        {
            namespace banded
            {

                ///@addtogroup private
                ///@{

#define SW_F_INITIALIZED_TO_INF

                // initialize the zero-th row of the banded DP-matrix with Wfah scoring
                //
                template <uint32 BAND_LEN, AlignmentType TYPE, typename H_band_type, typename F_band_type, typename scoring_type, typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init_row_zero_wfah(
                    H_band_type &H_band,
                    F_band_type &F_band,
                    const scoring_type &scoring,
                    const score_type infimum)
                {
                    H_band[0] = 0;
#pragma unroll
                    for (uint32 j = 1; j < BAND_LEN; ++j)
                        H_band[j] = TYPE == GLOBAL ? scoring.text_gap_open() + (j - 1) * scoring.text_gap_extension() : 0;

#if defined(SW_F_INITIALIZED_TO_INF)
// F[0,*] = -inf
#pragma unroll
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        F_band[j] = infimum;
#else
// F[0,*] = 0
#pragma unroll
                    for (uint32 j = 0; j < BAND_LEN; ++j)
                        F_band[j] = 0;
#endif

                    // The following is the initialization code for the "post" F-loop, i.e. for the case in which
                    // F is updated _after_ updating H
                    // F[0,*] = -inf, so F[1,*] is the gap open cost, except F[1,BAND_LEN-1] remains -inf
                    // #pragma unroll
                    // for (uint32 j = 0; j < BAND_LEN-1; ++j)
                    //    F_band[j] = scoring.pattern_gap_open();
                }

                ///
                /// A helper scoring context class, which can be used to adapt the basic
                /// wfah_alignment_score_dispatch algorithm to various situations, such as:
                ///   scoring
                ///   scoring within a window (i.e. saving only the last band within the window)
                ///   computing checkpoints
                ///   computing a flow submatrix
                ///
                template <uint32 BAND_LEN, AlignmentType TYPE>
                struct WfahScoringContext
                {
                    template <typename H_band_type, typename F_band_type, typename scoring_type, typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                        const uint32 i,
                        H_band_type &H_band,
                        F_band_type &F_band,
                        const scoring_type &scoring,
                        const score_type infimum)
                    {
                        init_row_zero_wfah<BAND_LEN, TYPE>(H_band, F_band, scoring, infimum);
                    }

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band) {}

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band) {}

                    // dir, edir, fdir are backtracking arrows for H,E and F submatrices, respectively at this cell
                    template <typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                        const uint32 i,
                        const uint32 j,
                        const score_type score,
                        const DirectionVector dir,
                        const DirectionVector edir,
                        const DirectionVector fdir) {                          
                        }

                    template <
                        typename context_type,
                        typename scoring_type,
                        typename score_type,
                        typename wfa_type,
                        typename text_cache_type,
                        typename ref_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void create_back( 
                        context_type &context,                      
                        score_type last_s,
                        const ref_type ref,
                        const text_cache_type q_cache,
                        const score_type alignement_offset,
                        const score_type alignement_k,
                        const score_type block,
                        const score_type M,
                        const score_type N,
                        const score_type V,
                        const score_type S,
                        const score_type G_o,
                        const score_type G_e,
                        const score_type zero,
                        const scoring_type &scoring,
                        wfa_type &wfa)
                    {}
                };

                ///
                /// A helper scoring context class, instantiated to configure wfah_alignment_score_dispatch
                /// to perform scoring in windows
                ///
                template <uint32 BAND_LEN, AlignmentType TYPE, typename checkpoint_type>
                struct WfahCheckpointedScoringContext
                {
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    WfahCheckpointedScoringContext(checkpoint_type checkpoints) : m_checkpoints(checkpoints) {}

                    template <typename H_band_type, typename F_band_type, typename scoring_type, typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                        const uint32 i,
                        H_band_type &H_band,
                        F_band_type &F_band,
                        const scoring_type &scoring,
                        const score_type infimum)
                    {
                        // check whether this is the first row
                        /*if (i == 0)
                            init_row_zero<BAND_LEN, TYPE>(H_band, F_band, scoring, infimum);
                        else
                        {
// load from checkpoint
#pragma unroll
                            for (uint32 j = 0; j < BAND_LEN; ++j)
                            {
                                H_band[j] = m_checkpoints[j];
                                F_band[j] = m_checkpoints[j];
                            }
                        }*/
                    }

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band) {}

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band)
                    {
                       /* const short infimum = Field_traits<short>::min() + 32;

// save the last row
#pragma unroll
                        for (uint32 j = 0; j < BAND_LEN; ++j)
                            m_checkpoints[j] = make_short2(
                                nvbio::max(H_band[j], infimum),
                                nvbio::max(F_band[j], infimum));*/
                    }

                    // dir, edir, fdir are backtracking arrows for H,E and F submatrices, respectively at this cell
                    template <typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                        const uint32 i,
                        const uint32 j,
                        const score_type score,
                        const DirectionVector dir,
                        const DirectionVector edir,
                        const DirectionVector fdir) {
                                                   }

                    template <
                        typename context_type,
                        typename scoring_type,
                        typename score_type,
                        typename wfa_type,
                        typename text_cache_type,
                        typename ref_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void create_back( 
                        context_type &context,                      
                        score_type last_s,
                        const ref_type ref,
                        const text_cache_type q_cache,
                        const score_type alignement_offset,
                        const score_type alignement_k,
                        const score_type block,
                        const score_type M,
                        const score_type N,
                        const score_type V,
                        const score_type S,
                        const score_type G_o,
                        const score_type G_e,
                        const score_type zero,
                        const scoring_type &scoring,
                        wfa_type &wfa)
                    {}

                    checkpoint_type m_checkpoints;
                };

                ///
                /// A helper scoring context class, instantiated to configure wfah_alignment_score_dispatch
                /// to store checkpoints of the banded DP matrix along the pattern
                ///
                template <uint32 BAND_LEN, AlignmentType TYPE, uint32 CHECKPOINTS, typename checkpoint_type>
                struct WfahCheckpointContext
                {
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    WfahCheckpointContext(checkpoint_type checkpoints) : m_checkpoints(checkpoints) {}

                    template <typename H_band_type, typename F_band_type, typename scoring_type, typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                        const uint32 i,
                        H_band_type &H_band,
                        F_band_type &F_band,
                        const scoring_type &scoring,
                        const score_type infimum)
                    {
                        init_row_zero<BAND_LEN, TYPE>(H_band, F_band, scoring, infimum);
                    }

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band)
                    {
                        // save checkpoint
                        /*if ((i & (CHECKPOINTS - 1)) == 0u)
                        {
                            const short infimum = Field_traits<short>::min() + 32;

                            const uint32 r = BAND_LEN * (i / CHECKPOINTS);
                            for (uint32 j = 0; j < BAND_LEN; ++j)
                                m_checkpoints[r + j] = make_short2(
                                    nvbio::max(H_band[j], infimum),
                                    nvbio::max(F_band[j], infimum));
                        }*/
                    }

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band)
                    {
                        /*const short infimum = Field_traits<short>::min() + 32;

                        // save the last row
                        const uint32 checkpoint_row = BAND_LEN * (i + CHECKPOINTS - 1) / CHECKPOINTS;
                        for (uint32 j = 0; j < BAND_LEN; ++j)
                            m_checkpoints[checkpoint_row + j] = make_short2(
                                nvbio::max(H_band[j], infimum),
                                nvbio::max(F_band[j], infimum));*/
                    }

                    // dir, edir, fdir are backtracking arrows for H,E and F submatrices, respectively at this cell
                    template <typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                        const uint32 i,
                        const uint32 j,
                        const score_type score,
                        const DirectionVector dir,
                        const DirectionVector edir,
                        const DirectionVector fdir) {                          
                        }

                    template <
                        typename context_type,
                        typename scoring_type,
                        typename score_type,
                        typename wfa_type,
                        typename text_cache_type,
                        typename ref_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void create_back( 
                        context_type &context,                      
                        score_type last_s,
                        const ref_type ref,
                        const text_cache_type q_cache,
                        const score_type alignement_offset,
                        const score_type alignement_k,
                        const score_type block,
                        const score_type M,
                        const score_type N,
                        const score_type V,
                        const score_type S,
                        const score_type G_o,
                        const score_type G_e,
                        const score_type zero,
                        const scoring_type &scoring,
                        wfa_type &wfa)
                    {}

                    checkpoint_type m_checkpoints;
                };

                ///
                /// A helper scoring context class, instantiated to configure wfah_alignment_score_dispatch
                /// to keep track of the direction vectors of a DP submatrix between given checkpoints
                ///
                template <uint32 BAND_LEN, AlignmentType TYPE, uint32 CHECKPOINTS, typename checkpoint_type, typename submatrix_type>
                struct WfahSubmatrixContext
                {
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    WfahSubmatrixContext(
                        const checkpoint_type checkpoints,
                        const uint32 checkpoint_id,
                        submatrix_type submatrix) : m_checkpoints(checkpoints),
                                                    m_checkpoint_id(checkpoint_id),
                                                    m_submatrix(submatrix) {}

                    template <typename H_band_type, typename F_band_type, typename scoring_type, typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                        const uint32 i,
                        H_band_type &H_band,
                        F_band_type &F_band,
                        const scoring_type &scoring,
                        const score_type infimum)
                    {
                        /*const uint32 r = m_checkpoint_id * BAND_LEN;

// restore the checkpoint
#pragma unroll
                        for (uint32 j = 0; j < BAND_LEN; ++j)
                        {
                            const short2 f = m_checkpoints[r + j];
                            H_band[j] = f.x;
                            F_band[j] = f.y;
                        }*/
                    }

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band) {}

                    template <typename H_band_type, typename F_band_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_row(
                        const uint32 i,
                        const H_band_type &H_band,
                        const F_band_type &F_band) {}

                    template <typename score_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                        const uint32 i,
                        const uint32 j,
                        const score_type score,
                        const DirectionVector hdir,
                        const DirectionVector edir,
                        const DirectionVector fdir)
                    {
                        // save the direction vectors for H,E,F
                        const uint8 cdir =
                            (TYPE == LOCAL ? (score == 0 ? SINK : hdir) : hdir) | edir | fdir;

                        const uint32 offset = m_checkpoint_id * CHECKPOINTS;

                        /*unsigned  long int c = (i - offset) * BAND_LEN + j;
                       

                        if (c > 100000)
                        {
                             const int  a = BAND_LEN;
                        const int  b = CHECKPOINTS;
                            c = c;
                        }*/

                       

                        m_submatrix[(i - offset) * BAND_LEN + j] = cdir;
                    }

                    template <
                        typename context_type,
                        typename scoring_type,
                        typename score_type,
                        typename wfa_type,
                        typename text_cache_type,
                        typename ref_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void create_back( 
                        context_type &context,                      
                        score_type last_s,
                        const ref_type ref,
                        const text_cache_type q_cache,
                        const score_type alignement_offset,
                        const score_type alignement_k,
                        const score_type block,
                        const score_type M,
                        const score_type N,
                        const score_type V,
                        const score_type S,
                        const score_type G_o,
                        const score_type G_e,
                        const score_type zero,
                        const scoring_type &scoring,
                        wfa_type &wfa)
                    {
                        score_type i, j, h, v;

                        score_type score = last_s;

                        int32 x = M;
                        int32 y = N;
                        
                        DirectionVector hdir = SUBSTITUTION;
                        DirectionVector hdir2 = SUBSTITUTION;

                        score_type offset = alignement_offset;

                        score_type k = alignement_k;

                        h = offset - k;
                        v = offset;

                        if (h < N)
                        {
                            score_type i = N - h;

                            while (i > 0)
                            {
                                // verif_back += 'D';

                                wfa.add_backtrace_etape(
                                    context,
                                    x,
                                    y,
                                    h,
                                    v,
                                    block,
                                    ref,
                                    q_cache,
                                    DELETION,
                                    SUBSTITUTION,
                                    DELETION_EXT,
                                    scoring);

                                i--;

                                y--;
                            }
                        }

                        if (v < M)
                        {
                            score_type i = N - h;

                            while (i > 0)
                            {
                                // verif_back += 'I';

                                wfa.add_backtrace_etape(
                                    context,
                                    x,
                                    y,
                                    h,
                                    v,
                                    block,
                                    ref,
                                    q_cache,
                                    INSERTION,
                                    SUBSTITUTION,
                                    INSERTION_EXT,
                                    scoring);

                                i--;

                                x--;
                                y--;
                            }
                        }

                        v = M;
                        h = N;

                        int32 pos = 1;

                        while (score > 0 && v > 0 && h > 0 && pos > 0)
                        {
                            score_type mismatch = wfa.f(score - S, wfa.Point_H_BAND);
                            score_type gap_open1 = wfa.f(score - G_o - G_e, wfa.Point_H_BAND);
                            score_type gap_extend1 = wfa.f(score - G_e, wfa.Point_H_BAND);

                            if (hdir == SUBSTITUTION)
                            {
                                hdir2 = SUBSTITUTION;

                                score_type pos_del, pos_ins;
                                DirectionVector hdir_del, hdir_ins;

                                // NVBIO_CUDA_DEBUG_ASSERT(deltaHk >= 0, "Problem indice K >= 0, DELTAV=%d k=%d score=%d\n", DELTAV, k, score);

                                score_type pos_mismatch = wfa.H_Band.get(mismatch, k) + 1;

                                score_type pos_del_open = wfa.H_Band.get(gap_open1, k + 1);
                                score_type pos_del_ext = wfa.F_Band.get(gap_extend1, k + 1);

                                if (pos_del_ext >= pos_del_open)
                                {
                                    pos_del = pos_del_ext;
                                    hdir_del = DELETION_EXT;
                                }
                                else
                                {
                                    pos_del = pos_del_open;
                                    hdir_del = DELETION;
                                }

                                score_type pos_ins_open = wfa.H_Band.get(gap_open1, k - 1);
                                score_type pos_ins_ext = wfa.E_Band.get(gap_extend1, k - 1);

                                if (pos_ins_ext >= pos_ins_open)
                                {
                                    pos_ins = pos_ins_ext + 1;
                                    hdir_ins = INSERTION_EXT;
                                }
                                else
                                {
                                    pos_ins = pos_ins_open + 1;
                                    hdir_ins = INSERTION;
                                }

                                if (pos_del >= pos_ins)
                                {
                                    pos = pos_del;
                                    hdir2 = hdir_del;
                                }
                                else
                                {
                                    pos = pos_ins;
                                    hdir2 = hdir_ins;
                                }

                                if (pos_mismatch >= pos)
                                {
                                    pos = pos_mismatch;
                                    hdir2 = SUBSTITUTION;
                                }

                                /*if (pos <0)
                                {
                                    pos=pos;
                                }*/

                                NVBIO_CUDA_DEBUG_ASSERT(pos >= 0 && pos < 1000, "(BAND) Problem pos result SUBSTITUTION, pos=%d  score=%d last_score=%d h=%d v=%d\n", pos, score, last_s, h, v );
                            }
                            else if (hdir == INSERTION)
                            {
                                score_type E = wfa.E_Band.get(gap_extend1, k - 1);
                                score_type H = wfa.H_Band.get(gap_open1, k - 1);

                                if (E >= H)
                                {
                                    pos = E + 1;
                                    hdir2 = INSERTION_EXT;
                                }
                                else
                                {
                                    pos = H + 1;
                                    hdir2 = INSERTION;
                                }

                                // char ref_str[1024];
                                // char read_str[1024];

                                /*dna_text_cache_to_string( ref, N, ref_str );

                                dna_ref_type_to_string( q_cache, M, read_str );*/

                                // NVBIO_CUDA_DEBUG_ASSERT(pos >= 0, "Problem pos result INSERTION, pos=%d\nref=%s\ntext=%s\n", pos, ref_str, read_str );
                            }
                            else if (hdir == DELETION)
                            {
                                score_type F = wfa.F_Band.get(gap_extend1, k + 1);
                                score_type H = wfa.H_Band.get(gap_open1, k + 1);

                                if (F >= H)
                                {
                                    pos = F;
                                    hdir2 = DELETION_EXT;
                                }
                                else
                                {
                                    pos = H;
                                    hdir2 = DELETION;
                                }

                                // char ref_str[1024];
                                // char read_str[1024];

                                /*dna_text_cache_to_string( ref, N, ref_str );

                                dna_ref_type_to_string( q_cache, M, read_str );*/

                                // NVBIO_CUDA_DEBUG_ASSERT(pos >= 0, "Problem pos result DELETION, pos=%d\nref=%s\ntext=%s\n", pos, ref_str, read_str);
                            }

                            if (pos > 0)
                            {
                                if (hdir == SUBSTITUTION)
                                {
                                    int32 offset_max = pos;

                                    if (offset - offset_max > 0)
                                    {
                                        for (int i = 1; i <= offset - offset_max; i++)
                                        {
                                            // verif_back += 'M';

                                            wfa.add_backtrace_etape(
                                                context,
                                                x,
                                                y,
                                                h,
                                                v,
                                                block,
                                                ref,
                                                q_cache,
                                                SUBSTITUTION,
                                                SUBSTITUTION,
                                                SUBSTITUTION,
                                                scoring);

                                            x--;

                                            h--;
                                            v--;
                                        }

                                        offset = offset_max;
                                    }
                                    
                                    h = offset - k;
                                    v = offset;

                                    if (h <= 0 || v <= 0)
                                        break;
                                }

                                if (hdir2 == SUBSTITUTION)
                                {
                                    // verif_back += 'X';

                                    wfa.add_backtrace_etape(
                                        context,
                                        x,
                                        y,
                                        h,
                                        v,
                                        block,
                                        ref,
                                        q_cache,
                                        SUBSTITUTION,
                                        SUBSTITUTION,
                                        SUBSTITUTION,
                                        scoring);

                                    x--;

                                    score -= S;
                                    hdir = SUBSTITUTION;
                                    offset--;
                                }
                                else if (hdir2 == INSERTION)
                                {
                                    // verif_back += 'I';

                                    wfa.add_backtrace_etape(
                                        context,
                                        x,
                                        y,
                                        h,
                                        v,
                                        block,
                                        ref,
                                        q_cache,
                                        INSERTION,
                                        SUBSTITUTION,
                                        INSERTION_EXT,
                                        scoring);

                                    x--;
                                    y--;

                                    score -= G_o + G_e;
                                    hdir = SUBSTITUTION;
                                    k--;
                                    offset--;
                                }
                                else if (hdir2 == INSERTION_EXT)
                                {
                                    //  verif_back += 'I';

                                    wfa.add_backtrace_etape(
                                        context,
                                        x,
                                        y,
                                        h,
                                        v,
                                        block,
                                        ref,
                                        q_cache,
                                        INSERTION,
                                        SUBSTITUTION,
                                        INSERTION_EXT,
                                        scoring);

                                    x--;
                                    y--;

                                    score -= G_e;
                                    hdir = INSERTION;
                                    k--;
                                    offset--;
                                }
                                else if (hdir2 == DELETION)
                                {
                                    //  verif_back += 'D';

                                    score -= G_o + G_e;
                                    hdir = SUBSTITUTION;
                                    k++;

                                    wfa.add_backtrace_etape(
                                        context,
                                        x,
                                        y,
                                        h,
                                        v,
                                        block,
                                        ref,
                                        q_cache,
                                        DELETION,
                                        SUBSTITUTION,
                                        DELETION_EXT,
                                        scoring);

                                    y--;
                                }
                                else if (hdir2 == DELETION_EXT)
                                {
                                    //  verif_back += 'D';

                                    score -= G_e;
                                    hdir = DELETION;
                                    k++;

                                    wfa.add_backtrace_etape(
                                        context,
                                        x,
                                        y,
                                        h,
                                        v,
                                        block,
                                        ref,
                                        q_cache,
                                        DELETION,
                                        SUBSTITUTION,
                                        DELETION_EXT,
                                        scoring);

                                    y--;
                                }

                                h = offset - k;
                                v = offset;
                            }
                        }

                        score = last_s;

                        /*if (TYPE == LOCAL)
                        {
                            if (CHECK_M)
                            {
                                // during local alignment we save the best score across all bands
                                for (int32 j = 1; j <= (int32)BAND_LEN; j++)
                                {
                                    if (block + j <= M)
                                        sink.report(H_band[j], make_uint2(line, block + j));
                                }
                            }
                            else
                            {
                                // during local alignment we save the best score across all bands
                                for (int32 j = 1; j <= (int32)BAND_LEN; j++)
                                    sink.report(H_band[j], make_uint2(line, block + j));
                            }
                        }
                        else if (CHECK_M)
                        {
                            if (TYPE == SEMI_GLOBAL)
                            {
                                // during semi-global alignment we save the best score across the last column H[*][M], at each row
                                save_boundary<BAND_LEN>(block, M, H_band, line, sink);
                            }
                        }*/

                        /////////////////
 
                        if (v > 0 && h > 0)
                        {
                            int32 num_matches = min(h, v);

                            for (int32 i = 1; i <= num_matches; i++)
                            {
                                // verif_back += 'M';

                                wfa.add_backtrace_etape(
                                    context,
                                    x,
                                    y,
                                    h,
                                    v,
                                    block,
                                    ref,
                                    q_cache,
                                    SUBSTITUTION,
                                    SUBSTITUTION,
                                    SUBSTITUTION,
                                    scoring);

                                    h--;
                                    v--;

                                    x--;
                            }
                        }

                        while (h > 0)
                        {
                            // verif_back += 'D';                           

                            wfa.add_backtrace_etape(
                                context,
                                x,
                                y,
                                h,
                                v,
                                block,
                                ref,
                                q_cache,
                                DELETION,
                                SUBSTITUTION,
                                DELETION_EXT,
                                scoring);

                            h--;

                            y--;
                        }

                        while (v > 0)
                        {
                            // verif_back += 'I';

                            wfa.add_backtrace_etape(
                                context,
                                x,
                                y,
                                h,
                                v,
                                block,
                                ref,
                                q_cache,
                                INSERTION,
                                INSERTION_EXT,
                                SUBSTITUTION,
                                scoring);

                            v--;

                            x--;
                            y--; 
                        }

                        while (y > 0)
                        {
                            // verif_back += 'D';                           

                            wfa.add_backtrace_etape(
                                context,
                                x,
                                y,
                                h,
                                v,
                                block,
                                ref,
                                q_cache,
                                SUBSTITUTION,
                                SUBSTITUTION,
                                SUBSTITUTION,
                                scoring);                          

                            y--;
                        }

                        while (x > 0)
                        {
                            // verif_back += 'D';                           

                            wfa.add_backtrace_etape(
                                context,
                                x,
                                y,
                                h,
                                v,
                                block,
                                ref,
                                q_cache,
                                SUBSTITUTION,
                                SUBSTITUTION,
                                SUBSTITUTION,
                                scoring);                          

                            x--;
                        }
                    }


                    checkpoint_type m_checkpoints;
                    uint32 m_checkpoint_id;
                    submatrix_type m_submatrix;
                };

                ///
                /// The template struct to dispatch calls to banded_alignment_score.
                ///
                template <uint32 BAND_LEN, AlignmentType TYPE>
                struct wfah_alignment_score_dispatch
                {                    
                    // update a DP row
                    //
                    template <
                        bool CHECK_M,
                        typename context_type,
                        typename text_cache_type,
                        typename score_type,
                        typename ref_type,
                        typename sink_type,
                        typename scoring_type,
                        typename wfa_type>
                    NVBIO_FORCEINLINE
#ifdef WFA_TESTS
                        NVBIO_HOST_DEVICE
                    #else
                        // NVBIO_DEVICE
                        NVBIO_HOST_DEVICE
                    #endif
                        static score_type
                        update(
                            context_type &context,
                            const score_type block,
                            const score_type M,
                            score_type i,
                            const score_type N,
                            const ref_type ref,
                            const text_cache_type q_cache,
                            score_type &temp_i,
                            score_type *H_band,
                            score_type *F_band,
                            sink_type &sink,
                            const score_type min_score,
                            const score_type V,
                            const score_type S,
                            const score_type G_o,
                            const score_type G_e,
                            const score_type zero,
                            const scoring_type &scoring,
                            wfa_type &wfa)
                    {
#ifndef WFA_TESTS
// #define SHARED __shared__
#define SHARED
#else
#define SHARED
#endif
                        
                    score_type s = 0;
                    score_type ss = 0;
                    score_type last_s = 0;
                    score_type indice = 0;

                    score_type mismatch = -1;
                    score_type gap_open1 = -1;
                    score_type gap_extend1 = -1;

                    score_type MN = M - N;
                    score_type alignment_k = MN;
                    score_type alignment_offset = M;
                    score_type absMN = (MN >= 0) ? MN : -MN;

                    score_type h = 0;
                    score_type v = 0;

                    score_type score = 0;
                    score_type lo = 0;
                    score_type hi = 0;
                    score_type k = 0;
                    score_type delta = 0;

                    int16 *H_Band_gap_open1_ptr = nullptr;
                    int16 *H_Band_gap_mismatch_ptr = nullptr;
                    int16 *E_Band_gap_extend1_ptr = nullptr;
                    int16 *F_Band_gap_extend1_ptr = nullptr;
                    int16 *H_Band_ptr = nullptr;
                    int16 *E_Band_ptr = nullptr;
                    int16 *F_Band_ptr = nullptr;

                    int32 pos_gap_extend1, pos_gap_open1, pos_mismatch;
                    int32 pos_ss = 0;

                    // SHARED int16 H[2 * DIM_SHARED + 1];
                    // SHARED int16 I[2 * DIM_SHARED + 1];
                    // SHARED int16 D[2 * DIM_SHARED + 1];
                    
                    wfa.H_Band.set(0, 0, 0);
                    wfa.E_Band.set(0, 0, WFA_MIN1);
                    wfa.F_Band.set(0, 0, WFA_MIN1);
                    wfa.H_Band.set_null(s, false);
                    wfa.E_Band.set_null(s, true);
                    wfa.F_Band.set_null(s, true);
                    wfa.H_Band.set_hi_lo(0, true);
                    wfa.F_Band.set_hi_lo(0, false);
                    wfa.E_Band.set_hi_lo(0, false);
                    // H[delta] = 0;
                    // D[delta] = WFA_MIN1;
                    // I[delta] = WFA_MIN1;
                    H_Band_ptr = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(0, DIM_SHARED)];
                    /*#ifndef WFA_TESTS
                                                    __syncthreads();
                    #endif*/
                    SHARED score_type I1, I2, IM, D1, D2, DM, HM, offset;

                    do
                    {
                        if (H_Band_ptr)
                        {
#pragma unroll
                            for (k = lo; k <= hi; k++)
                            {
                                offset = *(H_Band_ptr + k);
                                // H[k + delta];

                                if (offset < 0)
                                    continue;

                                uint16 h = uint16(offset - k);

                                const score_type MAX_ITER = M;
#pragma unroll
                                for (int i = 0; i < MAX_ITER && offset < M && h < N && q_cache[offset].x == ref[h]; ++i)
                                {
                                    offset++;
                                    h++;
                                }

                                *(H_Band_ptr + k) = offset;

                                /*#if __CUDA_ARCH__ >= 800
                                    __stwt(&H_Band_ptr[k], offset);
                                #else
                                    H_Band_ptr[k] = offset;
                                #endif  */
                            }

                            /*#ifndef WFA_TESTS
                                                            __syncthreads();
                            #endif*/

                            wfa.set_HEF_hi_lo(ss, hi, lo);

                            wfa.trim_ends(&wfa.H_Band, ss, N, M);
                            wfa.trim_ends(&wfa.E_Band, ss, N, M);
                            wfa.trim_ends(&wfa.F_Band, ss, N, M);
                        }

                        if (wavefront_heuristic_cutoff(wfa, ss, s, N, M))
                            break;

                        if (ss >= nvbio::min(WFA_BAND_LEN2_LIMIT_BAND, WFA_MAX_SCORE))
                        {
                            /*char ref_str[1000];
                            char read_str[1000];

                            dna_ref_type_to_string(ref, N, ref_str);
                            dna_text_cache_to_string(q_cache, M, read_str);

                            NVBIO_CUDA_DEBUG_ASSERT(s < 0, "depassement wfah!\n1:%s\n2:%s\nlen1=%d\nlen2=%d\n", ref_str, read_str, N, M);
                            */
                            break;
                        }

                        if (wfa.testEnd(N, M, wfa.H_Band, ss, alignment_k, alignment_offset, true) ||
                            wfa.testEnd(N, M, wfa.F_Band, ss, alignment_k, alignment_offset, false) ||
                            wfa.testEnd(N, M, wfa.E_Band, ss, alignment_k, alignment_offset, false))
                        {
                            wfa.heuristic.alignement_end_pos_score = ss;
                            wfa.heuristic.alignement_end_pos_k = alignment_k;
                            wfa.heuristic.alignement_end_pos_offset = alignment_offset;

                            break;
                        }

                        if (wfa.H_Band.get_lo(ss) <= MN &&
                            MN <= wfa.H_Band.get_hi(ss) &&
                            wfa.H_Band.get(ss, MN) >= M)
                        {
                            break;
                        }

                        s++;

                        // wfa.Point_H_BAND[++indice] = s;

                        wfa.H_Band.set_null(s, false);
                        wfa.E_Band.set_null(s, false);
                        wfa.F_Band.set_null(s, false);
                        wfa.H_Band.set_hi_lo(s, true);
                        wfa.F_Band.set_hi_lo(s, false);
                        wfa.E_Band.set_hi_lo(s, false);
                        wfa.H_Band.set(s, 0, WFA_MIN1);
                        wfa.E_Band.set(s, 0, WFA_MIN1);
                        wfa.F_Band.set(s, 0, WFA_MIN1);                      
                        // H[delta] = WFA_MIN1;
                        // D[delta] = WFA_MIN1;
                        // I[delta] = WFA_MIN1;                        

                        ss = wfa.f(s, wfa.Point_H_BAND);                      

                        if (ss >= 0)
                        {
                            mismatch = wfa.f(s - S, wfa.Point_H_BAND);
                            gap_open1 = wfa.f(s - G_o - G_e, wfa.Point_H_BAND);
                            gap_extend1 = wfa.f(s - G_e, wfa.Point_H_BAND);

                            H_Band_gap_open1_ptr = nullptr;
                            H_Band_gap_mismatch_ptr = nullptr;
                            E_Band_gap_extend1_ptr = nullptr;
                            F_Band_gap_extend1_ptr = nullptr;
                            H_Band_ptr = nullptr;
                            E_Band_ptr = nullptr;
                            F_Band_ptr = nullptr;

                            if (mismatch >= 0)
                            {
                                pos_mismatch = COMPUTE_DIM_SHARED(mismatch, DIM_SHARED);
                                H_Band_gap_mismatch_ptr = wfa.H_Band.get_null(mismatch) ? nullptr : &wfa.H_Band.scores[pos_mismatch - wfa.H_Band.get_lo(mismatch)];
                            }
                            if (gap_open1 >= 0)
                            {
                                pos_gap_open1 = COMPUTE_DIM_SHARED(gap_open1, DIM_SHARED);
                                H_Band_gap_open1_ptr = wfa.H_Band.get_null(gap_open1) ? nullptr : &wfa.H_Band.scores[pos_gap_open1 - wfa.H_Band.get_lo(gap_open1)];
                            }
                            if (gap_extend1 >= 0)
                            {
                                pos_gap_extend1 = COMPUTE_DIM_SHARED(gap_extend1, DIM_SHARED);
                                E_Band_gap_extend1_ptr = wfa.E_Band.get_null(gap_extend1) ? nullptr : &wfa.E_Band.scores[pos_gap_extend1 - wfa.E_Band.get_lo(gap_extend1)];
                                F_Band_gap_extend1_ptr = wfa.F_Band.get_null(gap_extend1) ? nullptr : &wfa.F_Band.scores[pos_gap_extend1 - wfa.F_Band.get_lo(gap_extend1)];
                            }

                            if (H_Band_gap_mismatch_ptr || H_Band_gap_open1_ptr || E_Band_gap_extend1_ptr || F_Band_gap_extend1_ptr)
                            {
                                /*#ifndef WFA_TESTS
                                                                    __syncthreads();
                                #endif*/                                

                                score_type H_lo_m = (H_Band_gap_mismatch_ptr) ? wfa.H_Band.get_lo(mismatch) : 1;
                                score_type H_hi_m = (H_Band_gap_mismatch_ptr) ? wfa.H_Band.get_hi(mismatch) : -1;

                                score_type H_lo_o = (H_Band_gap_open1_ptr) ? wfa.H_Band.get_lo(gap_open1) : 1;
                                score_type H_hi_o = (H_Band_gap_open1_ptr) ? wfa.H_Band.get_hi(gap_open1) : -1;

                                score_type E_lo = (E_Band_gap_extend1_ptr) ? wfa.E_Band.get_lo(gap_extend1) : 1;
                                score_type E_hi = (E_Band_gap_extend1_ptr) ? wfa.E_Band.get_hi(gap_extend1) : -1;

                                score_type F_lo = (F_Band_gap_extend1_ptr) ? wfa.F_Band.get_lo(gap_extend1) : 1;
                                score_type F_hi = (F_Band_gap_extend1_ptr) ? wfa.F_Band.get_hi(gap_extend1) : -1;

                                lo = nvbio::min(nvbio::min(H_lo_o - 1, H_lo_m), nvbio::min(E_lo + 1, F_lo - 1));
                                hi = nvbio::max(nvbio::max(H_hi_o + 1, H_hi_m), nvbio::max(E_hi + 1, F_hi - 1));

                                pos_ss = COMPUTE_DIM_SHARED(ss, DIM_SHARED);
                                H_Band_ptr = &wfa.H_Band.scores[pos_ss - lo];
                                E_Band_ptr = &wfa.E_Band.scores[pos_ss - lo];
                                F_Band_ptr = &wfa.F_Band.scores[pos_ss - lo];

                                // delta = -lo;
#pragma unroll 1
                                for (score_type k = lo; k <= hi; ++k)
                                {
                                    I1 = (E_Band_gap_extend1_ptr && k - 1 >= E_lo && k - 1 <= E_hi) ? E_Band_gap_extend1_ptr[k - 1] : WFA_MIN1;
                                    I2 = (H_Band_gap_open1_ptr && k - 1 >= H_lo_o && k - 1 <= H_hi_o) ? H_Band_gap_open1_ptr[k - 1] : WFA_MIN1;
                                    IM = nvbio::max(I1, I2) + 1u;
// I[k + delta] = IM;
#if __CUDA_ARCH__ >= 800
                                    __stwt(&E_Band_ptr[k], IM);
#else
                                    E_Band_ptr[k] = IM;
#endif

                                    D1 = (F_Band_gap_extend1_ptr && k + 1 >= F_lo && k + 1 <= F_hi) ? F_Band_gap_extend1_ptr[k + 1] : WFA_MIN1;
                                    D2 = (H_Band_gap_open1_ptr && k + 1 >= H_lo_o && k + 1 <= H_hi_o) ? H_Band_gap_open1_ptr[k + 1] : WFA_MIN1;
                                    DM = nvbio::max(D1, D2);
                                    // D[k + delta] = DM;

#if __CUDA_ARCH__ >= 800
                                    __stwt(&F_Band_ptr[k], DM);
#else
                                    F_Band_ptr[k] = DM;
#endif

                                    HM = (H_Band_gap_mismatch_ptr && k >= H_lo_m && k <= H_hi_m) ? H_Band_gap_mismatch_ptr[k] : WFA_MIN1;
                                    HM = nvbio::max3(HM + 1, IM, DM);

                                    if ((uint16)HM > M || (uint16)(HM - k) > N)
                                        HM = WAVEFRONT_OFFSET_NULL;

                                    H_Band_ptr[k] = HM;
                                    // H[k + delta] = HM;
                                }
                                /*#ifndef WFA_TESTS
                                                                    __syncthreads();
                                #endif*/

                                last_s = s;
                            }
                            // else
                            //     indice--;
                        }
                    } while (true);

                    // starting backtracing

                    context.create_back(context, last_s, ref, q_cache, alignment_offset, alignment_k, block, M, N, V, S, G_o, G_e, zero, scoring, wfa);

#ifndef WFA_TESTS
                    last_s -= N - M;

                    if (last_s < 0)
                        last_s = 0;
#endif

                    return last_s;
                    }

                    ///
                    /// Calculate the banded alignment score between a string
                    /// and a reference, using the Smith-Waterman algorithm with Wfah's scoring
                    /// scheme.
                    ///
                    /// \tparam context_type    a context class used to configure the behaviour
                    ///                         of the algorithm, initializing and capturing the output DP matrix;
                    ///                         it must implement the following interface:
                    ///
                    /// \code
                    /// struct context_type
                    /// {
                    ///     // initialize the first DP band
                    ///     template <typename scoring_type>
                    ///     void init(
                    ///         const uint32        i,
                    ///               int32*        H_band,
                    ///               int32*        F_band,
                    ///         const scoring_type& scoring,
                    ///         const int32         infimum);
                    ///
                    ///     // do anything with the i-th band
                    ///     void previous_row(
                    ///         const uint32        i,
                    ///         const int32*        H_band,
                    ///         const int32*        F_band);
                    ///
                    ///     // do anything with the last computed band, at row i
                    ///     void last_row(
                    ///         const uint32        i,
                    ///         const int32*        H_band,
                    ///         const int32*        F_band);
                    ///
                    ///     // do anything with the cell (i,j)
                    ///     void new_cell(
                    ///         const uint32            i,
                    ///         const uint32            j,
                    ///         const int32             score,
                    ///         const DirectionVector   dir,
                    ///         const DirectionVector   edir,
                    ///         const DirectionVector   fdir);
                    /// };
                    /// \endcode
                    ///
                    /// \param scoring      scoring scheme
                    /// \param pattern      shorter string (horizontal)
                    /// \param text         longer string (vertical)
                    /// \param pos          offset in the reference string
                    /// \param context      the context class
                    /// \param sink         output alignment sink
                    ///
                    /// \return             false if the alignment didn't reach the
                    ///                     minimum score, true otherwise
                    ///
                    template <
                        typename pattern_type,
                        typename qual_type,
                        typename text_type,
                        typename scoring_type,
                        typename context_type,
                        typename sink_type,
                        typename wfa_type>
                    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool run(
                        const scoring_type &scoring,
                        pattern_type pattern,
                        qual_type quals,
                        text_type text,
                        const int32 window_begin,
                        const int32 window_end,
                        const int32 pos,
                        const int32 min_score,
                        context_type &context,
                        sink_type &sink,
                        wfa_type& wfa)
                    {
#ifdef WFA_TESTS
                    #define BAND_LEN_BAND 1000
#else
                    #define BAND_LEN_BAND 200
#endif
                                              
                        typedef int32 score_type;
   
                        score_type pattern_len = pattern.length();
                        score_type text_len = text.length();
                        const score_type start = pos;
                                                       
                        if (text_len < pattern_len)
                            return false;                       
                        
                        typedef uchar2 text_cache_type;
                        text_cache_type q_cache[BAND_LEN_BAND];
                        typedef typename Reference_cache<BAND_LEN_BAND>::type ref_type;
                        uint32 text_cache_storage[Reference_cache<BAND_LEN_BAND>::BAND_WORDS];
                        text_cache_storage[0]=0;
                        ref_type ref_cache(text_cache_storage);

                        score_type begin = 0;
                        score_type text_len_mem = text_len;

#ifndef WFA_TESTS
                        score_type end = text_len;
                        const score_type lim = 20;
                        const score_type threshold = 5;

                        if (pattern_len >= 150)
                        {
#pragma unroll
                            for (score_type i = 0; i < text_len - pattern_len + lim; i++)
                            {
                                score_type num = 0;

#pragma unroll
                                for (score_type j = 0; j < lim; j++)
                                {
                                    if (text[i + j] != pattern[j])
                                        num++;

                                    if (num >= threshold)
                                        break;
                                }
                                if (num < threshold)
                                {
                                    begin = i;
                                    break;
                                }
                            }

#pragma unroll
                            for (score_type i = text_len - 1; i > pattern_len - lim; i--)
                            {
                                score_type num = 0;
#pragma unroll
                                for (score_type j = pattern_len - 1; j > pattern_len - 1 - lim; j--)
                                {
                                    if (text[i - (pattern_len - 1 - j)] != pattern[j])
                                        num++;

                                    if (num >= threshold)
                                        break;
                                }
                                if (num < threshold)
                                {
                                    end = i;
                                    break;
                                }
                            }

                            if (begin >= end)
                            {
                                begin = 0;
                            }
                            else
                            {
                                text_len = nvbio::max(pattern_len, end - begin);
                                begin = nvbio::max(0, begin);
                            }
                        }
#endif

                        // load first band of text
#pragma unroll 1
                        for (int32 j = 0; j < text_len; j++)
                        {
                            if (start + window_begin + begin + j < text_len_mem)
                                ref_cache[j] = text[start + window_begin + begin + j];
                        }


                        const score_type G_o = abs(scoring.pattern_gap_open());
                        const score_type G_e = abs(scoring.pattern_gap_extension());
                        const score_type S = abs(scoring.mismatch());
                        const score_type V = abs(scoring.match());
                        const score_type zero = score_type(0);
                        const score_type infimum = WFA_MIN;/*Field_traits<short>::min() -
                                                   nvbio::max(nvbio::max(G_o, G_e),
                                                             nvbio::max(scoring.text_gap_open(), scoring.text_gap_extension()));*/

                        score_type max_score = Field_traits<score_type>::min();

                        // initialize the first band (corresponding to the 0-th row of the DP matrix)
                        score_type H_band[BAND_LEN];
                        score_type F_band[BAND_LEN];

                        // wfa
#ifndef WFA_TESTS
                        wavefront_heuristic_set_wfadaptive(wfa, DIM_SHARED2 * 2 , DIM_SHARED2 * 2, 1);
                        //wavefront_heuristic_set_wfmash(wfa, DIM_SHARED2, DIM_SHARED2 , 1);
                        //wavefront_heuristic_set_xdrop(wfa, 15, 1);
                        //wavefront_heuristic_set_zdrop(wfa, 15, 3);                        
                        //wavefront_heuristic_set_banded_static(wfa, -DIM_SHARED2, DIM_SHARED2);
                        wavefront_heuristic_set_banded_adaptive(wfa, -DIM_SHARED2, DIM_SHARED2, 1);
#endif

                        // initialize bands
                        context.init(
                            window_begin,
                            H_band,
                            F_band,
                            scoring,
                            infimum);

                                             
                        // initialize extra bands
                        #pragma unroll 1
                        for (int32 i = 0; i < WFA_BAND_LEN2_Y; i++)
                        {                           
                            wfa.Point_H_BAND[i] = 0;              
                        }

                        #pragma unroll 1
                        for (int32 t = 0; t < window_end - window_begin; t++)
                        {
                            q_cache[t] = make_uchar2(pattern[t + window_begin], quals[t + window_begin]);
                        }

                        pattern_len = window_end - window_begin;

                        score_type temp_i = H_band[0];

                        max_score = update<true>(
                            context,
                            (score_type)0,
                            pattern_len,
                            (score_type)0,
                            text_len,
                            ref_cache,
                            q_cache,
                            temp_i,
                            H_band,
                            F_band,
                            sink,
                            (score_type)min_score,
                            V,
                            S,
                            G_o,
                            G_e,
                            zero,
                            scoring,
                            wfa);

#ifndef WFA_TESTS
/*#pragma unroll
                        for (score_type i = text_len_mem - begin -1- 1; i > pattern_len - lim; i--)
                        {
                            score_type num = 0;
#pragma unroll
                            for (score_type j = pattern_len - 1; j > pattern_len - 1 - lim; j--)
                            {
                                if (text[i - (pattern_len - 1 - j)] != pattern[j])
                                    num++;

                                if (num >= 1)
                                    break;
                            }
                            if (num < 1)
                            {
                                text_len_mem = i + 1;
                                break;
                            }
                        }*/
#endif

                        // save the last column
                        // context.last_column(window_end, text_len, pattern_len, temp);

                        // if (TYPE == GLOBAL)
                        // save_Mth<BAND_LEN>( M, H_band, N-1, sink );

                        if (TYPE == GLOBAL)
                            sink.report( -max_score, make_uint2( text_len_mem, pattern_len ) );

                        if (TYPE == SEMI_GLOBAL)
                            sink.report( -max_score, make_uint2( text_len_mem, pattern_len ) );

                        return true;
                    }
                };

                /// @} // end of private group

            } // namespace banded

            ///@addtogroup private
            ///@{

            ///
            /// Calculate the banded alignment score between a pattern and a text string
            /// using the Smith-Waterman algorithm.
            ///
            /// \param pattern      shorter string (horizontal)
            /// \param quals        qualities string
            /// \param text         longer string (vertical)
            /// \param pos          offset in the reference string
            /// \param sink         output alignment sink
            ///
            /// \return             false if the minimum score was not reached, true otherwise
            ///
            template <
                uint32 BAND_LEN,
                AlignmentType TYPE,
                typename scoring_type,
                typename pattern_type,
                typename qual_type,
                typename text_type,
                typename sink_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE bool banded_alignment_score(
                const WfahAligner<TYPE, scoring_type> &aligner,
                pattern_type pattern,
                qual_type quals,
                text_type text,
                const int32 min_score,
                sink_type &sink,
                wfa_type &wfa)
            {
                if (1)
                {
                    priv::banded::WfahScoringContext<BAND_LEN, TYPE> context;

                    return priv::banded::wfah_alignment_score_dispatch<BAND_LEN, TYPE>::run(aligner.scheme, pattern, quals, text, 0u, pattern.length(), 0u, min_score, context, sink, wfa);
                }
                else
                {
                    /*priv::banded::GotohScoringContext<BAND_LEN, TYPE> context;

                    aln::SimpleGotohScheme scoring;
                    scoring.m_match = 0;
                    scoring.m_mismatch = -1;
                    scoring.m_gap_open = -1;
                    scoring.m_gap_ext = -1;

                    aln::GotohAligner<TYPE, aln::SimpleGotohScheme> aligner2(scoring);

                    return priv::banded::gotoh_alignment_score_dispatch<BAND_LEN, TYPE>::run(aligner2.scheme, pattern, quals, text, 0u, pattern.length(), 0u, min_score, context, sink, wfa);*/
                    
                    return banded_alignment_score<BAND_LEN>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        pattern,
                        quals,
                        text,
                        min_score,
                        sink,
                        wfa);
                }
            }          
            

            ///
            /// Calculate a window of the banded alignment matrix between a pattern and a text strings,
            /// using the Smith-Waterman algorithm. A checkpoint is used to pass the initial row
            /// and store the final one at the end of the window.
            ///
            /// \param pattern      shorter string (horizontal)
            /// \param quals        qualities string
            /// \param text         longer string (vertical)
            /// \param pos          offset in the reference string
            /// \param sink         output alignment sink
            ///
            /// \return             false if the minimum score was not reached, true otherwise
            ///
            template <
                uint32 BAND_LEN,
                AlignmentType TYPE,
                typename scoring_type,
                typename pattern_type,
                typename qual_type,
                typename text_type,
                typename sink_type,
                typename checkpoint_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE bool banded_alignment_score(
                const WfahAligner<TYPE, scoring_type> &aligner,
                pattern_type pattern,
                qual_type quals,
                text_type text,
                const int32 min_score,
                const uint32 window_begin,
                const uint32 window_end,
                sink_type &sink,
                checkpoint_type checkpoint,
                wfa_type& wfa)
            {
#ifdef         WFA_TESTS
                {
                    priv::banded::WfahCheckpointedScoringContext<BAND_LEN, TYPE, checkpoint_type> context(checkpoint);

                    return priv::banded::wfah_alignment_score_dispatch<BAND_LEN, TYPE>::run(aligner.scheme, pattern, quals, text, window_begin, window_end, 0u, min_score, context, sink, wfa);
                }
#else
                {
                    return banded_alignment_score<BAND_LEN>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        pattern,
                        quals,
                        text,
                        min_score,
                        window_begin,
                        window_end,
                        sink,
                        checkpoint,
                        wfa);
                }
#endif
            }
            
            ///
            /// Calculate the banded Smith-Waterman between a pattern and a text string
            /// while saving "checkpoints" along the way, i.e. saving one band
            /// of the DP matrix every CHECKPOINTS rows.
            ///
            /// \tparam BAND_LEN            size of the DP band
            ///
            /// \tparam CHECKPOINTS         number of DP rows between each checkpoint
            ///
            /// \tparam checkpoint_type     a class to represent the collection of checkpoints,
            ///                             represented as a linear array storing each checkpointed
            ///                             band contiguously.
            ///                             The class has to provide the non-const indexing operator[].
            ///
            /// \param pattern              pattern to be aligned
            /// \param quals                qualities string
            /// \param text                 text to align the pattern to
            /// \param pos                  offset in the text where the pattern should start aligning
            /// \param min_score            minimum tolerated score
            /// \param sink                 output alignment sink
            /// \param checkpoints          output checkpoints
            ///
            template <
                uint32 BAND_LEN,
                uint32 CHECKPOINTS,
                AlignmentType TYPE,
                typename scoring_type,
                typename pattern_type,
                typename qual_type,
                typename text_type,
                typename sink_type,
                typename checkpoint_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE bool banded_alignment_checkpoints(
                const WfahAligner<TYPE, scoring_type> &aligner,
                pattern_type pattern,
                qual_type quals,
                text_type text,
                const int32 min_score,
                sink_type &sink,
                checkpoint_type checkpoints,
                wfa_type& wfa)
            {
#ifdef         WFA_TESTS
                {
                    priv::banded::WfahCheckpointContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type> context(checkpoints);

                    return priv::banded::wfah_alignment_score_dispatch<BAND_LEN, TYPE>::run(aligner.scheme, pattern, quals, text, 0u, pattern.length(), 0u, min_score, context, sink, wfa);
                }
#else
                {
                    return banded_alignment_checkpoints<BAND_LEN, CHECKPOINTS>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        pattern,
                        quals,
                        text,
                        min_score,
                        sink,
                        checkpoints,
                        wfa);
                }
#endif
            }
            
            
            ///
            /// Compute the banded Dynamic Programming submatrix between two given checkpoints,
            /// storing its flow at each cell.
            /// The function returns the submatrix height.
            ///
            /// \tparam BAND_LEN            size of the DP band
            ///
            /// \tparam CHECKPOINTS         number of DP rows between each checkpoint
            ///
            /// \tparam checkpoint_type     a class to represent the collection of checkpoints,
            ///                             represented as a linear array storing each checkpointed
            ///                             band contiguously.
            ///                             The class has to provide the const indexing operator[].
            ///
            /// \tparam submatrix_type      a class to store the flow H, E and F submatrix, represented
            ///                             as a linear array of size (BAND_LEN*CHECKPOINTS).
            ///                             The class has to provide the non-const indexing operator[].
            ///                             Note that the H submatrix entries can assume only 3 values,
            ///                             while the E and F only 2 - hence the aggregate needs 4 bits
            ///                             per cell.
            ///
            /// \param pattern              pattern to be aligned
            /// \param quals                pattern quality scores
            /// \param text                 text to align the pattern to
            /// \param pos                  offset in the text where the pattern should start aligning
            /// \param min_score            minimum tolerated score
            /// \param checkpoints          precalculated checkpoints
            /// \param checkpoint_id        index of the first checkpoint defining the submatrix location
            /// \param submatrix            output flow submatrix
            ///
            /// \return                     submatrix height
            ///
            template <
                uint32 BAND_LEN,
                uint32 CHECKPOINTS,
                AlignmentType TYPE,
                typename scoring_type,
                typename pattern_string,
                typename qual_string,
                typename text_string,
                typename checkpoint_type,
                typename submatrix_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                uint32
                banded_alignment_submatrix(
                    const WfahAligner<TYPE, scoring_type> &aligner,
                    pattern_string pattern,
                    qual_string quals,
                    text_string text,
                    const int32 min_score,
                    checkpoint_type checkpoints,
                    const uint32 checkpoint_id,
                    submatrix_type submatrix,
                    wfa_type& wfa)
            {
#ifdef         WFA_TESTS
                {
                    priv::banded::WfahSubmatrixContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type, submatrix_type> context(checkpoints, checkpoint_id, submatrix);
                    const uint32 window_begin = checkpoint_id * CHECKPOINTS;
                    const uint32 window_end = nvbio::min(window_begin + CHECKPOINTS, uint32(pattern.length()));
                    NullSink sink;

                    priv::banded::wfah_alignment_score_dispatch<BAND_LEN, TYPE>::run(aligner.scheme, pattern, quals, text, window_begin, window_end, 0u, min_score, context, sink, wfa);
                    return window_end - window_begin;
                }
#else
                {
                      return banded_alignment_submatrix<BAND_LEN, CHECKPOINTS>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        pattern,
                        quals,
                        text,
                        min_score,
                        checkpoints,
                        checkpoint_id,
                        submatrix,
                        wfa);
                }
#endif
            }

            ///
            /// Given the Dynamic Programming submatrix between two checkpoints,
            /// backtrack from a given destination cell.
            /// The function returns the resulting source cell.
            ///
            /// \tparam BAND_LEN            size of the DP band
            ///
            /// \tparam CHECKPOINTS         number of DP rows between each checkpoint
            ///
            /// \tparam checkpoint_type     a class to represent the collection of checkpoints,
            ///                             represented as a linear array storing each checkpointed
            ///                             band contiguously.
            ///                             The class has to provide the const indexing operator[].
            ///
            /// \tparam submatrix_type      a class to store the flow submatrix, represented
            ///                             as a linear array of size (BAND_LEN*CHECKPOINTS).
            ///                             The class has to provide the const indexing operator[].
            ///                             Note that the submatrix entries can assume only 3 values,
            ///                             and could hence be packed in 2 bits.
            ///
            /// \tparam backtracer_type     a class to store the resulting list of backtracking operations.
            ///                             A model of \ref Backtracer.
            ///
            /// \param checkpoints          precalculated checkpoints
            /// \param checkpoint_id        index of the first checkpoint defining the DP submatrix,
            ///                             storing all bands between checkpoint_id and checkpoint_id+1.
            /// \param submatrix            precalculated flow submatrix
            /// \param submatrix_height     submatrix height
            /// \param sink                 in/out sink of the DP solution
            /// \param backtracer           backtracking output handler
            /// \param state                which matrix (H/E/F) the backtracker is at, starting where it ended
            ///                             from the previous checkpoint/submatrix
            ///
            /// \return                     true if the alignment source has been found, false otherwise
            ///

            template <
                uint32 BAND_LEN,
                uint32 CHECKPOINTS,
                AlignmentType TYPE,
                typename scoring_type,
                typename checkpoint_type,
                typename submatrix_type,
                typename backtracer_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE bool banded_alignment_traceback(
                const WfahAligner<TYPE, scoring_type> &aligner,
                checkpoint_type checkpoints,
                const uint32 checkpoint_id,
                submatrix_type submatrix,
                const uint32 submatrix_height,
                uint8 &state,
                uint2 &sink,
                backtracer_type &backtracer,
                wfa_type& wfa)
            {
 #ifdef         WFA_TESTS
                {
                    int32 current_entry = sink.x;                             //- sink.y;
                    int32 current_row = sink.y - checkpoint_id * CHECKPOINTS; // - 1u;

                    /*NVBIO_CUDA_DEBUG_ASSERT(current_entry >= 0 &&
                                                current_entry < (int32)BAND_LEN,
                                            "sw::banded_alignment_backtrack(): sink (%u,%u) -> local x coordinate %d not in [0,%d[\n", sink.x, sink.y, current_entry, BAND_LEN);
                    NVBIO_CUDA_DEBUG_ASSERT(current_row >= 0, "sw::banded_alignment_backtrack(): sink (%u,%u) -> local y coordinate %d not in [0,%u[\n", sink.x, sink.y, current_row, submatrix_height);
                    NVBIO_CUDA_DEBUG_ASSERT(current_row < (int32)submatrix_height, "sw::banded_alignment_backtrack(): sink (%u,%u) -> local y coordinate %d not in [0,%u[\n", sink.x, sink.y, current_row, submatrix_height);
    */
                    // TODO: enlever && current_entry >= 0
                    while (current_row >= 0 && current_entry > 0)
                    {
                        const int32 submatrix_cell = current_row * BAND_LEN + current_entry;
                        NVBIO_CUDA_DEBUG_ASSERT(submatrix_cell >= 0, "submatrix_cell error submatrix_cell=%d current_row=%d BAND_LEN=%d current_entry=%d !\n", submatrix_cell, current_row, BAND_LEN, current_entry);
                        const uint8 op = submatrix[submatrix_cell];
                        const uint8 h_op = op & HMASK;

                        if (TYPE == LOCAL)
                        {
                            if (state == HSTATE && h_op == SINK)
                            {
                                sink.y = current_row + checkpoint_id * CHECKPOINTS + 1u;
                                sink.x = current_entry + sink.y;
                                return true;
                            }
                        }

                        if (state == ESTATE)
                        {
                            if ((op & INSERTION_EXT) == 0u)
                                state = HSTATE;
                            --current_entry;
                            if (current_entry > 0 && current_row >= 0)
                                backtracer.push(DELETION);
                        }
                        else if (state == FSTATE)
                        {
                            if ((op & DELETION_EXT) == 0u)
                                state = HSTATE;
                            --current_entry;
                            --current_row;
                            if (current_row >= 0)
                                backtracer.push(INSERTION);
                        }
                        else
                        {
                            if (h_op == DELETION)
                                state = ESTATE;
                            else if (h_op == INSERTION)
                                state = FSTATE;
                            else
                            {
                                --current_row;
                                if (current_row >= 0)
                                    backtracer.push(SUBSTITUTION);
                            }
                        }

                        // NVBIO_CUDA_ASSERT(current_entry >= 0 && current_entry < BAND_LEN);
                    }

                    sink.y = checkpoint_id * CHECKPOINTS;
                    sink.x = current_entry + sink.y;

                    return false;
                }
#else
                {                    
                    return banded_alignment_traceback<BAND_LEN, CHECKPOINTS>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        checkpoints,
                        checkpoint_id,
                        submatrix,
                        submatrix_height,
                        state,
                        sink,
                        backtracer,
                        wfa);
                }
#endif         
            }

            /// @} // end of private group

        } // namespace priv

    } // namespace sw
} // namespace nvbio

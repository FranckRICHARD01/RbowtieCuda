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

#include <nvbio/alignment/sink.h>
#include <nvbio/alignment/utils.h>
#include <nvbio/alignment/alignment_base_inl.h>
#include <nvbio/basic/iterator.h>
#include <nvbio/strings/vectorized_string.h>
#include <nvbio/alignment/wfa.h>
#include <iostream>
#include <string>

namespace nvbio
{
    namespace aln
    {

        // ----------------------------- Wfah functions ---------------------------- //

        // #define NVBIO_SW_VECTOR_LOADING_WFA

        namespace priv
        {

            template <typename string_type>
            struct wfah_use_vectorization
            {
                static const bool VALUE = false;
            };
            template <typename T>
            struct wfah_use_vectorization<vector_view<T *>>
            {
                static const bool VALUE = true;
            };
            template <typename T>
            struct wfah_use_vectorization<vector_view<const T *>>
            {
                static const bool VALUE = true;
            };

            //
            // A helper scoring context class, which can be used to adapt the basic
            // algorithm to various situations, such as:
            //   scoring
            //   scoring within a window (i.e. saving only the last band within the window)
            //   computing checkpoints
            //   computing a flow submatrix
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, typename algorithm_tag>
            struct WfahScoringContext
            {
                // initialize the j-th column of the DP matrix
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                // \param scoring   scoring scheme
                // \param zero      zero value
                // \param infimum   infimum value
                template <typename column_type, typename scoring_type, typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                    const uint32 j,
                    const uint32 N,
                    column_type &column,
                    const scoring_type &scoring,
                    const score_type zero,
                    const score_type infimum)
                {
                    if (j == 0) // ensure this context can be used for multi-pass scoring
                    {
                        for (uint32 i = 0; i < N; ++i)
                        {
                            column[i] = 0;
                            //column[i].x = equal<algorithm_tag, PatternBlockingTag>() ? TYPE == GLOBAL ? scoring.text_gap_open() + scoring.text_gap_extension() * i : zero : TYPE != LOCAL ? scoring.text_gap_open() + scoring.text_gap_extension() * i
                            //                                                                                                                                                              : zero;
                            //column[i].y = TYPE == LOCAL ? zero : infimum;
                        }
                    }
                }

                // do something with the previous column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_column(
                    const uint32 j,
                    const uint32 N,
                    const column_type column) {}

                // do something with the last column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_column(
                    const uint32 j,
                    const uint32 M,
                    const uint32 N,
                    const column_type column) {}

                // do something with the newly computed cell
                //
                // \param i         row index
                // \param N         number of rows (column size)
                // \param j         column index
                // \param M         number of columns (row size)
                // \param score     computed score
                // \param dir       direction flow
                template <typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                    const uint32 i,
                    const uint32 j,
                    const score_type score,
                    const DirectionVector dir,
                    const DirectionVector edir,
                    const DirectionVector fdir) {}

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

            //
            // A helper checkpointed-scoring context class which allows to perform scoring in multiple
            // passes, saving & restoring a checkpoint each time.
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, typename algorithm_tag, typename checkpoint_type>
            struct WfahCheckpointedScoringContext
            {
                // constructor
                //
                // \param checkpoints       input checkpoints array
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                WfahCheckpointedScoringContext(checkpoint_type checkpoint) : m_checkpoint(checkpoint) {}

                // initialize the j-th column of the DP matrix
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                // \param scoring   scoring scheme
                // \param zero      zero value
                // \param infimum   infimum value
                template <typename column_type, typename scoring_type, typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                    const uint32 j,
                    const uint32 N,
                    column_type &column,
                    const scoring_type &scoring,
                    const score_type zero,
                    const score_type infimum)
                {
                    if (j == 0)
                    {
                        for (uint32 i = 0; i < N; ++i)
                        {
                            column[i].y = 0;
                            
                            //column[i].x = equal<algorithm_tag, PatternBlockingTag>() ? TYPE == GLOBAL ? scoring.text_gap_open() + scoring.text_gap_extension() * i : zero : TYPE != LOCAL ? scoring.text_gap_open() + scoring.text_gap_extension() * i
                            //                                                                                                                                                              : zero;
                            //column[i].y = TYPE == LOCAL ? zero : infimum;
                        }
                    }
                    else
                    {
                        for (uint32 i = 0; i < N; ++i)
                        {
                            column[i] = m_checkpoint[i];
                            //column[i].x = m_checkpoint[i].x;
                            //column[i].y = m_checkpoint[i].y;
                        }
                    }
                }

                // do something with the previous column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_column(
                    const uint32 j,
                    const uint32 N,
                    const column_type column) {}

                // do something with the last column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_column(
                    const uint32 j,
                    const uint32 M,
                    const uint32 N,
                    const column_type column)
                {
                    for (uint32 i = 0; i < N; ++i)
                    {
                        m_checkpoint[i].x = column[i];
                        //m_checkpoint[i].x = column[i].x;
                        //m_checkpoint[i].y = column[i].y;
                    }
                }

                // do something with the newly computed cell
                //
                // \param i         row index
                // \param N         number of rows (column size)
                // \param j         column index
                // \param M         number of columns (row size)
                // \param score     computed score
                // \param dir       direction flow
                template <typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                    const uint32 i,
                    const uint32 N,
                    const uint32 j,
                    const uint32 M,
                    const score_type score,
                    const DirectionVector dir,
                    const DirectionVector edir,
                    const DirectionVector fdir) {}

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

                checkpoint_type m_checkpoint;
            };

            //
            // A helper scoring context class, instantiated to keep track of checkpoints
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, uint32 CHECKPOINTS, typename checkpoint_type>
            struct WfahCheckpointContext
            {
                // constructor
                //
                // \param checkpoints       input checkpoints array
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                WfahCheckpointContext(checkpoint_type checkpoints) : m_checkpoints(checkpoints) {}

                // initialize the j-th column of the DP matrix
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                // \param scoring   scoring scheme
                // \param zero      zero value
                // \param infimum   infimum value
                template <typename column_type, typename scoring_type, typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                    const uint32 j,
                    const uint32 N,
                    column_type &column,
                    const scoring_type &scoring,
                    const score_type zero,
                    const score_type infimum)
                {
                    for (uint32 i = 0; i < N; ++i)
                    {
                        column[i] = 0;                        
                        //column[i].x = TYPE == GLOBAL ? scoring.text_gap_open() + scoring.text_gap_extension() * i : zero;
                        //column[i].y = TYPE == LOCAL ? zero : infimum;
                    }
                }

                // do something with the previous column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_column(
                    const uint32 j,
                    const uint32 N,
                    const column_type column)
                {
                    typedef typename std::iterator_traits<column_type>::value_type vector_type;
                    typedef typename vector_traits<vector_type>::value_type value_type;

                    // save checkpoint
                    if ((j & (CHECKPOINTS - 1)) == 0u)
                    {
                        const uint32 checkpoint_id = j / CHECKPOINTS;

                        for (uint32 i = 0; i < N; ++i)
                            m_checkpoints[checkpoint_id * N + i] = 0;//make_vector<value_type>(column[i].x, column[i].y);
                    }
                }

                // do something with the last column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_column(
                    const uint32 j,
                    const uint32 M,
                    const uint32 N,
                    const column_type column) {}

                // do something with the newly computed cell
                //
                // \param i         row index
                // \param N         number of rows (column size)
                // \param j         column index
                // \param M         number of columns (row size)
                // \param score     computed score
                // \param dir       direction flow
                template <typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                    const uint32 i,                  
                    const uint32 j,                  
                    const score_type score,
                    const DirectionVector dir,
                    const DirectionVector edir,
                    const DirectionVector fdir) {}

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

            //
            // A helper scoring context class, instantiated to keep track of the direction vectors
            // of a DP submatrix between given checkpoints
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, uint32 CHECKPOINTS, typename checkpoint_type, typename submatrix_type>
            struct WfahSubmatrixContext
            {
                // constructor
                //
                // \param checkpoints       input checkpoints array
                // \param checkpoint_id     id of the checkpoint defining the first column of the submatrix
                // \param submatrix         submatrix output storage
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                WfahSubmatrixContext(
                    const checkpoint_type checkpoints,
                    const uint32 checkpoint_id,
                    const submatrix_type submatrix) : m_checkpoints(checkpoints),
                                                      m_checkpoint_id(checkpoint_id),
                                                      m_submatrix(submatrix) {}

                // initialize the j-th column of the DP matrix
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                // \param scoring   scoring scheme
                // \param zero      zero value
                // \param infimum   infimum value
                template <typename column_type, typename scoring_type, typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void init(
                    const uint32 j,
                    const uint32 N,
                    column_type &column,
                    const scoring_type &scoring,
                    const score_type zero,
                    const score_type infimum)
                {
                    /*typedef typename std::iterator_traits<column_type>::value_type value_type;

                    const uint32 r = m_checkpoint_id * N;

                    // restore the checkpoint
                    for (uint32 i = 0; i < N; ++i)
                    {
                        const value_type f = m_checkpoints[r + i];
                        column[i] = f;
                        //column[i].y = f.x;
                        //column[i].y = f.y;
                    }*/
                }

                // do something with the previous column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void previous_column(
                    const uint32 j,
                    const uint32 N,
                    const column_type column) {}

                // do something with the last column
                //
                // \param j         column index
                // \param N         column size
                // \param column    column values
                template <typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void last_column(
                    const uint32 j,
                    const uint32 M,
                    const uint32 N,
                    const column_type column) {}

                // do something with the newly computed cell
                //
                // \param i         row index
                // \param N         number of rows (column size)
                // \param j         column index
                // \param M         number of columns (row size)
                // \param score     computed score
                // \param dir       direction flow
                template <typename score_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void new_cell(
                    const uint32 i,                    
                    const uint32 j,                   
                    const score_type score,
                    const DirectionVector hdir,
                    const DirectionVector edir,
                    const DirectionVector fdir)
                {
                    // save the direction vector
                    const uint32 offset = m_checkpoint_id * CHECKPOINTS;
                                       
                    if (TYPE == LOCAL)
                        m_submatrix[j * CHECKPOINTS + (i - offset)] = score ? hdir : SINK;
                    else
                        m_submatrix[j * CHECKPOINTS + (i - offset)] = hdir;
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

                            /*if (h<0 || h >WFA_BAND_LEN2_X)
                            {
                                char ref_str[1000];
                                char read_str[1000];

                                dna_ref_type_to_string(ref, N, ref_str);
                                dna_text_cache_to_string(q_cache, M, read_str);
                            }*/

                            NVBIO_CUDA_DEBUG_ASSERT(score >= 0 && score < 1000, "Problem score=%d\n", score);
                            //NVBIO_CUDA_DEBUG_ASSERT(h >= 0 && h < WFA_BAND_LEN2_X, "Problem h=%d M=%d N=%d score=%d\n", h, M, N, score);
                            NVBIO_CUDA_DEBUG_ASSERT(v >= 0 && v < WFA_BAND_LEN2_X, "Problem v=%d\n", v);

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

                                NVBIO_CUDA_DEBUG_ASSERT(pos >= 0 && pos < 1000, "Problem pos result SUBSTITUTION, pos=%d  score=%d last_score=%d h=%d v=%d\n", pos, score, last_s, h, v );
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
                                            y--;

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
                                    y--;

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
                                    y--;
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

            template <uint32 BAND_LEN, AlignmentType TYPE, typename algorithm_tag, typename symbol_type>
            struct wfah_alignment_score_dispatch
            {
            };

            //
            // A template struct used to possibly specialize the implementation of the Wfah-based alignment based on
            // the template parameters.
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, typename symbol_type>
            struct wfah_alignment_score_dispatch<BAND_LEN, TYPE, PatternBlockingTag, symbol_type>
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
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type update(
                    context_type &context,
                    const score_type block,
                    const score_type M,
                    const score_type N,
                    const ref_type ref,
                    const text_cache_type q_cache,
                    score_type &temp_i,
                    score_type *H_band,
                    score_type *F_band,
                    sink_type &sink,
                    const score_type min_score,
                    score_type &max_score,
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
                    
                    wfa.H_Band.set_null(s, false);
                    wfa.E_Band.set_null(s, true);
                    wfa.F_Band.set_null(s, true);
                    wfa.H_Band.set_hi_lo(0, true);
                    wfa.F_Band.set_hi_lo(0, false);
                    wfa.E_Band.set_hi_lo(0, false);
                    wfa.H_Band.set(0, 0, 0);
                    wfa.E_Band.set(0, 0, WFA_MIN1);
                    wfa.F_Band.set(0, 0, WFA_MIN1);
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

                //
                // Calculate the alignment score between a string and a reference, using the Smith-Waterman algorithm,
                // using a templated column storage.
                //
                // This function is templated over:
                //   - a context that is passed the computed DP matrix values, and can be
                //     used to specialize its behavior.
                //   - a sink that is used to report successful alignments
                //
                // Furthermore, the function can be called on a window of the pattern, assuming that the context
                // will provide the proper initialization for the first column of the corresponding DP matrix window.
                //
                // \param context       template context class, used to specialize the behavior of the aligner
                // \param query         input pattern (horizontal string)
                // \param quals         input pattern qualities (horizontal string)
                // \param ref           input text (vertical string)
                // \param scoring       scoring scheme
                // \param min_score     minimum output score
                // \param sink          alignment sink
                // \param window_begin  beginning of pattern window
                // \param window_end    end of pattern window
                //
                // \return              false if early-exited, true otherwise
                //
                
                template <
                    typename context_type,
                    typename query_type,
                    typename qual_type,
                    typename ref_type,
                    typename scoring_type,
                    typename sink_type,
                    typename column_type,
                    typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool run(
                    const scoring_type &scoring,
                    context_type &context,
                    query_type query,
                    qual_type quals,
                    ref_type ref,
                    const int32 min_score,
                    sink_type &sink,
                    const int32 window_begin,
                    const int32 window_end,
                    column_type temp,
                    wfa_type& wfa)
                {
                    //
                    // This function breaks the DP matrix in vertical stripes of BAND_LEN cells,
                    // so as to keep the current reduced row (a band of coefficients) in registers
                    // throughout the entire computa[BAND_LEN+1]

#ifdef WFA_TESTS
                    #define BAND_LEN2 1000
#else
                    #define BAND_LEN2 600
#endif
                    typedef int32 score_type;

                    score_type M = query.length();
                    score_type N = ref.length();
                    score_type NN = N;

                    score_type begin = 0;
                    

#ifndef WFA_TESTS
                    score_type end = N;
                    score_type lim = 20;
                    const score_type threshold = 5;

                    if (M >= 150)
                    {
#pragma unroll
                        for (score_type i = 0; i < N - M + lim; i++)
                        {
                            score_type num = 0;

#pragma unroll
                            for (score_type j = 0; j < lim; j++)
                            {
                                if (ref[i + j] != query[j])
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
                        for (score_type i = N - 1; i > M - lim; i--)
                        {
                            score_type num = 0;
#pragma unroll
                            for (score_type j = M - 1; j > M - 1 - lim; j--)
                            {
                                if (ref[i - (M - 1 - j)] != query[j])
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
                            N = nvbio::max(M, end - begin);
                            begin = nvbio::max(0, begin);
                        }
                    }
#endif

                    score_type H_band[BAND_LEN + 1];
                    score_type F_band[BAND_LEN + 1];

                    const score_type G_o = abs(scoring.pattern_gap_open());
                    const score_type G_e = abs(scoring.pattern_gap_extension());
                    const score_type S = abs(scoring.mismatch());
                    const score_type V = abs(scoring.match());
                    const score_type zero = score_type(0);
                    const score_type infimum = WFA_MIN;//Field_traits<temp_scalar_type>::min() - nvbio::min(G_o, G_e);

                    typedef uchar2 text_cache_type;
                   // typedef typename Reference_cache<BAND_LEN>::type ref_type;
                    unsigned char  q_cache_ref[BAND_LEN2];
                    text_cache_type q_cache[BAND_LEN2];                         

                    // initialize the first column
                    context.init(window_begin, N, temp, scoring, zero, infimum);

                    /*for (score_type i = 0; i < N; i++)
                    {
                        for (score_type j = 0; j < M; j++)
                        {
                            context.new_cell(
                                i,
                                j,
                                WFA_MIN,
                                SUBSTITUTION,
                                SUBSTITUTION,
                                SUBSTITUTION);
                        }
                    }*/

                    // wfa
#ifndef WFA_TESTS
                    wavefront_heuristic_set_wfadaptive(wfa, DIM_SHARED2 * 2, DIM_SHARED2 * 2, 1);
                    // wavefront_heuristic_set_wfmash(wfa, 15, DIM_SHARED2 / 2, 1);
                    // wavefront_heuristic_set_xdrop(wfa, 15, 1);
                    //wavefront_heuristic_set_zdrop(wfa, 15, 3);
                    // wavefront_heuristic_set_banded_static(wfa, -DIM_SHARED2, DIM_SHARED2);
                    wavefront_heuristic_set_banded_adaptive(wfa, -DIM_SHARED2, DIM_SHARED2, 1);
#endif

// initialize extra bands
#pragma unroll 1
                    for (int32 j = 0; j < WFA_BAND_LEN2_Y; j++)
                    {
                        wfa.Point_H_BAND[j] = 0;
                    }

                    //const uint32 end_block = (window_end == (uint32)M) ? nvbio::max(BAND_LEN, BAND_LEN * (((uint32)M + BAND_LEN - 1) / BAND_LEN)) : window_end + BAND_LEN;

                    score_type max_score = Field_traits<score_type>::min();

#pragma unroll 1
                    for (int32 t = 0; t < window_end - window_begin; t++)
                    {
                        q_cache[t] = make_uchar2(query[t + window_begin], quals[t + window_begin]);
                    }

                    const int32 MM = window_end - window_begin;

                    // load a block of entries from eah query
#pragma unroll 1
                    for (int32 t = 0; t < N; t++)
                    {
                        /*if (window_begin + t < N)
                        {
                            q_cache_ref[t] = ref[window_begin + t];
                        }
                        else
                        {
                            q_cache_ref[t] = 0;
                        }*/

                        if (t + begin < NN) q_cache_ref[t] = ref[t + begin];
                    }     

                   
                    score_type temp_i = 0;

                    max_score = update<true>(
                        context,
                        0,
                        MM,
                        N,
                        q_cache_ref,
                        q_cache,
                        temp_i,
                        H_band,
                        F_band,
                        sink,
                        min_score,
                        max_score,
                        V,
                        S,
                        G_o,
                        G_e,
                        zero,
                        scoring,
                        wfa);

                    // save the last column
                    // context.last_column(window_end, text_len, pattern_len, temp);

                    // if (TYPE == GLOBAL)
                    // save_Mth<BAND_LEN>( M, H_band, N-1, sink );

#ifndef WFA_TESTS
                    
                    lim = 30;
                    //bool test = false;

#pragma unroll
                    for (score_type i = NN - 1; i > M - lim; i--)
                    {
                        score_type num = 0;
#pragma unroll
                        for (score_type j = M - 1; j > M - 1 - lim; j--)
                        {
                            if (ref[i - (M - 1 - j)] != query[j])
                                num++;

                            if (num >= 5)
                                break;
                        }
                        if (num < 5)
                        {
                            // num_sauv2 = num;
                            NN = i;
                            //test = true;
                            break;
                        }
                    }

                    /*if (!test)
                    {
                        NN = 500;
                        max_score = 255;
                    }*/
#endif

                    if (TYPE == GLOBAL)
                        sink.report(-max_score, make_uint2(N, MM));
#ifdef WFA_TESTS
                    if (TYPE == SEMI_GLOBAL)
                        sink.report(-max_score, make_uint2(N, MM));
#else
                    if (TYPE == SEMI_GLOBAL)
                        sink.report(-max_score, make_uint2(NN, MM));
#endif

                        //  sink.report( m, make_uint2( i+1, M ) );

                    return true;
                }

                //
                // Calculate the alignment score between a string and a reference, using the Smith-Waterman algorithm,
                // using local memory storage for the boundary columns.
                //
                // This function is templated over:
                //   1. a context that is passed the computed DP matrix values, and can be
                //      used to specialize its behavior.
                //   2. a sink that is used to report successful alignments
                //
                // Furthermore, the function can be called on a window of the pattern, assuming that the context
                // will provide the proper initialization for the first column of the corresponding DP matrix window.
                //
                // \param context       template context class, used to specialize the behavior of the aligner
                // \param query         input pattern (horizontal string)
                // \param quals         input pattern qualities (horizontal string)
                // \param ref           input text (vertical string)
                // \param scoring       scoring scheme
                // \param min_score     minimum output score
                // \param sink          alignment sink
                // \param window_begin  beginning of pattern window
                // \param window_end    end of pattern window
                //
                // \return              false if early-exited, true otherwise
                //
                template <
                    uint32 MAX_REF_LEN,
                    typename context_type,
                    typename query_type,
                    typename qual_type,
                    typename ref_type,
                    typename scoring_type,
                    typename sink_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool run(
                    const scoring_type &scoring,
                    context_type &context,
                    query_type query,
                    qual_type quals,
                    ref_type ref,
                    const int32 min_score,
                    sink_type &sink,
                    const uint32 window_begin,
                    const uint32 window_end)
                {
                    // instantiated a local memory array
                    short2 temp[MAX_REF_LEN];
                    short2 *temp_ptr = temp;

                    return run(
                        context,
                        query,
                        quals,
                        ref,
                        scoring,
                        min_score,
                        sink,
                        window_begin,
                        window_end,
                        temp_ptr);
                }
            };

            //
            // A template struct used to possibly specialize the implementation of the Wfah-based alignment based on
            // the template parameters.
            //
            template <uint32 BAND_LEN, AlignmentType TYPE, typename symbol_type>
            struct wfah_alignment_score_dispatch<BAND_LEN, TYPE, TextBlockingTag, symbol_type>
            {
                // update a DP row
                //
                template <
                    bool CHECK_N,
                    typename context_type,
                    typename ref_cache,
                    typename score_type,
                    typename temp_iterator,
                    typename sink_type,
                    typename scoring_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void update_row(
                    context_type &context,
                    const uint32 block,
                    const uint32 N,
                    const uint32 i,
                    const uint32 M,
                    const symbol_type q_i,
                    const uint8 qq_i,
                    // const score_type     m_i,
                    const ref_cache r_cache,
                    temp_iterator temp,
                    score_type &temp_i,
                    score_type *H_band,
                    score_type *F_band,
                    sink_type &sink,
                    const score_type min_score,
                    score_type &max_score,
                    const score_type G_o,
                    const score_type G_e,
                    const score_type zero,
                    const scoring_type scoring)
                {
                    typedef typename std::iterator_traits<temp_iterator>::value_type temp_cell_type;
                    typedef typename vector_traits<temp_cell_type>::value_type temp_scalar_type;

                    //
                    // NOTE:
                    // It might look as if we were going to make lots of out-of-bounds accesses here,
                    // as we loop across BAND_LEN cells irrespectively of whether the band straddles
                    // the end of the pattern.
                    // However, these accesses don't cause any page faults as they refer to properly
                    // sized temporary arrays (that in practice are placed in registers), and don't
                    // affect the results as they contribute to unused portions of the DP matrices.
                    // The reporting of the scores is properly guarded.
                    //

                    // set the 0-th coefficient in the band to be equal to the (i-1)-th row of the left column (diagonal term)
                    score_type H_diag = temp_i;

                    H_band[0] = temp_i = temp[i].x;
                    score_type E = temp[i].y;

#pragma unroll
                    for (uint32 j = 1; j <= BAND_LEN; ++j)
                    {
                        // update F
                        const score_type ftop = F_band[j] + G_e;
                        const score_type htop = H_band[j] + G_o;
                        F_band[j] = nvbio::max(ftop, htop);
                        const DirectionVector fdir = ftop > htop ? INSERTION_EXT : SUBSTITUTION;

                        // update E
                        const score_type eleft = E + G_e;
                        const score_type hleft = H_band[j - 1] + G_o;
                        E = nvbio::max(eleft, hleft);
                        const DirectionVector edir = eleft > hleft ? DELETION_EXT : SUBSTITUTION;

                        const symbol_type r_j = r_cache[j - 1];
                        // const score_type S_ij     = (r_j == q_i) ? m_i : scoring.mismatch( r_j, q_i, qq_i );
                        const score_type S_ij = scoring.substitution(block + j, i, r_j, q_i, qq_i);
                        const score_type diagonal = H_diag + S_ij;
                        const score_type top = F_band[j];
                        const score_type left = E;
                        score_type hi = nvbio::max3(left, top, diagonal);
                        hi = (TYPE == LOCAL) ? nvbio::max(hi, zero) : hi; // clamp to zero
                        H_diag = H_band[j];
                        H_band[j] = hi;

                        if ((CHECK_N == false) || (block + j <= N))
                        {
                            context.new_cell(
                                i, M,
                                block + j - 1, N,
                                hi,
                                top > left ? (top > diagonal ? INSERTION : SUBSTITUTION) : (left > diagonal ? DELETION : SUBSTITUTION),
                                edir,
                                fdir);
                        }
                    }

                    // save the last entry of the band
                    temp[i] = make_vector<temp_scalar_type>(H_band[BAND_LEN], E);

                    NVBIO_CUDA_ASSERT(H_band[BAND_LEN] >= Field_traits<temp_scalar_type>::min());
                    NVBIO_CUDA_ASSERT(H_band[BAND_LEN] <= Field_traits<temp_scalar_type>::max());
                    NVBIO_CUDA_ASSERT(E >= Field_traits<temp_scalar_type>::min());
                    NVBIO_CUDA_ASSERT(E <= Field_traits<temp_scalar_type>::max());

                    max_score = nvbio::max(max_score, H_band[BAND_LEN]);

                    if (TYPE == LOCAL)
                    {
                        if (CHECK_N)
                        {
                            // during local alignment we save the best score across all bands
                            for (uint32 j = 1; j <= BAND_LEN; ++j)
                            {
                                if (block + j <= N)
                                    sink.report(H_band[j], make_uint2(block + j, i + 1));
                            }
                        }
                        else
                        {
                            // during local alignment we save the best score across all bands
                            for (uint32 j = 1; j <= BAND_LEN; ++j)
                                sink.report(H_band[j], make_uint2(block + j, i + 1));
                        }
                    }
                }

                //
                // Calculate the alignment score between a string and a reference, using the Smith-Waterman algorithm,
                // using a templated column storage.
                //
                // This function is templated over:
                //   - a context that is passed the computed DP matrix values, and can be
                //     used to specialize its behavior.
                //   - a sink that is used to report successful alignments
                //
                // Furthermore, the function can be called on a window of the pattern, assuming that the context
                // will provide the proper initialization for the first column of the corresponding DP matrix window.
                //
                // \param context       template context class, used to specialize the behavior of the aligner
                // \param query         input pattern (horizontal string)
                // \param quals         input pattern qualities (horizontal string)
                // \param ref           input text (vertical string)
                // \param scoring       scoring scheme
                // \param min_score     minimum output score
                // \param sink          alignment sink
                // \param window_begin  beginning of pattern window
                // \param window_end    end of pattern window
                //
                // \return              false if early-exited, true otherwise
                //
                template <
                    typename context_type,
                    typename query_type,
                    typename qual_type,
                    typename ref_type,
                    typename scoring_type,
                    typename sink_type,
                    typename column_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool run(
                    const scoring_type &scoring,
                    context_type &context,
                    query_type query,
                    qual_type quals,
                    ref_type ref,
                    const int32 min_score,
                    sink_type &sink,
                    const uint32 window_begin,
                    const uint32 window_end,
                    column_type temp)
                {
                    //
                    // This function breaks the DP matrix in vertical stripes of BAND_LEN cells,
                    // so as to keep the current reduced row (a band of coefficients) in registers
                    // throughout the entire computation.
                    // Within each stripe, the matrix is updated top-to-bottom, and the right-most
                    // border of the stripe is saved to a local memory array so as to allow resuming
                    // its values when processing the next stripe.
                    //
                    const uint32 M = query.length();
                    const uint32 N = ref.length();

                    typedef int32 score_type;
                    score_type H_band[BAND_LEN + 1];
                    score_type F_band[BAND_LEN + 1];

                    typedef typename std::iterator_traits<column_type>::value_type temp_cell_type;
                    typedef typename vector_traits<temp_cell_type>::value_type temp_scalar_type;

                    const score_type G_o = scoring.pattern_gap_open();
                    const score_type G_e = scoring.pattern_gap_extension();
                    const score_type zero = score_type(0);
                    const score_type infimum = Field_traits<temp_scalar_type>::min() - nvbio::min(G_o, G_e);

                    uint8 r_cache[BAND_LEN];

                    // initialize the first column
                    context.init(window_begin, M, temp, scoring, zero, infimum);

                    const uint32 end_block = (window_end == N) ? nvbio::max(BAND_LEN, BAND_LEN * ((N + BAND_LEN - 1) / BAND_LEN)) : window_end + BAND_LEN;

                    // loop across the short edge of the DP matrix (i.e. the columns)
                    for (uint32 block = window_begin; block + BAND_LEN < end_block; block += BAND_LEN)
                    {
                        // save the previous column
                        context.previous_column(block, M, temp);

// load a block of entries from the reference
#pragma unroll
                        for (uint32 t = 0; t < BAND_LEN; ++t)
                            r_cache[t] = ref[block + t];

// initialize the first band
#pragma unroll
                        for (uint32 j = 0; j <= BAND_LEN; ++j)
                        {
                            H_band[j] = (TYPE == GLOBAL) ? (block + j > 0 ? G_o + G_e * (block + j - 1u) : zero) : zero;
                            F_band[j] = infimum;
                        }

                        score_type max_score = Field_traits<score_type>::min();

                        score_type temp_i = H_band[0];

                        // loop across the short edge of the DP matrix (i.e. the rows)
                        for (uint32 i = 0; i < M; ++i)
                        {
                            // load the new character from the query
                            const uint8 q_i = query[i];
                            const uint8 qq_i = quals[i];

                            // const int32 m_i = scoring.match(qq_i);

                            update_row<false>(
                                context,
                                block, N,
                                i, M,
                                q_i,
                                qq_i,
                                // m_i,
                                r_cache,
                                temp,
                                temp_i,
                                H_band,
                                F_band,
                                sink,
                                min_score,
                                max_score,
                                G_o, G_e,
                                zero,
                                scoring);
                        }
                        if (TYPE == SEMI_GLOBAL)
                        {
                            // during semi-global alignment we save the best score across the last row
                            for (uint32 j = 1; j <= BAND_LEN; ++j)
                                sink.report(H_band[j], make_uint2(block + j, M));
                        }

                        // we are now (N - block - BAND_LEN) columns from the last one: check whether
                        // we could theoretically reach the minimum score
                        const score_type missing_cols = score_type(N - block - BAND_LEN);
                        if (max_score + missing_cols * scoring.match(255) < score_type(min_score))
                            return false;
                    }

                    // process the very last stripe
                    if (window_end == N)
                    {
                        const uint32 block = end_block - BAND_LEN;

                        // save the previous column
                        context.previous_column(block, M, temp);

                        // load a block of entries from each query
                        const uint32 block_end = nvbio::min(block + BAND_LEN, N);
#pragma unroll
                        for (uint32 t = 0; t < BAND_LEN; ++t)
                        {
                            if (block + t < block_end)
                                r_cache[t] = ref[block + t];
                        }

// initialize the first band
#pragma unroll
                        for (uint32 j = 0; j <= BAND_LEN; ++j)
                        {
                            H_band[j] = (TYPE == GLOBAL) ? (block + j > 0 ? G_o + G_e * (block + j - 1u) : zero) : zero;
                            F_band[j] = infimum;
                        }

                        score_type max_score = Field_traits<score_type>::min();

                        score_type temp_i = H_band[0];

#if defined(NVBIO_SW_VECTOR_LOADING_WFA)
                        // check whether we should use vectorized loads
                        if (wfah_use_vectorization<query_type>::VALUE)
                        {
                            //
                            // loop across the long edge of the DP matrix (i.e. the rows)
                            //
                            const uint32 QUERY_VECTOR_WIDTH = vectorized_string<query_type>::VECTOR_WIDTH;

                            const uint2 vec_range = vectorized_string_range(query);

                            for (uint32 i = 0; i < vec_range.x; ++i)
                            {
                                // load the new character from the query
                                const uint8 q_i = query[i];
                                const uint8 qq_i = quals[i];

                                update_row<true>(
                                    context,
                                    block, N,
                                    i, M,
                                    q_i,
                                    qq_i,
                                    r_cache,
                                    temp,
                                    temp_i,
                                    H_band,
                                    F_band,
                                    sink,
                                    min_score,
                                    max_score,
                                    G_o, G_e,
                                    zero,
                                    scoring);
                            }
                            for (uint32 i = vec_range.x; i < vec_range.y; i += QUERY_VECTOR_WIDTH)
                            {
                                uint8 q[QUERY_VECTOR_WIDTH];
                                uint8 qq[QUERY_VECTOR_WIDTH];

                                // load QUERY_VECTOR_WIDTH new characters from the query
                                vectorized_string_load(query, i, q);

                                for (uint32 j = 0; j < QUERY_VECTOR_WIDTH; ++j)
                                    qq[j] = quals[i + j];

                                for (uint32 j = 0; j < QUERY_VECTOR_WIDTH; ++j)
                                {
                                    // load the new character from the reference
                                    const uint8 q_i = q[j];
                                    const uint8 qq_i = qq[j];

                                    update_row<true>(
                                        context,
                                        block, N,
                                        i + j, M,
                                        q_i,
                                        qq_i,
                                        r_cache,
                                        temp,
                                        temp_i,
                                        H_band,
                                        F_band,
                                        sink,
                                        min_score,
                                        max_score,
                                        G_o, G_e,
                                        zero,
                                        scoring);
                                }
                            }
                            for (uint32 i = vec_range.y; i < M; ++i)
                            {
                                // load the new character from the query
                                const uint8 q_i = query[i];
                                const uint8 qq_i = quals[i];

                                update_row<true>(
                                    context,
                                    block, N,
                                    i, M,
                                    q_i,
                                    qq_i,
                                    r_cache,
                                    temp,
                                    temp_i,
                                    H_band,
                                    F_band,
                                    sink,
                                    min_score,
                                    max_score,
                                    G_o, G_e,
                                    zero,
                                    scoring);
                            }
                        }
                        else
#endif
                        {
                            // typedef typename string_traits<query_type>::forward_iterator forward_query_iterator;
                            // typedef typename string_traits<qual_type>::forward_iterator  forward_qual_iterator;

                            // forward_query_iterator query_it( query.begin() );
                            // forward_qual_iterator  quals_it( quals.begin() );

                            //
                            // loop across the short edge of the DP matrix (i.e. the query)
                            //
                            for (uint32 i = 0; i < M; ++i)
                            {
                                // load the new character from the query
                                const uint8 q_i = query[i];
                                const uint8 qq_i = quals[i];
                                // const uint8 q_i  = *query_it; ++query_it;
                                // const uint8 qq_i = *quals_it; ++quals_it;

                                update_row<true>(
                                    context,
                                    block, N,
                                    i, M,
                                    q_i,
                                    qq_i,
                                    r_cache,
                                    temp,
                                    temp_i,
                                    H_band,
                                    F_band,
                                    sink,
                                    min_score,
                                    max_score,
                                    G_o, G_e,
                                    zero,
                                    scoring);
                            }
                        }

                        if (TYPE == SEMI_GLOBAL)
                        {
                            // during semi-global alignment we save the best score across the last row
                            for (uint32 j = 1; j <= BAND_LEN; ++j)
                            {
                                if (block + j <= N)
                                    sink.report(H_band[j], make_uint2(block + j, M));
                            }
                        }
                        else if (TYPE == GLOBAL)
                        {
                            // during global alignment we save the best score at cell [N][M]
                            for (uint32 j = 1; j <= BAND_LEN; ++j)
                            {
                                if (block + j == N)
                                    sink.report(H_band[j], make_uint2(block + j, M));
                            }
                        }
                    }

                    // save the last column
                    context.last_column(window_end, N, M, temp);
                    return true;
                }

                //
                // Calculate the alignment score between a string and a reference, using the Smith-Waterman algorithm,
                // using local memory storage for the boundary columns.
                //
                // This function is templated over:
                //   1. a context that is passed the computed DP matrix values, and can be
                //      used to specialize its behavior.
                //   2. a sink that is used to report successful alignments
                //
                // Furthermore, the function can be called on a window of the pattern, assuming that the context
                // will provide the proper initialization for the first column of the corresponding DP matrix window.
                //
                // \param context       template context class, used to specialize the behavior of the aligner
                // \param query         input pattern (horizontal string)
                // \param quals         input pattern qualities (horizontal string)
                // \param ref           input text (vertical string)
                // \param scoring       scoring scheme
                // \param min_score     minimum output score
                // \param sink          alignment sink
                // \param window_begin  beginning of pattern window
                // \param window_end    end of pattern window
                //
                // \return              false if early-exited, true otherwise
                //
                template <
                    uint32 MAX_PATTERN_LEN,
                    typename context_type,
                    typename string_type,
                    typename qual_type,
                    typename ref_type,
                    typename scoring_type,
                    typename sink_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool run(
                    const scoring_type &scoring,
                    context_type &context,
                    string_type query,
                    qual_type quals,
                    ref_type ref,
                    const int32 min_score,
                    sink_type &sink,
                    const uint32 window_begin,
                    const uint32 window_end)
                {
                    // instantiate a local memory array
                    short2 temp[MAX_PATTERN_LEN];
                    short2 *temp_ptr = temp;

                    return run(
                        context,
                        query,
                        quals,
                        ref,
                        scoring,
                        min_score,
                        sink,
                        window_begin,
                        window_end,
                        temp_ptr);
                }
            };

            template <AlignmentType TYPE, uint32 DIM, typename symbol_type>
            struct wfah_bandlen_selector
            {
                static const uint32 BAND_LEN = 16u / DIM;
            };

            template <AlignmentType TYPE, uint32 DIM>
            struct wfah_bandlen_selector<TYPE, DIM, simd4u8>
            {
#if __CUDA_ARCH__ >= 300
                static const uint32 BAND_LEN = 8u;
#else
                static const uint32 BAND_LEN = 1u;
#endif
            };

            //
            // Calculate the alignment score between a pattern and a text, using the Wfah algorithm.
            //
            // \param aligner      scoring scheme
            // \param pattern      pattern string (horizontal
            // \param quals        pattern qualities
            // \param text         text string (vertical)
            // \param min_score    minimum score
            //
            template <
                AlignmentType TYPE,
                typename scoring_type,
                typename algorithm_tag,
                typename pattern_string,
                typename qual_string,
                typename text_string,
                typename column_type>
            struct alignment_score_dispatch<
                WfahAligner<TYPE, scoring_type, algorithm_tag>,
                pattern_string,
                qual_string,
                text_string,
                column_type>
            {
                typedef WfahAligner<TYPE, scoring_type, algorithm_tag> aligner_type;

                /// dispatch scoring across the whole pattern
                ///
                ///
                /// \param aligner      scoring scheme
                /// \param pattern      pattern string (horizontal)
                /// \param quals        pattern qualities
                /// \param text         text string (vertical)
                /// \param min_score    minimum score
                /// \param sink         output alignment sink
                /// \param column       temporary column storage
                ///
                /// \return             true iff the minimum score was reached
                ///
                template <typename sink_type, typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool dispatch(
                    const aligner_type aligner,
                    const pattern_string pattern,
                    const qual_string quals,
                    const text_string text,
                    const int32 min_score,
                    sink_type &sink,
                    column_type column,
                    wfa_type& wfa)
                {
#ifdef WFA_TESTS
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = wfah_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::WfahScoringContext<BAND_LEN, TYPE, algorithm_tag> context;

                        const uint32 length = equal<algorithm_tag, PatternBlockingTag>() ? pattern.length() : text.length();

                        return wfah_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(aligner.scheme, context, pattern, quals, text, min_score, sink, 0, length, column, wfa);
                    }
#else
                    {              
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = sw_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::SWScoringContext<BAND_LEN, TYPE, algorithm_tag> context;

                        const uint32 length = equal<algorithm_tag, PatternBlockingTag>() ? pattern.length() : text.length();

                        return sw_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(EditDistanceSWScheme(), context, pattern, quals, text, min_score, sink, 0, length, column, wfa);
                    }
#endif
                }

                /// dispatch scoring in a window of the pattern
                ///
                /// \tparam checkpoint_type     a class to represent the checkpoint: an array of size equal to the text,
                ///                             that has to provide the const indexing operator[].
                ///
                /// \param aligner      scoring scheme
                /// \param pattern      pattern string (horizontal)
                /// \param quals        pattern qualities
                /// \param text         text string (vertical)
                /// \param min_score    minimum score
                /// \param sink         output alignment sink
                /// \param column       temporary column storage
                ///
                /// \return             true iff the minimum score was reached
                ///
                template <
                    typename sink_type,
                    typename checkpoint_type,
                    typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool dispatch(
                    const aligner_type aligner,
                    const pattern_string pattern,
                    const qual_string quals,
                    const text_string text,
                    const int32 min_score,
                    const uint32 window_begin,
                    const uint32 window_end,
                    sink_type &sink,
                    checkpoint_type checkpoint,
                    column_type column,
                    wfa_type& wfa)
                {
#ifdef WFA_TESTS
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = wfah_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::WfahCheckpointedScoringContext<BAND_LEN, TYPE, algorithm_tag, checkpoint_type> context(checkpoint);

                        return wfah_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(aligner.scheme, context, pattern, quals, text, min_score, sink, window_begin, window_end, column, wfa);
                    }
#else
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = sw_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::SWCheckpointedScoringContext<BAND_LEN, TYPE, algorithm_tag, checkpoint_type> context(checkpoint);

                        return sw_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(EditDistanceSWScheme(), context, pattern, quals, text, min_score, sink, window_begin, window_end, column, wfa);
                    }
#endif
                }

                /// dispatch scoring in a window of the pattern, retaining the intermediate results in the column
                /// vector, essentially used as a continuation
                ///
                /// \param aligner      scoring scheme
                /// \param pattern      pattern string (horizontal)
                /// \param quals        pattern qualities
                /// \param text         text string (vertical)
                /// \param min_score    minimum score
                /// \param sink         output alignment sink
                /// \param column       temporary column storage
                ///
                /// \return             true iff the minimum score was reached
                ///
                template <typename sink_type, typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool dispatch(
                    const aligner_type aligner,
                    const pattern_string pattern,
                    const qual_string quals,
                    const text_string text,
                    const int32 min_score,
                    const uint32 window_begin,
                    const uint32 window_end,
                    sink_type &sink,
                    column_type column,
                    wfa_type& wfa)
                {
#ifdef WFA_TESTS
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = wfah_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::WfahScoringContext<BAND_LEN, TYPE, algorithm_tag> context;

                        return wfah_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(aligner.scheme, context, pattern, quals, text, min_score, sink, window_begin, window_end, column, wfa);
                    }
#else
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = sw_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::SWScoringContext<BAND_LEN, TYPE, algorithm_tag> context;

                        return sw_alignment_score_dispatch<BAND_LEN, TYPE, algorithm_tag, symbol_type>::run(EditDistanceSWScheme(), context, pattern, quals, text, min_score, sink, window_begin, window_end, column, wfa);
                    }
#endif
                }
            };

            //
            // Calculate the alignment score between a pattern and a text, using the Wfah algorithm.
            //
            // \tparam CHECKPOINTS  number of columns between each checkpoint
            //
            // \param aligner      scoring scheme
            // \param pattern      pattern string (horizontal
            // \param quals        pattern qualities
            // \param text         text string (vertical)
            // \param min_score    minimum score
            //
            template <
                uint32 CHECKPOINTS,
                AlignmentType TYPE,
                typename scoring_type,
                typename pattern_string,
                typename qual_string,
                typename text_string,
                typename column_type>
            struct alignment_checkpointed_dispatch<
                CHECKPOINTS,
                WfahAligner<TYPE, scoring_type>,
                pattern_string,
                qual_string,
                text_string,
                column_type>
            {
                typedef WfahAligner<TYPE, scoring_type> aligner_type;

                //
                // Calculate a set of checkpoints of the DP matrix for the alignment between a pattern
                // and a text, using the Wfah-Smith-Waterman algorithm.
                //
                // \tparam checkpoint_type     a class to represent the collection of checkpoints,
                //                             represented as a linear array storing each checkpointed
                //                             band contiguously.
                //                             The class has to provide the const indexing operator[].
                //
                template <
                    typename sink_type,
                    typename checkpoint_type,
                    typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void dispatch_checkpoints(
                    const aligner_type aligner,
                    const pattern_string pattern,
                    const qual_string quals,
                    const text_string text,
                    const int32 min_score,
                    sink_type &sink,
                    checkpoint_type checkpoints,
                    column_type column,
                    wfa_type& wfa
                    )
                {
#ifdef WFA_TESTS
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = wfah_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::WfahCheckpointContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type> context(checkpoints);

                        wfah_alignment_score_dispatch<BAND_LEN, TYPE, PatternBlockingTag, symbol_type>::run(aligner.scheme, context, pattern, quals, text, min_score, sink, 0, pattern.length(), column, wfa);
                    }
#else
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = sw_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::SWCheckpointContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type> context(checkpoints);

                        sw_alignment_score_dispatch<BAND_LEN, TYPE, PatternBlockingTag, symbol_type>::run(EditDistanceSWScheme(), context, pattern, quals, text, min_score, sink, 0, pattern.length(), column, wfa);
                    }
#endif
                }

                //
                // Compute the banded Dynamic Programming submatrix between two given checkpoints,
                // storing its flow at each cell.
                // The function returns the submatrix width.
                //
                // \tparam BAND_LEN            size of the DP band
                //
                // \tparam checkpoint_type     a class to represent the collection of checkpoints,
                //                             represented as a linear array storing each checkpointed
                //                             band contiguously.
                //                             The class has to provide the const indexing operator[].
                //
                // \tparam submatrix_type      a class to store the flow H, E and F submatrix, represented
                //                             as a linear array of size (BAND_LEN*CHECKPOINTS).
                //                             The class has to provide the non-const indexing operator[].
                //                             Note that the H submatrix entries can assume only 3 values,
                //                             while the E and F only 2 - hence the aggregate needs 4 bits
                //                             per cell.
                //
                // \param checkpoints          the set of checkpointed rows
                // \param checkpoint_id        the starting checkpoint used as the beginning of the submatrix
                // \param submatrix            the output submatrix
                //
                // \return                     the submatrix width
                //
                template <
                    typename checkpoint_type,
                    typename submatrix_type,
                    typename wfa_type>
                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static uint32 dispatch_submatrix(
                    const aligner_type aligner,
                    const pattern_string pattern,
                    const qual_string quals,
                    const text_string text,
                    const int32 min_score,
                    checkpoint_type checkpoints,
                    const uint32 checkpoint_id,
                    submatrix_type submatrix,
                    column_type column,
                    wfa_type& wfa)
                {
#ifdef WFA_TESTS
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = wfah_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::WfahSubmatrixContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type, submatrix_type>
                            context(checkpoints, checkpoint_id, submatrix);

                        const uint32 window_begin = checkpoint_id * CHECKPOINTS;
                        const uint32 window_end = nvbio::min(window_begin + CHECKPOINTS, uint32(pattern.length()));

                        NullSink null_sink;
                        wfah_alignment_score_dispatch<BAND_LEN, TYPE, PatternBlockingTag, symbol_type>::run(aligner.scheme, context, pattern, quals, text, min_score, null_sink, window_begin, window_end, column, wfa);

                        return window_end - window_begin;
                    }
#else
                    {
                        typedef typename pattern_string::value_type symbol_type;

                        NVBIO_VAR_UNUSED const uint32 BAND_LEN = sw_bandlen_selector<TYPE, 1u, symbol_type>::BAND_LEN;

                        priv::SWSubmatrixContext<BAND_LEN, TYPE, CHECKPOINTS, checkpoint_type, submatrix_type>
                            context(checkpoints, checkpoint_id, submatrix);

                        const uint32 window_begin = checkpoint_id * CHECKPOINTS;
                        const uint32 window_end = nvbio::min(window_begin + CHECKPOINTS, uint32(pattern.length()));

                        NullSink null_sink;
                        sw_alignment_score_dispatch<BAND_LEN, TYPE, PatternBlockingTag, symbol_type>::run(EditDistanceSWScheme(), context, pattern, quals, text, min_score, null_sink, window_begin, window_end, column, wfa);

                        return window_end - window_begin;
                    }
#endif
                }
            };

            //
            // Given the Dynamic Programming submatrix between two checkpoints,
            // backtrace from a given destination cell, using Wfah's algorithm.
            // The function returns the resulting source cell.
            //
            // \tparam BAND_LEN            size of the DP band
            //
            // \tparam CHECKPOINTS         number of DP rows between each checkpoint
            //
            // \tparam checkpoint_type     a class to represent the collection of checkpoints,
            //                             represented as a linear array storing each checkpointed
            //                             band contiguously.
            //                             The class has to provide the const indexing operator[].
            //
            // \tparam submatrix_type      a class to store the flow submatrix, represented
            //                             as a linear array of size (BAND_LEN*CHECKPOINTS).
            //                             The class has to provide the const indexing operator[].
            //                             Note that the submatrix entries can assume only 3 values,
            //                             and could hence be packed in 2 bits.
            //
            // \tparam backtracer_type     a class to store the resulting list of backtracking operations.
            //                             A model of \ref Backtracer.
            //
            // \param checkpoints          precalculated checkpoints
            // \param checkpoint_id        index of the first checkpoint defining the DP submatrix,
            //                             storing all bands between checkpoint_id and checkpoint_id+1.
            // \param submatrix            precalculated flow submatrix
            // \param submatrix_height     submatrix width
            // \param submatrix_height     submatrix height
            // \param sink                 in/out sink of the DP solution
            // \param output               backtracking output handler
            //
            // \return                     true if the alignment source has been found, false otherwise
            //
            template <
                uint32 CHECKPOINTS,
                AlignmentType TYPE,
                typename scoring_type,
                typename checkpoint_type,
                typename submatrix_type,
                typename backtracer_type,
                typename wfa_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE bool alignment_traceback(
                const WfahAligner<TYPE, scoring_type> aligner,
                checkpoint_type checkpoints,
                const uint32 checkpoint_id,
                submatrix_type submatrix,
                const uint32 submatrix_width,
                const uint32 submatrix_height,
                uint8 &state,
                uint2 &sink,
                backtracer_type &output,
                wfa_type &wfa)
            {
#ifdef WFA_TESTS
                {
                    //
                    // Backtrack from the second checkpoint to the first looking up the flow submatrix.
                    //
                    int32 current_row = sink.x;
                    int32 current_col = sink.y - checkpoint_id * CHECKPOINTS; //- 1u;

                    /*NVBIO_CUDA_DEBUG_ASSERT(current_row > 0 &&
                                                current_row <= (int32)submatrix_height,
                                            "sw::alignment_backtrack(): sink (%u,%u) -> local x coordinate %d not in [0,%d[\n", sink.x, sink.y, current_row, submatrix_height);
                    NVBIO_CUDA_DEBUG_ASSERT(current_col >= 0, "sw::alignment_backtrack(): sink (%u,%u) -> local y coordinate %d not in [0,%u[ (checkpt %u)\n", sink.x, sink.y, current_col, submatrix_width, checkpoint_id);
                    NVBIO_CUDA_DEBUG_ASSERT(current_col < (int32)submatrix_width, "sw::alignment_backtrack(): sink (%u,%u) -> local y coordinate %d not in [0,%u[ (checkpt %u)\n", sink.x, sink.y, current_col, submatrix_width, checkpoint_id);
                    */
                    while (current_row > 0 &&
                           current_col >= 0)
                    {
                        const uint32 submatrix_cell = (current_row)*CHECKPOINTS + current_col;
                        const uint8 op = submatrix[submatrix_cell];

                        if (TYPE == LOCAL)
                        {
                            if (op == SINK)
                            {
                                sink.x = current_row;
                                sink.y = current_col + checkpoint_id * CHECKPOINTS; //+ 1u;
                                return true;
                            }
                        }

                        if (op != DELETION)
                            --current_col; // move to the previous column

                        if (op != INSERTION)
                            --current_row; // move to the previous row

                        output.push(op);
                    }

                    sink.x = current_row;

                    sink.y = current_col + checkpoint_id * CHECKPOINTS; //+ 1u;
                    return true;
                    // return current_row ? false : true; // if current_row == 0 we reached the end of the alignment along the text
                }
#else
                {
                    return alignment_traceback<CHECKPOINTS>(
                        make_smith_waterman_aligner<TYPE>(EditDistanceSWScheme()),
                        checkpoints,
                        checkpoint_id,
                        submatrix,
                        submatrix_width,
                        submatrix_height,
                        state,
                        sink,
                        output,
                        wfa);
                }
#endif
            }

        } // namespace priv

    } // namespace aln
} // namespace nvbio

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

#include <nvbio/basic/cuda/arch.h>
#include <nvbio/basic/cuda/scan.h>
#include <nvbio/alignment/sink.h>
#include <nvbio/alignment/utils.h>
#include <nvbio/alignment/warp_utils.h>
#include <nvbio/alignment/alignment_base_inl.h>


namespace nvbio {
namespace aln {
namespace priv {

template <
    uint32          BLOCKDIM,
    AlignmentType   TYPE,
    typename        scoring_type,
    typename        string_type,
    typename        qual_type,
    typename        ref_type,
    typename        column_type>
NVBIO_FORCEINLINE NVBIO_DEVICE
int32 wfah_alignment_score(
    const scoring_type&         scoring,
    const string_type           str,
    const qual_type             quals,
    const ref_type              ref,
    const int32                 min_score,
    uint2*                      sink,
    column_type                 temp)
{
#if __CUDA_ARCH__ >= 350
    typedef int32                        score_type;
    typedef alignment_result<score_type> alignment;


    const uint32 WARP_SIZE = 1u << cuda::Arch::LOG_WARP_SIZE;
    const uint32 NUM_WARPS = BLOCKDIM >> cuda::Arch::LOG_WARP_SIZE;

    const uint32 M = str.length();
    const uint32 N = ref.length();

    const score_type SCORE_TEXT_GAP_OPEN      = scoring.text_gap_open();
    const score_type SCORE_TEXT_GAP_EXTEND    = scoring.text_gap_extension();
    const score_type SCORE_PATTERN_GAP_OPEN   = scoring.pattern_gap_open();
    const score_type SCORE_PATTERN_GAP_EXTEND = scoring.pattern_gap_extension();

    // local scores
    score_type h_top, h_left, h_diag, hi;
    score_type e_top, ei;
    score_type f_left, fi;

    // local maximum score and corresponding sink
    alignment best_alignment = alignment::minimum_value();

    // current reference string character
    uint8 r_j;

    // per-thread cache for temp values and reference string characters
    // each thread loads a different value; cache values are shuffled down the warp at each iteration
    score_type temp_cache_h, temp_cache_f;
    uint8 reference_cache;

    // width of the current warp-block stripe of the DP matrix (always WARP_SIZE except for the last stripe)
    uint32 warp_block_width;

    // compute warp-block horizontal coordinate in DP matrix for this thread
    const uint32 wi = warp_tid() + 1;

    // initialize the leftmost matrix column
    for(uint32 i = warp_tid(); i < N; i += WARP_SIZE)
    {
        temp[i] = make_short2(
            (TYPE == GLOBAL ? SCORE_TEXT_GAP_OPEN + SCORE_TEXT_GAP_EXTEND * i : 0),
            (TYPE != LOCAL ? Field_traits<int16>::min() : 0) );
    }

    for(uint32 warp_block = 0; warp_block < M; warp_block += WARP_SIZE)
    {
        // width of this block
        warp_block_width = (warp_block + WARP_SIZE >= M ? M % WARP_SIZE : WARP_SIZE);
        // compute the horizontal coordinate of the current thread in the DP matrix (including border column)
        const uint32 i = wi + warp_block;

        // set top boundary values
        h_top = (TYPE != LOCAL ? SCORE_PATTERN_GAP_OPEN + SCORE_PATTERN_GAP_EXTEND * (i - 1) : 0);
        e_top = (TYPE != LOCAL ? Field_traits<int16>::min() : 0);

        // initialize diagonal
        h_diag = (TYPE != LOCAL ? SCORE_PATTERN_GAP_OPEN + SCORE_PATTERN_GAP_EXTEND * (i - 2) : 0);

        // load the query string character for the current thread
        const uint8 s_i = (i <= M ? str[i - 1] : 0);
        const uint8 q_i = (i <= M ? quals[i - 1] : 0);

        // initialize the best score for this stripe
        score_type max_score = Field_traits<score_type>::min();

        // loop over all DP anti-diagonals, excluding the border row/column
        for(uint32 block_diag = 2; block_diag <= warp_block_width + N; block_diag += WARP_SIZE)
        {
            // reload caches every WARP_SIZE diagonals
            const uint32 thread_j = (block_diag - 2) + warp_tid();

            if (thread_j < N)
            {
                temp_cache_h = temp[thread_j].x;
                temp_cache_f = temp[thread_j].y;
                reference_cache = ref[(block_diag - 2) + warp_tid()];
            } else {
                temp_cache_h = 0;
                temp_cache_f = 0;
                reference_cache = 0;
            }

            for(uint32 diag = block_diag; diag < block_diag + WARP_SIZE; diag++)
            {
                // compute the length of this anti-diagonal (excluding border row/column)
                const uint32 diag_len = nvbio::min3(diag - 1, WARP_SIZE, warp_block_width);
                // compute vertical coordinate of the current cell in the DP matrix (including border column)
                const uint32 j = diag - wi;

                // is the current cell inside the DP matrix?
                if (wi <= diag_len && j <= N)
                {
                    if (wi == 1)
                    {
                        // load new temp and reference values
                        r_j = reference_cache;
                        // initialize cell to the left of the current cell
                        h_left = temp_cache_h;
                        f_left = temp_cache_f;
                    }

                    // compute the match/mismatch score
                    const score_type S_ij = (r_j == s_i) ? scoring.match(q_i) : scoring.mismatch(q_i);

                    ei = nvbio::max(e_top + SCORE_TEXT_GAP_EXTEND,
                                    h_top + SCORE_TEXT_GAP_OPEN);
                    fi = nvbio::max(f_left + SCORE_PATTERN_GAP_EXTEND,
                                    h_left + SCORE_PATTERN_GAP_OPEN);
                    hi = nvbio::max3(h_diag + S_ij,
                                     ei,
                                     fi);

                    if (TYPE == LOCAL)
                    {
                        // clamp score to zero
                        hi = nvbio::max(hi, score_type(0));
                    }

                    // save off the last column
                    if (wi == WARP_SIZE)
                    {
                        temp[j - 1] = make_short2( hi, fi );

                        // keep track of the best score in this stripe
                        max_score = nvbio::max( max_score, hi );
                    }

                    // save the best score across the entire matrix for local scoring
                    // save the best score across the last column for semi-global scoring
                    if (TYPE == LOCAL ||
                        (TYPE == SEMI_GLOBAL && i == M))
                    {
                        if (hi > best_alignment.score)
                            best_alignment = alignment(hi, make_uint2(j, i));
                    }

                    // current left becomes diagonal for next iteration on this lane
                    h_diag = h_left;

                    // current value becomes h_top for next iteration on this lane
                    h_top = hi;
                    e_top = ei;
                }

                // move previous cell reference value across the warp
                r_j = __shfl_up(r_j, 1);
                // hi becomes h_left on the next lane
                h_left = __shfl_up(hi, 1);
                f_left = __shfl_up(fi, 1);

                // push temp_cache and reference_cache values down the warp
                temp_cache_h = __shfl_down(temp_cache_h, 1);
                temp_cache_f = __shfl_down(temp_cache_f, 1);
                reference_cache = __shfl_down(reference_cache, 1);
            }
        }

        if (warp_block + WARP_SIZE < M)
        {
            // we are now (M - warp_block - WARP_SIZE) columns away from the last one: check whether
            // we could theoretically reach the minimum score
            max_score = __shfl( max_score, WARP_SIZE - 1 );

            const score_type missing_cols = score_type(M - warp_block - WARP_SIZE);
            if (max_score + missing_cols * scoring.match(255) < score_type( min_score ))
                return Field_traits<int32>::min();
        }
    }

    if (TYPE == LOCAL || TYPE == SEMI_GLOBAL)
    {
        // do a warp-wide max-scan to find the largest score (TODO: use a reduction instead)
        __shared__ volatile alignment sm_red [WARP_SIZE * NUM_WARPS * 2];
        volatile alignment *sm_warp_red = sm_red + WARP_SIZE * warp_id() * 2;
        cuda::scan<32>(best_alignment, alignment::max_operator(), alignment::minimum_value(), sm_warp_red);
        best_alignment = cuda::scan_total<32>(sm_warp_red);
    }

    if (TYPE == GLOBAL)
    {
        best_alignment.score = __shfl(hi, warp_block_width - 1);
        best_alignment.sink = make_uint2(N,M);
    }

    *sink = best_alignment.sink;
    return best_alignment.score;
#else
    // unsupported on compute capability < 3.5
    return 0;
#endif
}

// private dispatcher for the warp-parallel version of the Wfah aligner
template <
    uint32          BLOCKDIM,
    AlignmentType   TYPE,
    typename        scoring_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        column_type>
NVBIO_FORCEINLINE NVBIO_DEVICE
int32 alignment_score(
    const WfahAligner<TYPE,scoring_type>       aligner,
    const pattern_string                        pattern,
    const qual_string                           quals,
    const text_string                           text,
    const  int32                                min_score,
          uint2*                                sink,
          column_type                           column)
{
#if defined(NVBIO_DEVICE_COMPILATION)
    return wfah_alignment_score<BLOCKDIM,TYPE>(
        aligner.scheme,
        pattern,
        quals,
        text,
        min_score,
        sink,
        column );
#else
    return Field_traits<int32>::min();
#endif
}

} // namespace priv
} // namespace aln
} // namespace nvbio

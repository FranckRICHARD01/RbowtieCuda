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

//
// Myers bit-vector scoring for the semi-global edit-distance path used by the
// opposite-mate rescue stage.
//
// Semantics replicated bit-for-bit from the legacy full-DP path
// (EditDistanceAligner<SEMI_GLOBAL>, PatternBlockingTag, BAND_LEN = 16):
//   - free start at every text column, sinks at pattern column M for each text row;
//   - sink reported through BestSink::report(score, uint2(i+1, M)), i.e. rightmost
//     tie-break on equality;
//   - stripe pruning: after each stripe of 16 pattern columns except the last one
//     (end_block = max(16, 16*ceil(M/16))), if the maximum over all text rows of the
//     boundary column value is below min_score the whole alignment fails and nothing
//     is reported (sinks are only produced by the final stripe).
//
// Validated exhaustively against the legacy path in /tmp/kilo/spike_opposite.
// Falls back to the exact legacy dispatch whenever the fast path cannot apply
// (pattern longer than 256 symbols, or any symbol outside the 0..3 DNA alphabet).
//

#pragma once

#include <nvbio/basic/types.h>
#include <nvbio/basic/numbers.h>
#include <nvbio/alignment/alignment_base.h>
#include <nvbio/alignment/utils.h>
#include <nvbio/alignment/utils_inl.h>
#include <nvbio/alignment/sw/sw_inl.h>
#include <nvbio/alignment/ed/ed_inl.h>

namespace nvbio {
namespace aln {

///
/// An edit-distance aligner scored with the Myers bit-vector algorithm, replicating
/// exactly the semi-global sink/pruning semantics of EditDistanceAligner<SEMI_GLOBAL>.
///
template <AlignmentType T_TYPE = SEMI_GLOBAL, typename AlgorithmType = PatternBlockingTag>
struct MyersOppositeAligner
{
    static const AlignmentType TYPE =   T_TYPE;         ///< the AlignmentType
    typedef EditDistanceTag             aligner_tag;    ///< scored as plain edit distance
    typedef AlgorithmType               algorithm_tag;  ///< the \ref AlgorithmTag "Algorithm Tag"
};

/// column storage: keep the very same footprint as the edit-distance path so that
/// batched buffer sizing is unaffected
///
template <AlignmentType TYPE, typename algorithm_tag>
struct column_storage_type< MyersOppositeAligner<TYPE,algorithm_tag> > { typedef int16 type; };

//
/// maximum number of reference gaps within a score boundary: identical to edit distance
///
template <AlignmentType TYPE, typename algorithm_tag>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
uint32 max_text_gaps(
    const MyersOppositeAligner<TYPE,algorithm_tag>&   aligner,
    int32                                             min_score,
    int32                                             pattern_len)
{
    return -min_score;
}

namespace priv {
namespace myers {

//
// Myers multi-word bit-vector core for the semi-global ED sink-anywhere problem.
// Reports exactly one sink (the BestSink max with <= tie-break is reproduced by a
// single report of the argmax over all text rows) and emulates the legacy stripe
// pruning; when pruned, nothing is reported and false is returned.
//
// \tparam W     number of 64-bit words covering the pattern (1..4)
//
template <uint32 W, typename pattern_string, typename text_string, typename sink_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool score(
    const pattern_string   pattern,
    const uint32           M,
    const text_string      text,
    const uint32           N,
    const int32            min_score,
          sink_type&       sink)
{
    uint64 Pv[W], Mv[W], Peq[4][W];
    #pragma unroll
    for (uint32 w = 0; w < W; ++w) { Pv[w] = ~uint64(0); Mv[w] = 0; }
    #pragma unroll
    for (uint32 s = 0; s < 4; ++s)
        #pragma unroll
        for (uint32 w = 0; w < W; ++w) Peq[s][w] = 0;

    for (uint32 j = 0; j < M; ++j) {
        const uint8 s = pattern[j];
        if (s < 4u) Peq[s][j >> 6] |= uint64(1) << (j & 63);
    }

    // legacy stripe pruning checkpoints: columns 16, 32, ... strictly before end_block,
    // with end_block = max(16, 16*ceil(M/16))
    const uint32 end_block = nvbio::max(uint32(16), uint32(16u * ((M + 15u) / 16u)));
    const uint32 ncp       = (end_block - 16u) / 16u;

    int32 cpmx[16];
    #pragma unroll
    for (uint32 k = 0; k < 16; ++k) cpmx[k] = Field_traits<int32>::min();

    int32 Dm = int32(M);
    bool  have = false;
    int32 best = 0; uint32 sinkx = 0;

    for (uint32 i = 0; i < N; ++i) {
        uint64 carry = 0, phc = 0, mhc = 0;
        const uint8 c = text[i];
        #pragma unroll
        for (uint32 w = 0; w < W; ++w) {
            const uint64 Eq   = (c < 4u) ? Peq[c][w] : 0;
            const uint64 Xv   = Eq | Mv[w];
            const uint64 a    = Eq & Pv[w];
            const uint64 s1   = a + Pv[w];
            const uint64 c1   = (s1 < a) ? 1u : 0u;
            const uint64 s2   = s1 + carry;
            const uint64 c2   = (s2 < s1) ? 1u : 0u;
            const uint64 Xh   = (s2 ^ Pv[w]) | Eq;
                  uint64 Ph   = Mv[w] | ~(Xh | Pv[w]);
                  uint64 Mh   = Pv[w] & Xh;
            if (w == W - 1) {
                const uint64 hb = uint64(1) << ((M - 1u) & 63u);
                if (Ph & hb) ++Dm;
                if (Mh & hb) --Dm;
            }
            const uint64 nphc = Ph >> 63;
            const uint64 nmhc = Mh >> 63;
                  Ph = (Ph << 1) | phc;
                  Mh = (Mh << 1) | mhc;
                  Pv[w] = Mh | ~(Xv | Ph);
                  Mv[w] = Ph & Xv;
                  carry = c1 | c2;
                  phc   = nphc;
                  mhc   = nmhc;
        }

        // sink at pattern column M for this text row (rightmost tie-break kept by argmax)
        const int32 sc = -Dm;
        if (!have || sc >= best) { have = true; best = sc; sinkx = i + 1; }

        // boundary-column values at each checkpoint, via vertical-delta popcounts:
        // D(pp) = sum over bits < pp of (+1 for Pv bit, -1 for Mv bit), score = -D(pp)
        if (ncp) {
            int32 cum[W+1];
            cum[0] = 0;
            #pragma unroll
            for (uint32 w = 0; w < W; ++w)
                cum[w+1] = cum[w] + int32(popc(Pv[w])) - int32(popc(Mv[w]));
            #pragma unroll
            for (uint32 k = 0; k < ncp; ++k) {
                const uint32 pp = 16u * (k + 1u);
                const uint32 wq = pp >> 6;
                const uint32 b  = pp & 63u;
                const uint64 mask = (uint64(1) << b) - 1;
                const int32 partial = int32(popc(Pv[wq] & mask)) - int32(popc(Mv[wq] & mask));
                const int32 Dp = cum[wq] + partial;
                if (-Dp > cpmx[k]) cpmx[k] = -Dp;
            }
        }
    }

    #pragma unroll
    for (uint32 k = 0; k < 16; ++k) {
        if (k < ncp && cpmx[k] < min_score)
            return false;   // legacy prunes before the final stripe: nothing reported
    }

    if (have) sink.report(best, make_uint2(sinkx, M));
    return true;
}

/// select the word count at runtime; only called when 1 <= M <= 256
///
template <typename pattern_string, typename text_string, typename sink_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool score_dispatch(
    const pattern_string   pattern,
    const uint32           M,
    const text_string      text,
    const uint32           N,
    const int32            min_score,
          sink_type&       sink)
{
    if (M <= 64u)  return score<1>(pattern, M, text, N, min_score, sink);
    if (M <= 128u) return score<2>(pattern, M, text, N, min_score, sink);
    if (M <= 192u) return score<3>(pattern, M, text, N, min_score, sink);
    return score<4>(pattern, M, text, N, min_score, sink);
}

} // namespace myers

///
/// dispatch scoring across the whole pattern with the Myers bit-vector core,
/// falling back to the exact legacy edit-distance path when unsupported
///
template <
    AlignmentType TYPE,
    typename      algorithm_tag,
    typename      pattern_string,
    typename      qual_string,
    typename      text_string,
    typename      column_type>
struct alignment_score_dispatch<
    MyersOppositeAligner<TYPE,algorithm_tag>,
    pattern_string,
    qual_string,
    text_string,
    column_type>
{
    typedef MyersOppositeAligner<TYPE,algorithm_tag> aligner_type;

    template <typename sink_type, typename wfa_type>
    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
    static bool dispatch(
        const aligner_type      aligner,
        const pattern_string    pattern,
        const qual_string       quals,
        const text_string       text,
        const  int32            min_score,
              sink_type&        sink,
              column_type       column,
              wfa_type&         wfa)
    {
        typedef EditDistanceAligner<TYPE,algorithm_tag> ed_aligner_type;

        const uint32 M = pattern.length();
        const uint32 N = text.length();

        // fast path only covers the DNA alphabet and patterns up to 4 words
        bool supported = (M >= 1u) && (M <= 256u);
        if (supported) {
            for (uint32 j = 0; j < M && supported; ++j)
                if (pattern[j] >= 4u) supported = false;
            for (uint32 i = 0; i < N && supported; ++i)
                if (text[i] >= 4u) supported = false;
        }

        if (!supported) {
            // exact legacy behavior, byte-for-byte identical to the previous opposite path
            return alignment_score_dispatch<ed_aligner_type,pattern_string,qual_string,text_string,column_type>::dispatch(
                ed_aligner_type(), pattern, quals, text, min_score, sink, column, wfa);
        }

        // supported: M in [1,256] and DNA-only alphabet; the core always runs and its
        // return value mirrors the legacy success/pruning status exactly
        return myers::score_dispatch(pattern, M, text, N, min_score, sink);
    }
};

} // namespace priv

} // namespace aln
} // namespace nvbio

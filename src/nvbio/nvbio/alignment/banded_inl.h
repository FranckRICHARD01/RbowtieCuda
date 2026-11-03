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
#include <nvbio/alignment/utils.h>
#include <nvbio/alignment/alignment_base_inl.h>
#include <nvbio/alignment/ed/ed_banded_inl.h>
#include <nvbio/alignment/sw/sw_banded_inl.h>
#include <nvbio/alignment/gotoh/gotoh_banded_inl.h>
#include <nvbio/alignment/wfah/wfah_banded_inl.h>
#include <nvbio/alignment/myers/myers_banded_inl.h>
namespace nvbio {
namespace aln {

#if !defined(__CUDA_ARCH__)
//
// A wavefront aligner answered on the host is answered by the affine banded DP carrying the very
// penalties the aligner declares -- the same cost model the wavefront charges, a gap of L symbols
// costing gap_open + L * gap_extension, and one running under no wavefront budget.
//
// The wavefront engines work by credits : a cell is credited with the wavefront it was reached on,
// and where the rows run out the traceback completes a tail up to that credit -- a tail the
// wavefront never proved.  Such a report costs, base by base, more than the number it advertises,
// and no alignment of the window costs the advertised number at all : the score has no witness.
// The banded DP of the same scheme traces a real path and costs exactly what that path costs, so
// on the host it is taken as the authority of the window, and score and alignment become one and
// the same object.  Device code is left to the wavefront engines, which carry the per-window
// budget the host does not have.
//
//
// A bounded local witness for the traceback: the emitted operations are written here, and the
// number reported is computed from them.  A path that overflows the record is left incomplete,
// which makes the computation fall back on the DP score rather than corrupt memory.
//
struct OpsBacktracker
{
    enum { CAPACITY = 4096u };

    OpsBacktracker() : aln_len(0u) { aln[0] = '\0'; }

    void clear() { aln_len = 0u; aln[0] = '\0'; }

    void push(const uint8 op)
    {
        if (op < 3u && aln_len + 1u < CAPACITY)
        {
            aln[aln_len++] = "MID"[op];
            aln[aln_len] = '\0';
        }
    }
    void clip(const uint32 len)
    {
        for (uint32 i = 0; i < len; ++i)
        {
            if (aln_len + 1u >= CAPACITY) break;
            aln[aln_len++] = 'S';
        }
        aln[aln_len] = '\0';
    }

    char   aln[CAPACITY];
    uint32 aln_len;
};

template <
    uint32          BAND_LEN,
    AlignmentType   TYPE,
    typename        scoring_type,
    typename        algorithm_tag,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        sink_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool banded_alignment_score_authority(
    WfahTag,
    const WfahAligner<TYPE, scoring_type, algorithm_tag> &aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
          sink_type&        sink,
          wfa_type&         wfa)
{
    const bool answered = banded_alignment_score<BAND_LEN>(
        GotohAligner<TYPE, WfahScoreScheme<scoring_type> >(
            WfahScoreScheme<scoring_type>(aligner.scheme)),
        pattern,
        quals,
        text,
        min_score,
        sink,
        wfa);

    // The number advertised must be the cost of the path the traceback is about to emit, not the
    // value of a cell: same aligner, same DP, a full traceback into a local witness, and the cost
    // of the operations it writes.
    if (answered && (uint32) pattern.size() <= 1024u)
    {
        OpsBacktracker ops;
        const Alignment<int32> traced = banded_alignment_traceback<BAND_LEN, 1024u, 1000u>(
            GotohAligner<TYPE, WfahScoreScheme<scoring_type> >(
                WfahScoreScheme<scoring_type>(aligner.scheme)),
            pattern, quals, text, min_score, ops, wfa);

        const int32 cost = backtracker_ops_cost_impl(ops, (int32) traced.source.x, pattern, text, aligner.scheme, 0);
        if (cost != Field_traits<int32>::min())
            sink.score = -cost;
    }
    return answered;
}

// every other aligner answers for itself, being DP based already
//
template <
    uint32          BAND_LEN,
    typename        aligner_tag,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        sink_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool banded_alignment_score_authority(
    aligner_tag,
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
          sink_type&        sink,
          wfa_type&         wfa)
{
    return priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        sink,
        wfa);
}
#endif

// a function to compute the alignment score between a pattern and a text string 
// with banded DP alignment.
//
// \param aligner             alignment algorithm
// \param pattern             pattern string
// \param quals               quality string
// \param text                text string
// \param min_score           threshold alignment score
// \param sink                output sink
//
template <
    uint32          BAND_LEN,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        sink_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool banded_alignment_score(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const  int32            min_score,
          sink_type&        sink,
          wfa_type&         wfa)
{
#if !defined(__CUDA_ARCH__)
    // on the host, a wavefront aligner is answered by the affine banded DP of its own scheme
    return banded_alignment_score_authority<BAND_LEN>(
        typename aligner_type::aligner_tag(),
        aligner,
        pattern,
        quals,
        text,
        min_score,
        sink,
        wfa );
#else
    return priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        sink,
        wfa );
#endif
}


// a function to compute the alignment score between a pattern and a text string 
// with banded DP alignment.
//
// \param aligner             alignment algorithm
// \param pattern             pattern string
// \param quals               quality string
// \param text                text string
// \param min_score           threshold alignment score
// \param sink                output sink
//
template <
    uint32          BAND_LEN,
    typename        aligner_type,
    typename        pattern_string,
    typename        text_string,
    typename        sink_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool banded_alignment_score(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const text_string       text,
    const  int32            min_score,
          sink_type&        sink,
          wfa_type&         wfa)
{
    return priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        trivial_quality_string(),
        text,
        min_score,
        sink,
        wfa);
}



// a function to compute the alignment score between a pattern and a text string 
// with banded DP alignment.
//
// \param aligner              alignment algorithm
// \param pattern              pattern string
// \param quals                quality string
// \param text                 text string
// \param min_score            threshold alignment score
//
// \return                      best alignment score
//
template <
    uint32          BAND_LEN,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
int32 banded_alignment_score(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const  int32            min_score,
    wfa_type&               wfa)
{
    BestSink<int32> sink; 
    priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        sink,
        wfa );
    return sink.score;
}


// a function to compute the alignment score between a pattern and a text string 
// with banded DP alignment.
//
// \param aligner              alignment algorithm
// \param pattern              pattern string
// \param text                 text string
// \param min_score            threshold alignment score
//
// \return                      best alignment score
//
template <
    uint32          BAND_LEN,
    typename        aligner_type,
    typename        pattern_string,
    typename        text_string,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
int32 banded_alignment_score(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const text_string       text,
    const int32             min_score,
    wfa_type&               wfa)
{
    BestSink<int32> sink;       
    priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        trivial_quality_string(),
        text,
        min_score,
        sink,
        wfa );
    return sink.score;
}

// a function to compute the alignment score between a pattern and a text string 
// with banded DP alignment in multiple passes, where each pass handles a window
// of the pattern and saves a checkpoint for the next.
//
// \param aligner             alignment algorithm
// \param pattern             pattern string
// \param quals               quality string
// \param text                text string
// \param min_score           threshold alignment score
// \param sink                output sink
//
template <
    uint32          BAND_LEN,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        sink_type,
    typename        checkpoint_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
bool banded_alignment_score(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const  int32            min_score,
    const uint32            window_begin,
    const uint32            window_end,
          sink_type&        sink,
          checkpoint_type   checkpoint,
          wfa_type&         wfa)
{
    return priv::banded_alignment_score<BAND_LEN>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        window_begin,
        window_end,
        sink,
        checkpoint,
        wfa );
}

//
// Calculate a set of checkpoints of the DP matrix for the alignment between a pattern
// and a text.
//
// \param aligner              scoring scheme
// \param pattern              pattern string (horizontal
// \param quals                pattern qualities
// \param text                 text string (vertical)
// \param min_score            minimum score
// \param sink                 output sink
// \param checkpoints          the set of output checkpointed rows
//
template <
    uint32          BAND_LEN,
    uint32          CHECKPOINTS,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        sink_type,
    typename        checkpoint_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
void banded_alignment_score_checkpoints(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const  int32            min_score,
          sink_type&        sink,
          checkpoint_type   checkpoints,
          wfa_type&         wfa)
{
    priv::banded_alignment_checkpoints<BAND_LEN,CHECKPOINTS>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        sink,
        checkpoints,
        wfa );
}

//
// Compute the banded Dynamic Programming submatrix between two given checkpoints,
// storing its flow at each cell.
// The function returns the submatrix width.
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
//                             The class has to provide the non-const indexing operator[].
//                             Note that the submatrix entries can assume only 3 values,
//                             and could hence be packed in 2 bits.
//
// \param aligner              scoring scheme
// \param pattern              pattern string (horizontal
// \param quals                pattern qualities
// \param text                 text string (vertical)
// \param min_score            minimum score
// \param checkpoints          the set of checkpointed rows
// \param checkpoint_id        the starting checkpoint used as the beginning of the submatrix
// \param submatrix            the output submatrix
//
// \return                     the submatrix width
//
template <
    uint32          BAND_LEN,
    uint32          CHECKPOINTS,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        checkpoint_type,
    typename        submatrix_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
uint32 banded_alignment_score_submatrix(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
    checkpoint_type         checkpoints,
    const uint32            checkpoint_id,
    submatrix_type          submatrix,
    wfa_type&               wfa)
{
    return priv::banded_alignment_submatrix<BAND_LEN,CHECKPOINTS>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        checkpoints,
        checkpoint_id,
        submatrix,
        wfa );
}


//
// Backtrace an optimal alignment using a full DP algorithm.
//
// \tparam CHECKPOINTS         number of DP rows between each checkpoint
// \tparam aligner_type        an \ref Aligner "Aligner" algorithm
// \tparam pattern_string      a string representing the pattern.
// \tparam qual_string         an array representing the pattern qualities.
// \tparam text_string         a string representing the text.
// \tparam checkpoints_type    an array-like class defining operator[], used to represent a reduced DP score matrix,
//                             containing all the matrix columns whose index is a multiple of CHECKPOINTS;
//                             the type of the matrix cells depends on the aligner, and can be obtained as
//                             typename checkpoints_storage_type<aligner_type>::type
// \tparam submatrix_type      an array-like class defining operator[], used to represent a temporary DP flow submatrix,
//                             containing all the matrix flow cells between two checkpointed columns.
//                             the number of bits needed for the submatrix cells depends on the aligner, and can be obtained as
//                             direction_vector_traits<aligner_type>::BITS
// \tparam backtracer_type     a model of \ref Backtracer.
//
// \param aligner              alignment algorithm
// \param pattern              pattern to be aligned
// \param quals                pattern quality scores
// \param text                 text to align the pattern to
// \param min_score            minimum accepted score
// \param backtracer           backtracking delegate
// \param checkpoints          temporary checkpoints storage
// \param submatrix            temporary submatrix storage
//
// \return                     reported alignment
// 
template <
    uint32          BAND_LEN,
    uint32          CHECKPOINTS,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        backtracer_type,
    typename        checkpoints_type,
    typename        submatrix_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
Alignment<int32> banded_alignment_traceback(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
    backtracer_type&        backtracer,
    checkpoints_type        checkpoints,
    submatrix_type          submatrix,
    wfa_type&               wfa)
{
    //
    // This function performs backtracing in a completely generic fashion that works
    // for any DP based aligner.
    // It does so using 2 module-wide generic functions:
    //      * alignment_score_checkpoints()
    //      * alignment_score_submatrix()
    //
    // and one aligner-specific function:
    //      * priv::alignment_traceback()
    //
    BestSink<int32> best;
    banded_alignment_score_checkpoints<BAND_LEN,CHECKPOINTS>(
        aligner, pattern, quals, text, /*pos,*/ min_score, best, checkpoints, wfa );

    const uint32 n_checkpoints = 1u + (pattern.length() + CHECKPOINTS-1)/CHECKPOINTS;

    uint2 sink = best.sink;
    if (sink.x == uint32(-1) ||
        sink.y == uint32(-1))
        return Alignment<int32>( best.score, make_uint2( uint32(-1), uint32(-1) ), make_uint2( uint32(-1), uint32(-1) ) );

    // clip the end of the alignment
    backtracer.clip( pattern.length() - sink.y );

    // find the checkpoint containing the sink
    int32 checkpoint_id = n_checkpoints-2;

    if (aligner_type::TYPE == LOCAL)
    {
        for (; checkpoint_id >= 0; --checkpoint_id)
        {
            if (checkpoint_id * CHECKPOINTS < sink.y)
                break;
        }
    }

    //store state (H, E, or F) between checkpoints
    uint8 state = HSTATE;

    // backtrack until needed
    for (; checkpoint_id >= 0; --checkpoint_id)
    {
        const uint32 submatrix_height = banded_alignment_score_submatrix<BAND_LEN,CHECKPOINTS>(
            aligner, pattern, quals, text, /*pos,*/ min_score, checkpoints, checkpoint_id, submatrix, wfa );

        if (priv::banded_alignment_traceback<BAND_LEN,CHECKPOINTS>(
            aligner, checkpoints, checkpoint_id, submatrix, submatrix_height, state, sink, backtracer, wfa ))
            break;
    }

    // clip the beginning of the alignment
    backtracer.clip( sink.y );

    return Alignment<int32>( best.score, sink, best.sink );
}

//
// Backtrace an optimal alignment using a full DP algorithm.
//
// \tparam MAX_PATTERN_LEN     maximum pattern length
// \tparam MAX_TEXT_LEN        maximum text length
// \tparam CHECKPOINTS         number of DP rows between each checkpoint
// \tparam aligner_type        an \ref Aligner "Aligner" algorithm
// \tparam pattern_string      a string representing the pattern.
// \tparam qual_string         an array representing the pattern qualities.
// \tparam text_string         a string representing the text.
// \tparam backtracer_type     a model of \ref Backtracer.
//
// \param aligner              alignment algorithm
// \param pattern              pattern to be aligned
// \param quals                pattern quality scores
// \param text                 text to align the pattern to
// \param min_score            minimum accepted score
// \param backtracer           backtracking delegate
//
// \return                     reported alignment
//
//
//
// The number the authority reports is the cost of the path it emits, not the value of a cell
// whose path diverges from it.  The operations are read back where the caller's witness keeps
// them (a table of M/I/D/S characters and its length), under the witness's own costing rule:
// the first symbol of a run costs gap_open + gap_extension, later symbols gap_extension, a
// pairing costs match or mismatch, a clip advances the pattern for free.  When the callee
// exposes no such members, the DP score is reported unchanged and the path is left untouched.
//
template <
    typename backtracer_type,
    typename scoring_type,
    typename pattern_string,
    typename text_string>
auto backtracker_ops_cost_impl(
    const backtracer_type&  bt,
    const int32             offset,
    const pattern_string    pattern,
    const text_string       text,
    const scoring_type&     scoring,
    int                     dispatch) -> decltype(bt.aln, bt.aln_len, int32())
{
    const int32  p_open  = scoring.pattern_gap_open();
    const int32  p_ext   = scoring.pattern_gap_extension();
    const int32  t_open  = scoring.text_gap_open();
    const int32  t_ext   = scoring.text_gap_extension();
    const int32  mat     = scoring.match();
    const int32  mm      = scoring.mismatch();
    const uint32 M       = (uint32) pattern.size();
    const uint32 N       = (uint32) text.size();
    (void) dispatch;
    const uint32 cap     = (uint32) (sizeof(bt.aln) / sizeof(bt.aln[0]));
    const uint32 n       = (bt.aln_len < cap) ? bt.aln_len : cap;

    // an empty path witnesses nothing: report the DP score as it stands
    if (n == 0u)
        return Field_traits<int32>::min();

    int32  score = 0;
    int32  state = -1;
    uint32 j     = 0u;
    uint32 k     = (uint32) offset;

    // the operations are pushed sink to source: the path is read backwards, as the witness reads it
    for (uint32 a = 0; a < n; ++a)
    {
        const char op = bt.aln[n - a - 1u];
        if (op == '\0')
            break;

        if (op == 'S')
        {
            if (j >= M) return Field_traits<int32>::min();
            ++j;
        }
        else if (op == 'I')
        {
            if (j >= M) return Field_traits<int32>::min();
            score += (state != 'I') ? p_open + p_ext : p_ext;
            ++j;
        }
        else if (op == 'D')
        {
            if (k >= N) return Field_traits<int32>::min();
            score += (state != 'D') ? t_open + t_ext : t_ext;
            ++k;
        }
        else if (op == 'M')
        {
            if (j >= M || k >= N) return Field_traits<int32>::min();
            score += (pattern[j] == text[k]) ? mat : mm;
            ++j; ++k;
        }
        else
            return Field_traits<int32>::min();

        state = (int32) op;
    }
    return score;
}

template <
    typename backtracer_type,
    typename scoring_type,
    typename pattern_string,
    typename text_string>
int32 backtracker_ops_cost_impl(
    const backtracer_type&,
    const int32,
    const pattern_string&,
    const text_string&,
    const scoring_type&,
    long                    dispatch)
{
    (void) dispatch;
    return Field_traits<int32>::min();
}

// A wavefront that declines a window reports nothing at all : the generic walk hands back a
// sink of (-1,-1), and the caller is left with an empty score and an empty path. On the host,
// hand the window to the affine banded DP carrying the very penalties the aligner declares --
// the same cost model the wavefront charges, gap_open + L * gap_extension, and one running
// under no wavefront bound, tracing a real path -- so that every reported score comes with the
// alignment proving it.
//
template <
    uint32          BAND_LEN,
    uint32          MAX_PATTERN_LEN,
    uint32          CHECKPOINTS,
    AlignmentType   TYPE,
    typename        scoring_type,
    typename        algorithm_tag,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        backtracer_type,
    typename        wfa_type,
    typename        checkpoints_type,
    typename        submatrix_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
Alignment<int32> banded_alignment_traceback_ondp(
    WfahTag,
    const WfahAligner<TYPE, scoring_type, algorithm_tag> &aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
    backtracer_type&        backtracer,
    wfa_type&               wfa,
    checkpoints_type,
    submatrix_type)
{
    backtracer.clear();
    const Alignment<int32> aln = banded_alignment_traceback<BAND_LEN, MAX_PATTERN_LEN, CHECKPOINTS>(
        GotohAligner<TYPE, WfahScoreScheme<scoring_type> >(
            WfahScoreScheme<scoring_type>(aligner.scheme)),
        pattern, quals, text, min_score, backtracer, wfa);

    // the number reported becomes the cost of the path as it is emitted, as the witness walks it
    Alignment<int32> reported = aln;
    const int32 cost = backtracker_ops_cost_impl(backtracer, (int32) aln.source.x, pattern, text, aligner.scheme, 0);
    if (cost != Field_traits<int32>::min())
        reported.score = -cost;

    return reported;
}

// every other aligner is DP based already : there is nothing under it to fall back to
//
template <
    uint32          BAND_LEN,
    uint32          MAX_PATTERN_LEN,
    uint32          CHECKPOINTS,
    typename        aligner_tag,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        backtracer_type,
    typename        wfa_type,
    typename        checkpoints_type,
    typename        submatrix_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
Alignment<int32> banded_alignment_traceback_ondp(
    aligner_tag,
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
    backtracer_type&        backtracer,
    wfa_type&               wfa,
    checkpoints_type        checkpoints,
    submatrix_type          submatrix)
{
    return banded_alignment_traceback<BAND_LEN,CHECKPOINTS>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        backtracer,
        checkpoints,
        submatrix,
        wfa );
}

template <
    uint32          BAND_LEN,
    uint32          MAX_PATTERN_LEN,
    uint32          CHECKPOINTS,
    typename        aligner_type,
    typename        pattern_string,
    typename        qual_string,
    typename        text_string,
    typename        backtracer_type,
    typename        wfa_type>
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
Alignment<int32> banded_alignment_traceback(
    const aligner_type      aligner,
    const pattern_string    pattern,
    const qual_string       quals,
    const text_string       text,
    const int32             min_score,
    backtracer_type&        backtracer,
    wfa_type&               wfa)
{
    typedef typename checkpoint_storage_type<aligner_type>::type checkpoint_type;
    typedef typename column_storage_type<aligner_type>::type     cell_type;

    const uint32 BITS = direction_vector_traits<aligner_type>::BITS;
    const uint32 ELEMENTS_PER_WORD = 32 / BITS;

    checkpoint_type checkpoints[ BAND_LEN*(1u + (MAX_PATTERN_LEN + CHECKPOINTS-1)/CHECKPOINTS) ];

    NVBIO_VAR_UNUSED const uint32 SUBMATRIX_WORDS = (BAND_LEN*CHECKPOINTS + ELEMENTS_PER_WORD-1) / ELEMENTS_PER_WORD;
    uint32 submatrix_base[ SUBMATRIX_WORDS ];
    PackedStream<uint32*,uint8,BITS,false> submatrix( submatrix_base );

#if !defined(__CUDA_ARCH__)
    // on the host, a wavefront aligner is answered by the affine banded DP of its own scheme,
    // which traces a real path and costs exactly what that path costs
    return banded_alignment_traceback_ondp<BAND_LEN, MAX_PATTERN_LEN, CHECKPOINTS>(
        typename aligner_type::aligner_tag(),
        aligner,
        pattern,
        quals,
        text,
        min_score,
        backtracer,
        wfa,
        checkpoints,
        submatrix );
#else
    return banded_alignment_traceback<BAND_LEN,CHECKPOINTS>(
        aligner,
        pattern,
        quals,
        text,
        min_score,
        backtracer,
        checkpoints,
        submatrix,
        wfa );
#endif
}



} // namespace alignment
} // namespace nvbio

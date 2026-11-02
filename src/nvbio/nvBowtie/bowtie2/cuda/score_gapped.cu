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

// NOTE: myers_opposite.h must precede every header pulling in score_opposite_inl.h:
// aln::max_text_gaps is called with a qualified name there, so overload lookup freezes
// at that definition point and would miss the Myers one.
#include <nvbio/alignment/myers_opposite.h>
#include <nvBowtie/bowtie2/cuda/score.h>
#include <nvBowtie/bowtie2/cuda/score_opposite_impl.h>
#include <nvBowtie/bowtie2/cuda/scoring.h>
#include <nvBowtie/bowtie2/cuda/params.h>
#include <nvbio/io/utils.h>
#include <type_traits>

namespace nvbio {
namespace bowtie2 {
namespace cuda {

namespace {

// tag-dispatch: only the edit-distance end-to-end aligner can be replaced by the
// Myers bit-vector core; every other scheme keeps the legacy path untouched
struct legacy_opposite_tag {};
struct myers_opposite_tag {};

template <typename aligner_type>
struct opposite_selector { typedef legacy_opposite_tag tag; };

template <typename algorithm_tag>
struct opposite_selector< aln::EditDistanceAligner<aln::SEMI_GLOBAL,algorithm_tag> > {
    typedef myers_opposite_tag tag;
};

template <typename pipeline_type, typename aligner_type>
void run_e2e_opposite(
    const pipeline_type&  pipeline,
    const aligner_type&   aligner,
    const ParamsPOD&      params,
    legacy_opposite_tag)
{
    detail::opposite_score_best(pipeline, aligner, params);
}

template <typename pipeline_type, typename aligner_type>
void run_e2e_opposite(
    const pipeline_type&  pipeline,
    const aligner_type&   aligner,
    const ParamsPOD&      params,
    myers_opposite_tag)
{
    if (params.opposite_myers) {
        typedef aln::MyersOppositeAligner<aln::SEMI_GLOBAL, typename aligner_type::algorithm_tag> myers_aligner;
        detail::opposite_score_best(pipeline, myers_aligner(), params);
    }
    else {
        detail::opposite_score_best(pipeline, aligner, params);
    }
}

} // anonymous namespace

//
// execute a batch of full-DP alignment score calculations for the opposite mates, best mapping
//
template <typename scheme_type>
void gapped_opposite_score_best_t(
    const BestApproxScoringPipelineState<scheme_type>&  pipeline,
    const ParamsPOD                                     params)
{
    if (params.alignment_type == LocalAlignment)
    {
        detail::opposite_score_best(
            pipeline,
            pipeline.scoring_scheme.local_aligner(),
            params );
    }
    else
    {
        typedef typename std::decay<decltype(pipeline.scoring_scheme.end_to_end_aligner())>::type e2e_aligner_type;
        run_e2e_opposite(
            pipeline,
            pipeline.scoring_scheme.end_to_end_aligner(),
            params,
            typename opposite_selector<e2e_aligner_type>::tag() );
    }
}

//
// execute a batch of full-DP alignment score calculations for the opposite mates, best mapping
//
// \b inputs:
//  - HitQueues::seed
//  - HitQueues::loc
//  - HitQueues::score
//  - HitQueues::sink
//
// \b outputs:
//  - HitQueues::opposite_score
//  - HitQueues::opposite_loc
//  - HitQueues::opposite_sink
//
void gapped_opposite_score_best(
    const BestApproxScoringPipelineState<EditDistanceScoringScheme>&    pipeline,
    const ParamsPOD&                                                    params)
{
    gapped_opposite_score_best_t( pipeline, params );
}

//
// execute a batch of full-DP alignment score calculations for the opposite mates, best mapping
//
// \b inputs:
//  - HitQueues::seed
//  - HitQueues::loc
//  - HitQueues::score
//  - HitQueues::sink
//
// \b outputs:
//  - HitQueues::opposite_score
//  - HitQueues::opposite_loc
//  - HitQueues::opposite_sink
//
void gapped_opposite_score_best(
    const BestApproxScoringPipelineState<SmithWatermanScoringScheme<> >&    pipeline,
    const ParamsPOD&                                                        params)
{
    gapped_opposite_score_best_t( pipeline, params );
}

//
// execute a batch of full-DP alignment score calculations for the opposite mates, best mapping
//
// \b inputs:
//  - HitQueues::seed
//  - HitQueues::loc
//  - HitQueues::score
//  - HitQueues::sink
//
// \b outputs:
//  - HitQueues::opposite_score
//  - HitQueues::opposite_loc
//  - HitQueues::opposite_sink
//
void gapped_opposite_score_best(
    const BestApproxScoringPipelineState<WfaScoringScheme<> >&              pipeline,
    const ParamsPOD&                                                        params)
{
    gapped_opposite_score_best_t( pipeline, params );
}


} // namespace cuda
} // namespace bowtie2
} // namespace nvbio

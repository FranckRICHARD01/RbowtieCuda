#pragma once

#include <nvbio/basic/types.h>
#include <nvbio/alignment/alignment_base.h>

// WFA
//#define WFA_TESTS
#define WFA_MIN (-32640)
#define WFA_MIN1 (-32639)
#define DPMATRIX_DIAGONAL_NULL 32640
#define WAVEFRONT_OFFSET_NULL (-20000)
#define DELTA 2
#ifdef WFA_TESTS
#define WFA_BAND_LEN_Y 600
#define DIM_SHARED 300
#define WFA_BAND_LEN2_LIMIT_BAND 450
#else
#define WFA_BAND_LEN_Y 60
#define DIM_SHARED 15
#define WFA_BAND_LEN2_LIMIT_BAND 50
#endif
#define DIM_SHARED2 (DIM_SHARED - DELTA)
#define COMPUTE_DIM_SHARED(score, band) ((score) * (2 * (band) + 1) + DELTA + 1)
// #define COMPUTE_DIM_SHARED(score, band) (score * score)
#define WFA_MAX_SCORE WFA_BAND_LEN_Y
#define WFA_BAND_LEN2_X WFA_BAND_LEN_Y
#define WFA_BAND_LEN2_Y WFA_BAND_LEN_Y
#define WFA_BAND_LEN2_CACHE (WFA_BAND_LEN_Y - 2)
#define WFA_BAND_LEN2_LIMIT (WFA_BAND_LEN_Y - 5)
#define DIM_Y_SHARED COMPUTE_DIM_SHARED(WFA_BAND_LEN_Y + 1, DIM_SHARED)
// #define DIM_Y_SHARED WFA_BAND_LEN2_X * WFA_BAND_LEN2_Y
#define WFA_LAUNCH 2
#define CORRECT_SIZE_MATRIX 1

typedef enum
{
    WFA_HEURISTIC_NONE = 0x0000000000000000ul,
    WFA_HEURISTIC_ADAPTIVE_BANDED_STATIC = 0x0000000000000001ul,
    WFA_HEURISTIC_ADAPTIVE_BANDED_ADAPTIVE = 0x0000000000000002ul,
    WFA_HEURISTIC_ADAPTIVE = 0x0000000000000004ul,
    WFA_HEURISTIC_XDROP = 0x0000000000000010ul,
    WFA_HEURISTIC_ZDROP = 0x0000000000000020ul,
    WFA_HEURISTIC_MASH = 0x0000000000000040ul,
} wfa_heuristic_type;

namespace nvbio
{
    namespace aln
    {
        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_wfadaptive(
            wfa_type &wfa,
            const score_type min_wavefront_length,
            const score_type max_distance_threshold,
            const score_type steps_between_cutoffs)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_ADAPTIVE;
            wfa.heuristic.min_wavefront_length = min_wavefront_length;
            wfa.heuristic.max_distance_threshold = max_distance_threshold;
            wfa.heuristic.steps_between_cutoffs = steps_between_cutoffs;
            wfa.heuristic.steps_wait = steps_between_cutoffs;

            wfa.heuristic.penalties_match = 0;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_wfmash(
            wfa_type &wfa,
            const score_type min_wavefront_length,
            const score_type max_distance_threshold,
            const score_type steps_between_cutoffs)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_MASH;
            wfa.heuristic.min_wavefront_length = min_wavefront_length;
            wfa.heuristic.max_distance_threshold = max_distance_threshold;
            wfa.heuristic.steps_between_cutoffs = steps_between_cutoffs;
            wfa.heuristic.steps_wait = steps_between_cutoffs;

            wfa.heuristic.penalties_match = 0;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_xdrop(
            wfa_type &wfa,
            const score_type xdrop,
            const score_type steps_between_cutoffs)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_XDROP;
            wfa.heuristic.xdrop = xdrop;
            wfa.heuristic.steps_between_cutoffs = steps_between_cutoffs;
            wfa.heuristic.steps_wait = steps_between_cutoffs;
            wfa.heuristic.max_sw_score = 0;
            wfa.heuristic.max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
            wfa.heuristic.max_sw_score_k = DPMATRIX_DIAGONAL_NULL;

            wfa.heuristic.penalties_match = 0;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_zdrop(
            wfa_type &wfa,
            const score_type zdrop,
            const score_type steps_between_cutoffs)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_ZDROP;
            wfa.heuristic.zdrop = zdrop;
            wfa.heuristic.steps_between_cutoffs = steps_between_cutoffs;
            wfa.heuristic.steps_wait = steps_between_cutoffs;
            wfa.heuristic.max_sw_score = 0;
            wfa.heuristic.max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
            wfa.heuristic.max_sw_score_k = DPMATRIX_DIAGONAL_NULL;

            wfa.heuristic.penalties_match = 0;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_banded_static(
            wfa_type &wfa,
            const score_type band_min_k,
            const score_type band_max_k)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_ADAPTIVE_BANDED_STATIC;
            wfa.heuristic.min_k = band_min_k;
            wfa.heuristic.max_k = band_max_k;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_set_banded_adaptive(
            wfa_type &wfa,
            const score_type band_min_k,
            const score_type band_max_k,
            const score_type steps_between_cutoffs)
        {
            wfa.heuristic.method |= WFA_HEURISTIC_ADAPTIVE_BANDED_ADAPTIVE;
            wfa.heuristic.min_k = band_min_k;
            wfa.heuristic.max_k = band_max_k;
            wfa.heuristic.steps_between_cutoffs = steps_between_cutoffs;
            // Internals
            wfa.heuristic.steps_wait = steps_between_cutoffs;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_clear(
            wfa_type &wfa)
        {
            wfa.heuristic.method = WFA_HEURISTIC_NONE;
            wfa.heuristic.steps_wait = wfa.heuristic.steps_between_cutoffs;
            wfa.heuristic.max_sw_score = 0;
            wfa.heuristic.max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
            wfa.heuristic.max_sw_score_k = DPMATRIX_DIAGONAL_NULL;

            wfa.heuristic.penalties_match = 0;
        }

        template <typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type wf_distance_end2end(
            const int16 offset,
            const score_type k,
            const score_type N,
            const score_type M)
        {
            const score_type left_v = N - (offset - k);
            const score_type left_h = M - offset;
            return (offset >= 0) ? nvbio::max(left_v, left_h) : -WAVEFRONT_OFFSET_NULL;
        }

        template <typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type wf_distance_end2end_weighted(
            const int16 offset,
            const score_type k,
            const score_type N,
            const score_type M,
            const score_type mfactor)
        {
            const score_type v = offset - k;
            const score_type h = offset;
            const score_type left_v = ((float)(N - v) / N * mfactor);
            const score_type left_h = ((float)(M - h) / M * mfactor);
            return (offset >= 0) ? nvbio::max(left_v, left_h) : -WAVEFRONT_OFFSET_NULL;
        }

        template <typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type wf_distance_endsfree(
            const int16 offset,
            const score_type k,
            const score_type N,
            const score_type M,
            const score_type pattern_end_free,
            const score_type text_end_free)
        {
            const score_type left_v = N - (offset - k);
            const score_type left_h = M - offset;
            const score_type left_v_endsfree = left_v - pattern_end_free;
            const score_type left_h_endsfree = left_h - text_end_free;
            const score_type dist_up = nvbio::max(left_h, left_v_endsfree);
            const score_type dist_down = nvbio::max(left_v, left_h_endsfree);
            return (offset >= 0) ? nvbio::min(dist_up, dist_down) : WAVEFRONT_OFFSET_NULL;
        }

        template <typename score_type, typename wfa_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type wf_compute_distance_end2end(
            wfa_type &wfa,
            const score_type score,
            const score_type N,
            const score_type M,
            int16 *distances)
        {
            const int16 *offsets = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            score_type k, min_distance = nvbio::max(N, M);
            const score_type lo = wfa.H_Band.get_lo(score);

#pragma unroll
            for (score_type k = lo; k <= wfa.H_Band.get_hi(score); ++k)
            {
                const score_type distance = wf_distance_end2end(offsets[k - lo], k, N, M);
                distances[k - lo + 1] = distance;
                min_distance = nvbio::min(min_distance, distance);
            }

            return min_distance;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type wf_compute_distance_end2end_weighted(
            wfa_type &wfa,
            const score_type score,
            const score_type N,
            const score_type M,
            int16 *distances)
        {
            const score_type mfactor = ((float)(N + M) / 2);
            const int16 *offsets = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            score_type k, min_distance = nvbio::max(N, M);
            const score_type lo = wfa.H_Band.get_lo(score);

#pragma unroll
            for (score_type k = lo; k <= wfa.H_Band.get_hi(score); ++k)
            {
                const score_type distance = wf_distance_end2end_weighted(offsets[k - lo], k, N, M, mfactor);
                distances[k - lo + 1] = distance;
                min_distance = nvbio::min(min_distance, distance);
            }

            return min_distance;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wf_heuristic_wfadaptice_reduce(
            wfa_type &wfa,
            const score_type score,
            const int16 *distances,
            const score_type min_distance,
            const score_type max_distance_threshold,
            const score_type min_k,
            const score_type max_k)
        {
            score_type k;

            const score_type top_limit = nvbio::min(max_k, wfa.H_Band.get_hi(score));
            const score_type lo = wfa.H_Band.get_lo(score);
            score_type lo_reduced = lo;
            for (k = lo; k < top_limit; ++k)
            {
                if (distances[k - lo + 1] - min_distance <= max_distance_threshold)
                    break;
                ++lo_reduced;
            }
            wfa.H_Band.set_lo(score, lo_reduced);

            const score_type bottom_limit = nvbio::max(min_k, wfa.H_Band.get_lo(score));
            score_type hi_reduced = wfa.H_Band.get_hi(score);
            for (k = wfa.H_Band.get_hi(score); k > bottom_limit; --k)
            {
                if (distances[k - lo + 1] - min_distance <= max_distance_threshold)
                    break;
                --hi_reduced;
            }
            wfa.H_Band.set_hi(score, hi_reduced);

            int32 pos = COMPUTE_DIM_SHARED(score, DIM_SHARED);

            if (lo != lo_reduced)
            {
                if (lo_reduced > lo)
                {
#pragma roll
                    for (score_type i = lo_reduced; i <= hi_reduced; i++)
                    {
                        wfa.H_Band.scores[pos + (i - lo_reduced)] = wfa.H_Band.scores[pos - (lo_reduced - lo) + (i - lo_reduced)];                      
                    }
                }
                else
                {
#pragma roll
                    for (score_type i = hi_reduced; i >= lo_reduced; i--)
                    {
                        wfa.H_Band.scores[pos + (i - lo_reduced)] = wfa.H_Band.scores[pos - (lo_reduced - lo) + (i - lo_reduced)]; 
                    }
                }
            }
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_wfadaptive(
            wfa_type &wfa,
            const score_type score,
            const score_type N,
            const score_type M,
            const bool wfmash_mode)
        {
            const score_type min_wavefront_length = wfa.heuristic.min_wavefront_length;
            const score_type max_distance_threshold = wfa.heuristic.max_distance_threshold;

            if (wfa.heuristic.steps_wait > 0)
                return;

            const score_type base_hi = wfa.H_Band.get_hi(score);
            const score_type base_lo = wfa.H_Band.get_lo(score);
            if (base_hi - base_lo + 1 < min_wavefront_length)
                return;

            int16 *H_Band_ptr = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            int16 distances[100]; //= &wfa.H_Band.scores[COMPUTE_DIM_SHARED(WFA_BAND_LEN2_CACHE, DIM_SHARED)];

#pragma unroll
            for (score_type k = base_lo - 1; k <= base_hi + 1; ++k)
            {
                distances[k - base_lo + 1] = H_Band_ptr[k - base_lo];
            }

            score_type min_distance;
            if (wfmash_mode)
            {
                min_distance = wf_compute_distance_end2end_weighted(wfa, score, N, M, distances);
            }
            else
            {
                min_distance = wf_compute_distance_end2end(wfa, score, N, M, distances);
            }

            const score_type alignment_k = M - N;
            wf_heuristic_wfadaptice_reduce(wfa, score, distances, min_distance, max_distance_threshold, alignment_k, alignment_k);

            wfa.heuristic.steps_wait = wfa.heuristic.steps_between_cutoffs;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wf_heuristic_compute_sw_scores(
            wfa_type &wfa,
            const score_type score,
            const score_type s,
            int16 *sw_scores,
            score_type &max_sw_score,
            score_type &max_k,
            score_type &max_offset)
        {
            const score_type wf_match = wfa.heuristic.penalties_match;
            const score_type swg_match = (wf_match != 0) ? (wfa.heuristic.penalties_match) : 1;

            const int16 *offsets = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            score_type k, cmax_sw_score = WFA_MIN1, cmax_k = 0, cmax_offset = 0;
            const score_type lo = wfa.H_Band.get_lo(score);

#pragma unroll
            for (score_type k = lo; k <= wfa.H_Band.get_hi(score); ++k)
            {
                const score_type offset = offsets[k - lo];
                if (offset < 0)
                    continue;
                const score_type v = offset - k;
                const score_type h = offset;
                const score_type sw_score = (swg_match * (v + h) - s) / 2;
                sw_scores[k - lo] = sw_score;
                if (cmax_sw_score < sw_score)
                {
                    cmax_sw_score = sw_score;
                    cmax_k = k;
                    cmax_offset = offset;
                }
            }

            max_sw_score = cmax_sw_score;
            max_k = cmax_k;
            max_offset = cmax_offset;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_xdrop(
            wfa_type &wfa,
            const score_type score,
            const score_type s)
        {
            // Check steps
            if (wfa.heuristic.steps_wait > 0)
                return;

            // Parameters
            const score_type base_hi = wfa.H_Band.get_hi(score);
            const score_type base_lo = wfa.H_Band.get_lo(score);

            int16 *offsets = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            int16 sw_scores[100];// = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(WFA_BAND_LEN2_CACHE, DIM_SHARED)];

            // Use victim as temporal buffer
#pragma unroll
            for (score_type k = base_lo - 1; k <= base_hi + 1; ++k)
            {
                sw_scores[k - base_lo + 1] = offsets[k - base_lo];
            }

            // Compute SW scores
            score_type cmax_sw_score, cmax_k, dummy;
            wf_heuristic_compute_sw_scores(
                wfa, score, s, sw_scores,
                cmax_sw_score, cmax_k, dummy);

            // Apply X-Drop
            const score_type xdrop = wfa.heuristic.xdrop;
            const score_type max_sw_score = wfa.heuristic.max_sw_score;
            if (wfa.heuristic.max_sw_score_k != DPMATRIX_DIAGONAL_NULL)
            {
                score_type lo = wfa.H_Band.get_lo(score);

                // Reduce from bottom
                score_type k;
                for (k = lo; k <= wfa.H_Band.get_hi(score); ++k)
                {
                    if (offsets[k] < 0)
                        continue;

                    if (max_sw_score - (score_type)sw_scores[k - lo] < xdrop)
                        break;
                }
                wfa.H_Band.set_lo(score, k);
                lo = k;
                // Reduce from top
                for (k = wfa.H_Band.get_hi(score); k >= lo; --k)
                {
                    if (offsets[k] < 0)
                        continue;

                    if (max_sw_score - (score_type)sw_scores[k - lo] < xdrop)
                        break;
                }
                wfa.H_Band.set_hi(score, k);
                // Update maximum score observed
                if (cmax_sw_score > wfa.heuristic.max_sw_score)
                {
                    wfa.heuristic.max_sw_score = cmax_sw_score;
                    wfa.heuristic.max_sw_score_k = cmax_k;
                }
            }
            else
            {
                // Update maximum score observed
                wfa.heuristic.max_sw_score = cmax_sw_score;
                wfa.heuristic.max_sw_score_k = cmax_k;
            }
            // Set wait steps (don't repeat this heuristic often)
            wfa.heuristic.steps_wait = wfa.heuristic.steps_between_cutoffs;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool wavefront_heuristic_zdrop(
            wfa_type &wfa,
            const score_type score,
            const score_type s)
        {
            if (wfa.heuristic.steps_wait > 0)
                return false;

            const score_type base_hi = wfa.H_Band.get_hi(score);
            const score_type base_lo = wfa.H_Band.get_lo(score);

            const int16 *H_Band_ptr = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)];
            int16 sw_scores[100];// = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(WFA_BAND_LEN2_CACHE, DIM_SHARED)];

#pragma unroll
            for (score_type k = base_lo - 1; k <= base_hi + 1; ++k)
            {
                sw_scores[k - base_lo + 1] = H_Band_ptr[k - base_lo];
            }

            score_type cmax_sw_score, cmax_k, cmax_offset;
            wf_heuristic_compute_sw_scores(wfa, score, s, sw_scores, cmax_sw_score, cmax_k, cmax_offset);

            const score_type zdrop = wfa.heuristic.zdrop;
            const score_type max_sw_score = wfa.heuristic.max_sw_score;
            const score_type max_k = wfa.heuristic.max_sw_score_k;
            const score_type max_offset = wfa.heuristic.max_sw_score_offset;
            if (max_k != DPMATRIX_DIAGONAL_NULL)
            {
                if (cmax_sw_score > wfa.heuristic.max_sw_score)
                {
                    wfa.heuristic.max_sw_score = cmax_sw_score;
                    wfa.heuristic.max_wf_score = score;
                    wfa.heuristic.max_sw_score_k = cmax_k;
                    wfa.heuristic.max_sw_score_offset = cmax_offset;
                }
                else
                {
                    if (max_sw_score - cmax_sw_score > zdrop)
                    {
                        wfa.heuristic.alignement_end_pos_score = wfa.heuristic.max_wf_score;
                        wfa.heuristic.alignement_end_pos_k = max_k;
                        wfa.heuristic.alignement_end_pos_offset = max_offset;

                        return true;
                    }
                }
            }
            else
            {
                wfa.heuristic.max_sw_score = cmax_sw_score;
                wfa.heuristic.max_wf_score = score;
                wfa.heuristic.max_sw_score_k = cmax_k;
                wfa.heuristic.max_sw_score_offset = cmax_offset;
            }

            wfa.heuristic.steps_wait = wfa.heuristic.steps_between_cutoffs;

            return false;
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_banded_static(
            wfa_type &wfa,
            const score_type score)
        {
            // Check wavefront limits
            if (wfa.H_Band.get_lo(score) < wfa.heuristic.min_k)
                wfa.H_Band.set_lo(score, wfa.heuristic.min_k);
            if (wfa.H_Band.get_hi(score) > wfa.heuristic.max_k)
                wfa.H_Band.set_hi(score, wfa.heuristic.max_k);
        }

        template <typename wfa_type, typename score_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wavefront_heuristic_banded_adaptive(
            wfa_type &wfa,
            const score_type score,
            const score_type pattern_length,
            const score_type text_length)
        {
            // Check steps
            if (wfa.heuristic.steps_wait > 0)
                return;
            // Check wavefront length
            const score_type hi = wfa.H_Band.get_hi(score);
            const score_type lo = wfa.H_Band.get_lo(score);
            const score_type wf_length = hi - lo + 1;
            if (wf_length < 4)
                return; // We cannot do anything here
            // Adjust the band
            uint32 pos = COMPUTE_DIM_SHARED(score, DIM_SHARED);
            const int16 *offsets = &wfa.H_Band.scores[pos] - lo;
            const score_type max_wf_length = wfa.heuristic.max_k - wfa.heuristic.min_k + 1;
            if (wf_length > max_wf_length)
            {
                // Sample wavefront
                const score_type leeway = (wf_length - max_wf_length) / 2;
                const score_type quarter = wf_length / 4;
                const score_type dist_p0 = wf_distance_end2end(
                    offsets[lo], lo, pattern_length, text_length);
                const score_type dist_p1 = wf_distance_end2end(
                    offsets[lo + quarter], (score_type)(lo + quarter), pattern_length, text_length);
                const score_type dist_p2 = wf_distance_end2end(
                    offsets[lo + 2 * quarter], (score_type)(lo + 2 * quarter), pattern_length, text_length);
                const score_type dist_p3 = wf_distance_end2end(
                    offsets[hi], hi, pattern_length, text_length);
                // Heuristically decide where to place the band
                score_type new_lo = lo;
                if (dist_p0 > dist_p3)
                    new_lo += leeway;
                if (dist_p1 > dist_p2)
                    new_lo += leeway;
                // Set wavefront limits
                wfa.H_Band.set_lo(score, new_lo);
                if (wfa.H_Band.get_lo(score) < lo)
                    wfa.H_Band.set_lo(score, lo);
                wfa.H_Band.set_hi(score, new_lo + max_wf_length - 1);
                if (wfa.H_Band.get_hi(score) > hi)
                    wfa.H_Band.set_hi(score, hi);

                new_lo = wfa.H_Band.get_lo(score);
                score_type new_hi = wfa.H_Band.get_hi(score);

                if (lo != new_lo)
                {
                    if (new_lo > lo)
                    {
#pragma roll
                        for (score_type i = new_lo; i <= new_hi; i++)
                        {
                            wfa.H_Band.scores[pos + (i - new_lo)] = wfa.H_Band.scores[pos + (new_lo - lo) + (i - new_lo)];
                        }
                    }
                    else
                    {
#pragma roll
                        for (score_type i = new_hi; i >= new_lo; i--)
                        {
                            wfa.H_Band.scores[pos + (i - new_lo)] = wfa.H_Band.scores[pos + (new_lo - lo) + (i - new_lo)];
                        }
                    }
                }
            }
            // Set wait steps (don't repeat this heuristic often)
            wfa.heuristic.steps_wait = wfa.heuristic.steps_between_cutoffs;
        }

        template <
            typename score_type,
            typename wfa_type>
        NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool wavefront_heuristic_cutoff(
            wfa_type &wfa,
            const score_type score,
            const score_type s,
            const score_type N,
            const score_type M)
        {
            int16 *H_Band_ptr = &wfa.H_Band.scores[COMPUTE_DIM_SHARED(score, DIM_SHARED)] - wfa.H_Band.get_lo(score);

            if (!H_Band_ptr || wfa.H_Band.get_null(score) || wfa.H_Band.get_lo(score) > wfa.H_Band.get_hi(score))
                return false;

            wfa.heuristic.steps_wait--;

            const score_type hi_base = wfa.H_Band.get_hi(score);
            const score_type lo_base = wfa.H_Band.get_lo(score);

            if (wfa.heuristic.method & WFA_HEURISTIC_ADAPTIVE)
            {
                wavefront_heuristic_wfadaptive(wfa, score, N, M, false);
            }

            if (wfa.heuristic.method & WFA_HEURISTIC_MASH)
            {
                wavefront_heuristic_wfadaptive(wfa, score, N, M, true);
            }

            if (wfa.heuristic.method & WFA_HEURISTIC_XDROP)
            {
                wavefront_heuristic_xdrop(wfa, score, s);
            }

            if (wfa.heuristic.method & WFA_HEURISTIC_ZDROP)
            {
                if (wavefront_heuristic_zdrop(wfa, score, s)) return true;
            }

            if (wfa.heuristic.method & WFA_HEURISTIC_ADAPTIVE_BANDED_STATIC)
            {
                wavefront_heuristic_banded_static(wfa, score);
            }

            if (wfa.heuristic.method & WFA_HEURISTIC_ADAPTIVE_BANDED_ADAPTIVE)
            {
                wavefront_heuristic_banded_adaptive(wfa, score, N, M);
            }

            if (lo_base == wfa.H_Band.get_lo(score) && hi_base == wfa.H_Band.get_hi(score))
                return false;

            if (wfa.H_Band.get_lo(score) > wfa.H_Band.get_hi(score))
            {
                wfa.H_Band.set_hi(score, 0);
                wfa.H_Band.set_lo(score, 0);
                wfa.H_Band.set_null(score, true);
            }

            wfa.wf_heuristic_equate(wfa.E_Band, wfa.H_Band, score);
            wfa.wf_heuristic_equate(wfa.F_Band, wfa.H_Band, score);

            return false;
        }

        template <typename score_type>
        struct wfa_type
        {
            struct defHeuristic
            {
                score_type min_wavefront_length;
                score_type max_distance_threshold;
                score_type steps_between_cutoffs;
                score_type steps_wait;
                score_type xdrop;
                score_type zdrop;
                score_type max_sw_score;
                score_type max_sw_score_k;
                score_type max_sw_score_offset;
                score_type alignement_end_pos_k;
                score_type alignement_end_pos_score;
                score_type alignement_end_pos_offset;
                score_type penalties_match;
                score_type max_wf_score;
                score_type min_k;
                score_type max_k;
                unsigned long method;

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                defHeuristic()
                {
                    min_wavefront_length = 0;
                    max_distance_threshold = 0;
                    steps_between_cutoffs = 0;
                    steps_wait = 1;
                    xdrop = 0;
                    zdrop = 0;
                    max_sw_score = 0;
                    max_sw_score_k = 0;
                    max_sw_score_offset = 0;
                    alignement_end_pos_k = 0;
                    alignement_end_pos_score = 0;
                    alignement_end_pos_offset = 0;
                    penalties_match = 0;
                    max_wf_score = 0;
                    min_k = 0;
                    max_k = 0;
                    method = WFA_HEURISTIC_NONE;
                }
            };

            struct defAlignment
            {
                int16 *scores = nullptr;
                int16 *hi = nullptr;
                int16 *lo = nullptr;
                bool  *null = nullptr;

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                defAlignment()
                {
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
                set_scores_data(
                    int16 *data)
                {
                    scores = data;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
                set_lo_data(
                    int16 *data)
                {
                    lo = data;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
                set_hi_data(
                    int16 *data)
                {
                    hi = data;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
                set_null_data(
                    bool *data)
                {
                    null = data;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    score_type
                    get_lo(
                        int32 i)
                {
                    //if (i < 0)
                    //    return 1;

                                       

                    NVBIO_CUDA_DEBUG_ASSERT(i >= 0 && i < WFA_BAND_LEN2_Y, "get_lo() indice problem i=%d\n", i);
                    NVBIO_CUDA_DEBUG_ASSERT(lo[i] >= -1000 && lo[i] < 1000, "get_lo() return problem  i=%i get_lo=%i\n", i, lo[i]);

                    return lo[i];
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    score_type
                    get_hi(
                        int32 i)
                {
                    //if (i < 0)
                    //    return -1;

                    NVBIO_CUDA_DEBUG_ASSERT(i >= 0 && i < WFA_BAND_LEN2_Y, "get_hi() indice problem i=%d\n", i);
                    NVBIO_CUDA_DEBUG_ASSERT(hi[i] >= -1000 && hi[i] < 1000, "get_hi() return problem  i=%i get_hi=%i\n", i, hi[i]);

                    return hi[i];
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    bool
                    get_null(
                         int32 i)
                {
                    if (i < 0)
                        return true;

                    NVBIO_CUDA_DEBUG_ASSERT(i >= 0 && i < WFA_BAND_LEN2_Y, "get_null() indice problem i=%d\n", i);

                    return null[i];
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void set_lo(
                    int32 i,
                    score_type v)
                {
#if __CUDA_ARCH__ >= 800
                    __stwt(&lo[i], v);
#else
                    lo[i] = v;
#endif
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void set_hi(
                    int32 i,
                    score_type v)
                {
#if __CUDA_ARCH__ >= 800
                    __stwt(&hi[i], v);
#else
                    hi[i] = v;
#endif
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void set_null(
                    int32 i,
                    bool v)
                {
                    null[i] = v;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
                    score_type
                    get(
                        int32 i,
                        int32 j)
                {
                    if (i < 0)
                        return WFA_MIN;

                    if (j > hi[i] || j < lo[i])
                        return WFA_MIN;

                    int32 pos = COMPUTE_DIM_SHARED(i, DIM_SHARED) - lo[i];

                    // NVBIO_CUDA_DEBUG_ASSERT(pos < WFA_BAND_LEN2_X * WFA_BAND_LEN2_Y, "get() indice problem pos < BAND_LEN2_X * BAND_LEN2_Y i=%i j=%i pos=%i ^=%d\n", i, j, pos, WFA_BAND_LEN2_X * WFA_BAND_LEN2_Y);

                    // NVBIO_CUDA_DEBUG_ASSERT(scores[pos] >= 0, "get() return negative value i=%i j=%i pos=%i val=%d\n", i, j, pos, scores[pos]);

                    /*if (scores[pos] > WFA_BAND_LEN2_Y || (scores[pos] < 0 && scores[pos] != WFA_MIN && scores[pos] != WFA_MIN + 1))
                    {
                          NVBIO_CUDA_DEBUG_ASSERT(pos < 0, "test i=%i j=%i pos=%i ^=%d\n", i, j, pos, scores[pos]);
                    }*/

                    return scores[pos + j];
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void set(
                    int32 i,
                    int32 j,
                    score_type s)
                {
                    int32 pos = COMPUTE_DIM_SHARED(i, DIM_SHARED) - lo[i];

                    // NVBIO_CUDA_DEBUG_ASSERT(i >= 0, "set() indice problem i=%d\n", i);
                    // NVBIO_CUDA_DEBUG_ASSERT(pos >= -(i + 1), "set() indice problem pos >= 0 i=%i j=%i pos=%i\n", i, j, pos);
                    // NVBIO_CUDA_DEBUG_ASSERT(pos < DIM_Y_SHARED, "set() indice problem pos < DIM_Y_SHARED i=%i j=%i pos=%i ^=%d\n", i, j, pos, DIM_Y_SHARED);

#if __CUDA_ARCH__ >= 800
                    __stwt(&scores[pos + j], s);
#else
                    scores[pos + j] = s;
#endif
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void incrementation(
                    int32 i,
                    int32 j)
                {
                    // int32 h = i;
                    int32 pos = COMPUTE_DIM_SHARED(i, DIM_SHARED) - lo[i];

                    //     NVBIO_CUDA_DEBUG_ASSERT(i >= 0, "incrementation() indice problem i >= 0, i=%d \n", i);
                    //    NVBIO_CUDA_DEBUG_ASSERT(pos + j >= -(h + 1), "incrementation() indice problem pos + j >= 0 i=%i j=%i pos=%i\n", i, j, pos);
                    //    NVBIO_CUDA_DEBUG_ASSERT(pos + j < WFA_BAND_LEN2_X * WFA_BAND_LEN2_Y, "incrementation() indice problem pos + j < BAND_LEN2_X * BAND_LEN2_Y j=%i pos=%i ^=%d\n", j, pos, WFA_BAND_LEN2_X * WFA_BAND_LEN2_Y);

                    scores[pos + j]++;
                }

                NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
                set_hi_lo(int32 i, bool H)
                {
                    if (H)
                    {

#if __CUDA_ARCH__ >= 800
                        __stwt(&lo[i], 0);
                        __stwt(&hi[i], 0);
#else
                        lo[i] = 0;
                        hi[i] = 0;
#endif
                    }
                    else
                    {

#if __CUDA_ARCH__ >= 800
                        __stwt(&lo[i], 0);
                        __stwt(&hi[i], 0);
#else
                        lo[i] = 0;
                        hi[i] = 0;
#endif
                    }
                }
            };

            // Scores
            defAlignment H_Band;
            defAlignment F_Band;
            defAlignment E_Band;

            int16 *Point_H_BAND;

            defHeuristic heuristic;

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
            set_pointH_data(
                int16 *data)
            {
                Point_H_BAND = data;
            }

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE void
            set_HEF_hi_lo(
                int32 ss,
                score_type hi,
                score_type lo)
            {
#if __CUDA_ARCH__ >= 800
                __stwt(&H_Band.lo[ss], lo);
                __stwt(&H_Band.hi[ss], hi);

                __stwt(&E_Band.lo[ss], lo);
                __stwt(&E_Band.hi[ss], hi);

                __stwt(&F_Band.lo[ss], lo);
                __stwt(&F_Band.hi[ss], hi);
#else

                H_Band.lo[ss] = lo;
                H_Band.hi[ss] = hi;

                E_Band.lo[ss] = lo;
                E_Band.hi[ss] = hi;

                F_Band.lo[ss] = lo;
                F_Band.hi[ss] = hi;
#endif
            }

            template <typename context_type, typename ref_type, typename text_cache_type, typename scoring_type>
            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void add_backtrace_etape(
                context_type context,
                score_type x,
                score_type y,
                score_type h,
                score_type v,
                const uint32 block,
                const ref_type ref,
                const text_cache_type q_cache,
                const DirectionVector d1,
                const DirectionVector d2,
                const DirectionVector d3,
                const scoring_type scoring)
            {
                // if (CHECK_M == false)

                NVBIO_CUDA_DEBUG_ASSERT(x >= 0 && x < WFA_BAND_LEN2_Y, "add_backtrace_etape indice problem x=%d block=%d (x=%d y=%d h=%d v=%d)\n", x, block, x, y, h, v);
                // NVBIO_CUDA_DEBUG_ASSERT(y >= 0 && y < WFA_BAND_LEN2_Y, "add_backtrace_etape indice problem y=%d block=%d (x=%d y=%d h=%d v=%d)\n", y, block, x, y, h, v);
                // NVBIO_CUDA_DEBUG_ASSERT(h >= 0 && h < WFA_BAND_LEN2_Y, "add_backtrace_etape indice problem h=%d block=%d (x=%d y=%d h=%d v=%d)\n", h, block, x, y, h, v);
                NVBIO_CUDA_DEBUG_ASSERT(v >= 0 && v < WFA_BAND_LEN2_Y, "add_backtrace_etape indice problem v=%d block=%d (x=%d y=%d h=%d v=%d)\n", v, block, x, y, h, v);
                NVBIO_CUDA_DEBUG_ASSERT(block < 1, "add_backtrace_etape indice problem block=%d\n", block);

                // if (x >= 0 && y >= 0 && h >= 0 && v >= 0)
                {
                    context.new_cell(
                        x,
                        y,
                        scoring.substitution(h, block + v, ref[h], q_cache[block + v].x, q_cache[block + v].y),
                        d1,
                        d2,
                        d3);
                }
            }

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static score_type f(
                score_type i,
                int16 *mat)
            {
                /*score_type k = i - 10;

                if (k < 0)
                    k = 0;

                while (mat[k] < i)
                    k++;

                return (mat[k] == i) ? k : WFA_MIN;*/

                return i;
            }

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void trim_ends(
                defAlignment *mat,
                score_type s,
                uint32 length_ref,
                uint32 length_str)
            {
                int16 lo = mat->get_lo(s);
                int16 lo_mem = lo;
                int16 k;

                int32 pos = COMPUTE_DIM_SHARED(s, DIM_SHARED);

                for (k = mat->get_hi(s); k >= lo; --k)
                {
                    uint16 offset = mat->scores[pos - lo + k];

                    if (offset <= length_str && (uint16)(offset - k) <= length_ref)
                        break;
                }

                mat->set_hi(s, k);
                int16 hi = k;

                for (k = mat->get_lo(s); k <= hi; ++k)
                {
                    uint16 offset = mat->scores[pos - lo + k];

                    if (offset <= length_str && (uint16)(offset - k) <= length_ref)
                        break;
                }

                mat->set_lo(s, k);

                if (mat->get_lo(s) > mat->get_hi(s))
                {
                    mat->set_lo(s, 0);
                    mat->set_hi(s, 0);
                    if (s == 0)
                        mat->set(s, 0, WFA_MIN);
                    else
                        mat->set_null(s, true);
                }

                lo = mat->get_lo(s);
                hi = mat->get_hi(s);

                if (lo != lo_mem)
                {
                    if (lo > lo_mem)
                    {
                        for (score_type i = lo; i <= hi; i++)
                        {
                            mat->scores[pos + (i - lo)] = mat->scores[pos + (lo - lo_mem) + (i - lo)];                            
                        }
                    }
                    else
                    {
                        for (score_type i = hi; i >= lo; i--)
                        {
                            mat->scores[pos + (i - lo)] = mat->scores[pos + (lo - lo_mem) + (i - lo)];
                        }
                    }
                }
            }

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static void wf_heuristic_equate(
                defAlignment &dst,
                defAlignment &src,
                const score_type score)
            {
                score_type lo_mem = dst.get_lo(score);
                score_type hi_mem = dst.get_hi(score);
                
                if (!dst.get_null(score))
                {
                    if (src.get_lo(score) > dst.get_lo(score))
                        dst.set_lo(score, src.get_lo(score));
                    if (src.get_hi(score) < dst.get_hi(score))
                        dst.set_hi(score, src.get_hi(score));
                    if (dst.get_lo(score) > dst.get_hi(score))
                    {
                        dst.set_lo(score, 0);
                        dst.set_hi(score, 0);
                        dst.set_null(score, true);
                    }

                    uint32 pos = COMPUTE_DIM_SHARED(score, DIM_SHARED);

                    score_type lo = dst.get_lo(score);
                    score_type hi = dst.get_hi(score);

                    if (lo != lo_mem)
                    {
                        if (lo > lo_mem)
                        {
#pragma roll
                            for (score_type i = lo; i <= hi; i++)
                            {
                                dst.scores[pos + (i - lo)] = dst.scores[pos + (lo - lo_mem) + (i - lo)];                 
                            }
                        }
                        else
                        {
#pragma roll
                            for (score_type i = hi; i >= lo; i--)
                            {
                                dst.scores[pos + (i - lo)] = dst.scores[pos + (lo - lo_mem) + (i - lo)]; 
                            }
                        }
                    }
                }
            }

            NVBIO_FORCEINLINE NVBIO_HOST_DEVICE static bool testEnd(
                score_type length_ref,
                score_type length_str,
                defAlignment mat,
                score_type s,
                score_type alignment_k,
                score_type alignment_offset,
                bool matrix_H)
            {
                //const bool test = !matrix_H && mat.get_null(s);

                if (/*!test ||*/ mat.get_lo(s) > alignment_k || alignment_k > mat.get_hi(s))
                    return false;

                if (mat.get(s, alignment_k) < alignment_offset)
                    return false;

                return true;
            }
        };
    };
};

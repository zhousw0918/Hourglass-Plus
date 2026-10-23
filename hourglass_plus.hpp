#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "hourglass.cpp"

enum class SegmentBranch {
    kUncorrelated,
    kBinary,
    kBitmap,
};

struct Segment {
    size_t start_idx = 0;
    size_t end_idx = 0;
    uint64_t min_key = 0;
    uint64_t max_key = 0;

    double local_corr = 0.0;
    size_t optimal_ls = 0;
    size_t fp_len = 0;
    double prefix_fpr = 1.0;
    double suffix_fpr = 0.0;
    double optimal_fpr = 1.0;
    size_t query_count = 0;
    size_t unique_prefix_count = 0;

    SegmentBranch branch = SegmentBranch::kUncorrelated;
    double curvature = 0.0;
    double ls_slope = 0.0;
    double marginal_gain = 0.0;

    std::vector<std::pair<uint64_t, uint64_t>> local_queries;

    size_t num_keys() const { return end_idx - start_idx; }
};

struct SegmentModelEstimate {
    size_t unique_prefix_count = 0;
    size_t optimal_ls = 0;
    size_t fp_len = 0;
    double prefix_fpr = 1.0;
    double suffix_fpr = 0.0;
    double estimated_fpr = 1.0;
    SegmentBranch branch = SegmentBranch::kUncorrelated;
    double curvature = 0.0;
    double ls_slope = 0.0;
    double marginal_gain = 0.0;
};

inline uint64_t SafeHourglassPayloadSizeInBits(Hourglass& hourglass) {
    uint64_t total_bits = hourglass.kBytesTags_ * 8;
    if (hourglass.correlated_) {
        total_bits += hourglass.kBytesKeys_ * 8 + hourglass.rm.bitCount();
    }
#ifdef ADAPT
    if (hourglass.fingerprint_extractor_ != nullptr) {
        total_bits += hourglass.fingerprint_extractor_->bitCount();
    }
#endif
    return total_bits;
}

inline uint64_t SafeHourglassSizeInBits(Hourglass& hourglass) {
    return SafeHourglassPayloadSizeInBits(hourglass) + sizeof(hourglass) * 8;
}

template <typename t_itr>
inline double EstimateLocalCorrValue(const t_itr begin, const t_itr end,
                                     const std::vector<std::pair<uint64_t, uint64_t>>& queries,
                                     size_t range_size) {
    if (queries.empty() || begin == end) return 0.0;

    const uint64_t interval = 1ULL << range_size;
    uint64_t corr = 0;

    for (const auto& query : queries) {
        const uint64_t left = query.first;
        const uint64_t right = query.second;
        uint64_t distance = 0;

        auto low = std::lower_bound(begin, end, left);
        if (low == end) {
            distance = left - *(end - 1);
        } else {
            const uint64_t dist_left =
                (low == begin) ? std::numeric_limits<uint64_t>::max() : (left - *(low - 1));
            const uint64_t dist_right = (*low <= right) ? 0 : (*low - right);
            distance = std::min(dist_left, dist_right);
        }

        if (distance <= interval) {
            ++corr;
        }
    }

    if (corr == 0) return 0.0;
    const double D = 1.0 - std::log2(queries.size() / static_cast<double>(corr) * interval) / 30.0;
    return std::max(0.0, std::min(1.0, D));
}

template <typename t_itr>
inline uint64_t CountCorrelatedQueries(const t_itr begin, const t_itr end,
                                       const std::vector<std::pair<uint64_t, uint64_t>>& queries,
                                       size_t range_size) {
    if (queries.empty() || begin == end) return 0;
    const uint64_t interval = 1ULL << range_size;
    uint64_t corr = 0;
    for (const auto& query : queries) {
        const uint64_t left = query.first, right = query.second;
        uint64_t distance = 0;
        auto low = std::lower_bound(begin, end, left);
        if (low == end) {
            distance = left - *(end - 1);
        } else {
            const uint64_t dist_left =
                (low == begin) ? std::numeric_limits<uint64_t>::max() : (left - *(low - 1));
            const uint64_t dist_right = (*low <= right) ? 0 : (*low - right);
            distance = std::min(dist_left, dist_right);
        }
        if (distance <= interval) ++corr;
    }
    return corr;
}

constexpr double kShrinkTestLevel = 0.01;
inline double PoissonLowerTail(double lambda, double c) {
    double term = std::exp(-lambda), sum = term;
    for (int k = 1; k <= static_cast<int>(c); ++k) { term *= lambda / k; sum += term; }
    return sum;
}
template <typename t_itr>
inline double EstimateLocalCorrValueShrunk(const t_itr begin, const t_itr end,
                                           const std::vector<std::pair<uint64_t, uint64_t>>& queries,
                                           size_t range_size, double global_rate) {
    if (begin == end) return 0.0;
    const double interval = static_cast<double>(1ULL << range_size);
    const double c = static_cast<double>(CountCorrelatedQueries(begin, end, queries, range_size));
    const double n = static_cast<double>(queries.size());
    double rate = 0.0;
    if (global_rate > 0.0 && PoissonLowerTail(n * global_rate, c) >= kShrinkTestLevel) {
        const double alpha = 1.0 / global_rate;
        rate = (c + alpha * global_rate) / (n + alpha);
    } else if (n > 0.0) {
        rate = c / n;
    }
    if (rate <= 0.0) return 0.0;
    const double D = 1.0 - std::log2(interval / rate) / 30.0;
    return std::max(0.0, std::min(1.0, D));
}

template <typename t_itr>
inline size_t CountUniquePrefixes(const t_itr begin, const t_itr end, size_t suffix_bits) {
    if (begin == end) return 0;

    uint64_t prev_prefix = *begin >> suffix_bits;
    size_t count = 1;
    for (auto it = begin + 1; it != end; ++it) {
        const uint64_t cur_prefix = *it >> suffix_bits;
        if (cur_prefix != prev_prefix) {
            prev_prefix = cur_prefix;
            ++count;
        }
    }
    return count;
}

inline size_t GetNextPow2(size_t x) {
    size_t power = 1;
    while (power < x) power <<= 1;
    return power;
}

inline size_t GetPhysicalInstances(size_t num_keys, size_t base_block_size) {
    if (base_block_size == 0 || num_keys == 0) return 0;
    const size_t full_blocks = num_keys / base_block_size;
    const size_t tail = num_keys % base_block_size;
    return static_cast<size_t>(__builtin_popcountll(full_blocks)) + (tail > 0 ? 1 : 0);
}

inline size_t ClampFingerprintLength(double fp_len) {
    size_t rounded = static_cast<size_t>(fp_len);
    if (rounded > 32) return 32;
    if (rounded <= 5) return 6;
    return rounded;
}

inline double EstimatePrefixFPR(double fp_len) {
    if (fp_len <= 0.0) return 1.0;
    return std::min(1.0, std::pow(2.0, 2.0 - fp_len));
}

inline double EstimateSuffixFPR(double D, size_t optimal_ls, size_t max_range_size) {
    if (optimal_ls >= max_range_size) return 0.0;
    const double exponent =
        static_cast<double>(max_range_size - optimal_ls) - 30.0 * (1.0 - D);
    return std::min(1.0, std::pow(2.0, exponent));
}

inline double EstimateAverageKeysPerPrefix(size_t num_keys, size_t unique_prefix_count) {
    if (num_keys == 0 || unique_prefix_count == 0) return 0.0;
    return static_cast<double>(num_keys) / static_cast<double>(unique_prefix_count);
}

inline double EstimateTheoryLoadRatioW() {
    return 0.95 / 4.0;
}

inline double EstimateBinaryBitmapThreshold(size_t optimal_ls) {
    if (optimal_ls == 0) return 0.0;
    return std::pow(2.0, static_cast<double>(optimal_ls))
         / static_cast<double>(optimal_ls);
}

inline bool UseBinaryEncoding(size_t optimal_ls, double avg_keys_per_prefix) {
    if (optimal_ls == 0) return false;
    return avg_keys_per_prefix <= EstimateBinaryBitmapThreshold(optimal_ls);
}

inline SegmentBranch ClassifySegmentBranch(size_t num_keys,
                                           size_t unique_prefix_count,
                                           size_t optimal_ls) {
    if (optimal_ls == 0 || num_keys == 0 || unique_prefix_count == 0) {
        return SegmentBranch::kUncorrelated;
    }
    return UseBinaryEncoding(optimal_ls, EstimateAverageKeysPerPrefix(num_keys, unique_prefix_count))
        ? SegmentBranch::kBinary
        : SegmentBranch::kBitmap;
}

inline size_t CalculateOriginalHourglassLs(double D, double bpk, size_t num_keys,
                                           size_t unique_prefix_count, size_t max_range_size) {
    if (num_keys == 0 || unique_prefix_count == 0) return 0;

    const int l_s = static_cast<int>(max_range_size);
    const int l_q = static_cast<int>(std::floor(30.0 * (1.0 - D + 0.01)));
    const double rounded_bpk = std::round(bpk);
    const double avg = static_cast<double>(num_keys) / static_cast<double>(unique_prefix_count);

    int l = static_cast<int>(l_s - l_q - 2 + 0.95 * avg * (rounded_bpk - 1.0)
                             - std::log2(0.95 * avg));
    l = std::max(0, std::min(l_s, l));

    if (l == 0 || avg <= std::pow(2.0, static_cast<double>(l)) / static_cast<double>(l)) {
        if (avg * (rounded_bpk - l) < 4.0) {
            const double fallback =
                rounded_bpk - static_cast<double>(unique_prefix_count)
                            / static_cast<double>(num_keys) * 4.0;
            l = fallback > 0.0 ? static_cast<int>(fallback) : 0;
        }
        return static_cast<size_t>(std::max(0, std::min(l_s, l)));
    }

    double gap = 2000.0;
    for (size_t i = 0; i <= max_range_size; ++i) {
        const double derivative = 0.95 * std::pow(2.0, static_cast<double>(i))
                                + 2.0 * static_cast<double>(i)
                                - static_cast<double>(l_s)
                                - static_cast<double>(l_q)
                                - 2.0
                                + 0.95 * rounded_bpk * avg
                                + 0.53;
        if (std::abs(derivative) < gap) {
            l = static_cast<int>(i);
            gap = std::abs(derivative);
        }
    }
    return static_cast<size_t>(std::max(0, std::min(l_s, l)));
}

inline size_t CalculateOptimalLs(double D, double bpk, size_t num_keys,
                                 size_t unique_prefix_count, size_t max_range_size) {
    return CalculateOriginalHourglassLs(D, bpk, num_keys, unique_prefix_count, max_range_size);
}

template <typename t_itr>
inline double EstimateFingerprintLengthRaw(const t_itr begin, const t_itr end,
                                           size_t optimal_ls, size_t max_range_size,
                                           size_t unique_prefix_count, double bpk) {
    const size_t num_keys = std::distance(begin, end);
    if (num_keys == 0) return 6.0;

    const size_t rounded_bpk = static_cast<size_t>(std::llround(bpk));
    if (optimal_ls == 0) {
        return static_cast<double>(rounded_bpk + 1);
    }

    const size_t pruned_layers = max_range_size - optimal_ls;
    const size_t threshold =
        static_cast<size_t>(std::pow(2.0, static_cast<double>(optimal_ls)))
        / optimal_ls + 1;
    uint64_t key_bits = 0;
    size_t max_suffix_count = 1;
    std::vector<size_t> num_ones_len(optimal_ls + 1, 0);

    uint64_t prev_prefix = *begin >> max_range_size;
    uint16_t prev_suffix =
        static_cast<uint16_t>((*begin >> pruned_layers) & ((1ULL << optimal_ls) - 1));
    size_t suffix_count = 1;

    auto finalize_prefix = [&](size_t count) {
        max_suffix_count = std::max(max_suffix_count, count);
        const size_t bits_needed =
            static_cast<size_t>(std::ceil(std::log2(static_cast<double>(count) + 1.0)));
        for (size_t i = 0; i < bits_needed; ++i) {
            num_ones_len[i] += (count >> i) & 0x1;
        }
        key_bits += std::min(count, threshold) * optimal_ls;
    };

    for (auto it = begin + 1; it != end; ++it) {
        const uint64_t key = *it;
        const uint64_t cur_prefix = key >> max_range_size;
        const uint16_t cur_suffix =
            static_cast<uint16_t>((key >> pruned_layers) & ((1ULL << optimal_ls) - 1));

        if (cur_prefix != prev_prefix) {
            finalize_prefix(suffix_count);
            prev_prefix = cur_prefix;
            prev_suffix = cur_suffix;
            suffix_count = 1;
        } else if (cur_suffix != prev_suffix) {
            ++suffix_count;
            prev_suffix = cur_suffix;
        }
    }
    finalize_prefix(suffix_count);

    const size_t assoc = 4;
    size_t num_buckets = GetNextPow2(std::max<size_t>(1, unique_prefix_count / assoc));
    const double frac =
        static_cast<double>(unique_prefix_count) / static_cast<double>(num_buckets * assoc);
    if (frac > 0.95) {
        num_buckets <<= 1;
    }

    const size_t num_entries = num_buckets * assoc;
    const size_t keys_tag_len =
        static_cast<size_t>(std::ceil(std::log2(static_cast<double>(max_suffix_count) + 1.0)));

    double len_bits = 0.0;
    for (size_t i = 0; i < keys_tag_len; ++i) {
        double proportion = num_ones_len[i] / static_cast<double>(num_entries);
        proportion = proportion < 0.7 ? proportion : 1.0 - proportion;
        if (proportion > 0.3) {
            len_bits += 1.25 * static_cast<double>(num_entries);
        } else {
            len_bits += (34.72 * std::pow(proportion, 3)
                      - 21.85 * std::pow(proportion, 2)
                      + 7.23 * proportion
                      + 0.02) * static_cast<double>(num_entries);
        }
    }

    const double total_budget_bits = static_cast<double>(rounded_bpk) * static_cast<double>(num_keys);
    double fp_len = 0.0;
    if (total_budget_bits >= static_cast<double>(key_bits) + len_bits) {
        const double prefix_bits = total_budget_bits - static_cast<double>(key_bits) - len_bits;
        fp_len = prefix_bits / static_cast<double>(num_entries) + 1.7;
    }
    return fp_len;
}

template <typename t_itr>
inline size_t EstimateFingerprintLength(const t_itr begin, const t_itr end,
                                        size_t optimal_ls, size_t max_range_size,
                                        size_t unique_prefix_count, double bpk) {
    return ClampFingerprintLength(EstimateFingerprintLengthRaw(begin, end, optimal_ls, max_range_size,
                                                               unique_prefix_count, bpk));
}

template <typename t_itr>
inline int SelectSuffixLengthBySamples(const t_itr begin, const t_itr end,
                                       const std::vector<std::pair<uint64_t, uint64_t>>& queries,
                                       size_t max_range_size, double bpk) {
    if (begin == end || queries.empty()) return -1;
    const size_t r = max_range_size;
    const size_t unique_prefix_count = CountUniquePrefixes(begin, end, r);
    const uint64_t kmin = *begin, kmax = *(end - 1);
    struct QInfo { uint64_t L, R, k1, k2; bool has1, has2; int absent; };
    std::vector<QInfo> info;
    info.reserve(queries.size());
    for (const auto& q : queries) {
        const uint64_t L = std::max(q.first, kmin), R = std::min(q.second, kmax);
        if (L > R) continue;
        auto it = std::lower_bound(begin, end, L);
        if (it != end && *it <= R) continue;
        QInfo x{L, R, 0, 0, it != begin, it != end, 0};
        if (x.has1) x.k1 = *(it - 1);
        if (x.has2) x.k2 = *it;
        const uint64_t bL = L >> r, bR = R >> r;
        int present = 0;
        for (uint64_t b = bL; b <= bR; ++b) {
            if ((x.has1 && (x.k1 >> r) == b) || (x.has2 && (x.k2 >> r) == b)) ++present;
        }
        x.absent = static_cast<int>(bR - bL + 1) - present;
        info.push_back(x);
    }
    if (info.empty()) return -1;
    double best_cost = std::numeric_limits<double>::infinity();
    int best_ls = -1;
    for (size_t ls = 0; ls <= r; ++ls) {

        const double fp = static_cast<double>(
            ClampFingerprintLength(EstimateFingerprintLengthRaw(begin, end, ls, r, unique_prefix_count, bpk)));
        const double eps = EstimatePrefixFPR(fp);
        const size_t p = r - ls;
        double cost = 0.0, truncation = 0.0;
        for (const auto& x : info) {
            double t = 0.0;
            if (ls < r) {
                if ((x.has1 && (x.k1 >> p) == (x.L >> p)) || (x.has2 && (x.k2 >> p) == (x.R >> p))) t = 1.0;
            }
            truncation += t;
            cost += std::min(1.0, t + x.absent * eps);
        }
        if (ls < r && truncation == 0.0) cost += 0.5;
        if (cost <= best_cost) { best_cost = cost; best_ls = static_cast<int>(ls); }
    }
    return best_ls;
}

inline double EstimateTheoreticalFingerprintLength(size_t optimal_ls, double bpk,
                                                   size_t num_keys,
                                                   size_t unique_prefix_count) {
    if (num_keys == 0 || unique_prefix_count == 0) return 0.0;

    constexpr double kLoadFactor = 0.95;
    constexpr double kBucketSize = 4.0;
    const double n = static_cast<double>(num_keys);
    const double np = static_cast<double>(unique_prefix_count);
    const double avg = EstimateAverageKeysPerPrefix(num_keys, unique_prefix_count);

    if (optimal_ls == 0) {
        return (kLoadFactor * n * bpk) / (np * kBucketSize);
    }
    if (UseBinaryEncoding(optimal_ls, avg)) {
        return (kLoadFactor * n * (bpk - static_cast<double>(optimal_ls)))
             / (np * kBucketSize);
    }
    return (kLoadFactor / kBucketSize)
         * ((n * bpk) / np - std::pow(2.0, static_cast<double>(optimal_ls)));
}

inline double EstimateFPRAtLs(double D, size_t ls, size_t max_range_size,
                              double bpk, size_t num_keys,
                              size_t unique_prefix_count) {
    const double fp_len =
        EstimateTheoreticalFingerprintLength(ls, bpk, num_keys, unique_prefix_count);
    return std::min(1.0, EstimatePrefixFPR(fp_len) + EstimateSuffixFPR(D, ls, max_range_size));
}

inline double EstimateBinaryMergeFitCoeff() {
    const double w = EstimateTheoryLoadRatioW();
    const double ln2 = std::log(2.0);
    return w * ln2 * ln2;
}

inline double EstimateBinaryMergeLsSlope() {
    const double w = EstimateTheoryLoadRatioW();
    return 30.0 / (w + 1.0);
}

inline double EstimateBinaryMergeEta() {
    const double w = EstimateTheoryLoadRatioW();
    return std::log(2.0) / (w * (w + 1.0));
}

template <typename t_itr>
inline SegmentModelEstimate EstimateSegmentModel(const t_itr begin, const t_itr end,
                                                 double D, size_t max_range_size,
                                                 double bpk) {
    SegmentModelEstimate estimate;
    const size_t num_keys = std::distance(begin, end);
    if (num_keys == 0) return estimate;

    estimate.unique_prefix_count = CountUniquePrefixes(begin, end, max_range_size);
    estimate.optimal_ls = CalculateOriginalHourglassLs(
        D, bpk, num_keys, estimate.unique_prefix_count, max_range_size);
    estimate.fp_len = EstimateFingerprintLength(
        begin, end, estimate.optimal_ls, max_range_size, estimate.unique_prefix_count, bpk);
    estimate.prefix_fpr = EstimatePrefixFPR(static_cast<double>(estimate.fp_len));
    estimate.suffix_fpr = EstimateSuffixFPR(D, estimate.optimal_ls, max_range_size);
    estimate.estimated_fpr = std::min(1.0, estimate.prefix_fpr + estimate.suffix_fpr);
    estimate.branch =
        ClassifySegmentBranch(num_keys, estimate.unique_prefix_count, estimate.optimal_ls);

    if (estimate.branch == SegmentBranch::kUncorrelated) {
        return estimate;
    }

    if (estimate.branch == SegmentBranch::kBinary) {
        estimate.curvature = EstimateBinaryMergeFitCoeff();
        estimate.ls_slope = EstimateBinaryMergeLsSlope();
        estimate.marginal_gain = EstimateBinaryMergeEta();
        return estimate;
    }

    const double center = std::max(
        EstimateFPRAtLs(D, estimate.optimal_ls, max_range_size,
                        bpk, num_keys, estimate.unique_prefix_count),
        1e-12);
    if (estimate.optimal_ls > 0 && estimate.optimal_ls < max_range_size) {
        const double left = EstimateFPRAtLs(D, estimate.optimal_ls - 1, max_range_size,
                                           bpk, num_keys, estimate.unique_prefix_count);
        const double right = EstimateFPRAtLs(D, estimate.optimal_ls + 1, max_range_size,
                                            bpk, num_keys, estimate.unique_prefix_count);
        estimate.curvature = std::max(0.0, (left - 2.0 * center + right) / center);
    }

    const double d_delta = 0.1;
    const double low_D = std::max(0.0, D - d_delta);
    const double high_D = std::min(1.0, D + d_delta);
    if (high_D > low_D) {
        const size_t low_ls = CalculateOriginalHourglassLs(
            low_D, bpk, num_keys, estimate.unique_prefix_count, max_range_size);
        const size_t high_ls = CalculateOriginalHourglassLs(
            high_D, bpk, num_keys, estimate.unique_prefix_count, max_range_size);
        estimate.ls_slope = std::max(
            0.0,
            (static_cast<double>(high_ls) - static_cast<double>(low_ls)) / (high_D - low_D));
    }

    const double b_delta = 1.0;
    const double low_bpk = std::max(1.0, bpk - b_delta);
    const double high_bpk = bpk + b_delta;
    if (high_bpk > low_bpk) {
        const size_t low_ls = CalculateOriginalHourglassLs(
            D, low_bpk, num_keys, estimate.unique_prefix_count, max_range_size);
        const size_t high_ls = CalculateOriginalHourglassLs(
            D, high_bpk, num_keys, estimate.unique_prefix_count, max_range_size);
        const double low_fpr = std::max(
            EstimateFPRAtLs(D, low_ls, max_range_size,
                            low_bpk, num_keys, estimate.unique_prefix_count),
            1e-12);
        const double high_fpr = std::max(
            EstimateFPRAtLs(D, high_ls, max_range_size,
                            high_bpk, num_keys, estimate.unique_prefix_count),
            1e-12);
        estimate.marginal_gain =
            std::max(0.0, (std::log(low_fpr) - std::log(high_fpr)) / (high_bpk - low_bpk));
    }
    return estimate;
}

inline double EstimateSegmentIndexBits() {
    return 64.0 + sizeof(void*) * 8.0 + 8.0;
}

inline double EstimateImplementationSegmentMetadataBits() {
    const size_t shared_lookup_table_bytes = sizeof(PermEncoding);
    const size_t hourglass_header_bytes =
        sizeof(Hourglass) > shared_lookup_table_bytes
            ? sizeof(Hourglass) - shared_lookup_table_bytes
            : sizeof(Hourglass);
    return (hourglass_header_bytes + sizeof(Segment) + sizeof(std::unique_ptr<Hourglass>)) * 8.0;
}

inline double EstimateLogicalSegmentMetadataBits(double metadata_bits_per_segment) {
    return metadata_bits_per_segment > 0.0
        ? metadata_bits_per_segment
        : EstimateImplementationSegmentMetadataBits();
}

inline size_t EstimateCuckooAllocatedSlots(size_t unique_prefix_count) {
    if (unique_prefix_count == 0) return 0;

    constexpr size_t kAssoc = 4;
    size_t num_buckets = GetNextPow2(std::max<size_t>(1, unique_prefix_count / kAssoc));
    const double load =
        static_cast<double>(unique_prefix_count) / static_cast<double>(num_buckets * kAssoc);
    if (load > 0.95) {
        num_buckets <<= 1;
    }
    return num_buckets * kAssoc;
}

inline double EstimateCuckooTableLoad(size_t unique_prefix_count) {
    const size_t allocated_slots = EstimateCuckooAllocatedSlots(unique_prefix_count);
    if (allocated_slots == 0) return 1.0;
    return static_cast<double>(unique_prefix_count) / static_cast<double>(allocated_slots);
}

inline double EstimateCuckooTableLoad(const Segment& segment) {
    return EstimateCuckooTableLoad(segment.unique_prefix_count);
}

template <typename t_itr>
inline void RefreshSegmentModel(Segment& segment,
                                const t_itr begin,
                                double local_corr,
                                size_t max_range_size,
                                double bpk) {
    auto seg_begin = begin + segment.start_idx;
    auto seg_end = begin + segment.end_idx;
    const auto estimate = EstimateSegmentModel(seg_begin, seg_end, local_corr, max_range_size, bpk);

    segment.local_corr = local_corr;
    segment.unique_prefix_count = estimate.unique_prefix_count;
    segment.optimal_ls = estimate.optimal_ls;
    segment.fp_len = estimate.fp_len;
    segment.prefix_fpr = estimate.prefix_fpr;
    segment.suffix_fpr = estimate.suffix_fpr;
    segment.optimal_fpr = estimate.estimated_fpr;
    segment.branch = estimate.branch;
    segment.curvature = estimate.curvature;
    segment.ls_slope = estimate.ls_slope;
    segment.marginal_gain = estimate.marginal_gain;
}

template <typename t_itr>
std::vector<Segment> BuildInitialSegments_StrictM(
    const t_itr begin, const t_itr end,
    const std::vector<std::pair<uint64_t, uint64_t>>& sample_queries,
    size_t max_range_size,
    double bpk,
    size_t block_size) {
    std::vector<Segment> segments;
    const size_t total_keys = std::distance(begin, end);
    if (total_keys == 0 || block_size == 0) return segments;

    for (size_t i = 0; i < total_keys; i += block_size) {
        Segment seg;
        seg.start_idx = i;
        seg.end_idx = std::min(i + block_size, total_keys);
        seg.min_key = *(begin + seg.start_idx);
        seg.max_key = *(begin + seg.end_idx - 1);
        segments.push_back(seg);
    }

    for (const auto& q : sample_queries) {
        auto it = std::upper_bound(
            segments.begin(), segments.end(), q.first,
            [](uint64_t value, const Segment& segment) { return value < segment.min_key; });
        if (it != segments.begin()) {
            std::prev(it)->local_queries.push_back(q);
        }
    }

    const double global_rate = sample_queries.empty() ? 0.0
        : static_cast<double>(CountCorrelatedQueries(begin, end, sample_queries, max_range_size))
          / static_cast<double>(sample_queries.size());
    for (auto& seg : segments) {
        auto seg_begin = begin + seg.start_idx;
        auto seg_end = begin + seg.end_idx;
        seg.query_count = seg.local_queries.size();
        const double local_corr = EstimateLocalCorrValueShrunk(
            seg_begin, seg_end, seg.local_queries, max_range_size, global_rate);
        RefreshSegmentModel(seg, begin, local_corr, max_range_size, bpk);
        seg.local_queries.clear();
        seg.local_queries.shrink_to_fit();
    }

    return segments;
}

template <typename t_itr>
inline std::vector<Segment> CapacityCalibrateSegments(
    const std::vector<Segment>& input_segments,
    const t_itr begin,
    double merge_tolerance,
    double metadata_bits_per_segment,
    double bpk,
    size_t max_range_size) {
    std::vector<Segment> calibrated = input_segments;
    if (calibrated.size() < 2) return calibrated;

    constexpr double kLoadGuard = 0.80;
    constexpr double kCorrGuard = 0.15;
    constexpr double kRelaxedMergeGamma = 2.0;
    constexpr size_t kMaxRepairPasses = 5;

    const double saved_metadata_bits =
        EstimateLogicalSegmentMetadataBits(metadata_bits_per_segment);
    size_t accepted_repairs = 0;

    auto build_merged_segment = [&](const Segment& left,
                                    const Segment& right,
                                    Segment& merged_segment,
                                    SegmentModelEstimate& estimate,
                                    double& variance_term,
                                    double& relaxed_rhs) {
        const double n1 = static_cast<double>(left.num_keys());
        const double n2 = static_cast<double>(right.num_keys());
        const double n = n1 + n2;
        if (n == 0.0) return false;

        const double Dmix = (n1 * left.local_corr + n2 * right.local_corr) / n;
        auto merged_begin = begin + left.start_idx;
        auto merged_end = begin + right.end_idx;
        estimate = EstimateSegmentModel(merged_begin, merged_end, Dmix, max_range_size, bpk);

        variance_term =
            (n1 * n2 / n) * (right.local_corr - left.local_corr)
          * (right.local_corr - left.local_corr);
        relaxed_rhs = 0.0;
        const bool same_config = left.optimal_ls == right.optimal_ls && estimate.optimal_ls == left.optimal_ls;

        if (variance_term > merge_tolerance && !same_config) {
            if (estimate.branch == SegmentBranch::kUncorrelated
                || estimate.curvature <= 0.0
                || estimate.ls_slope <= 0.0
                || estimate.marginal_gain <= 0.0) {
                return false;
            }
            relaxed_rhs =
                kRelaxedMergeGamma * (2.0 * saved_metadata_bits * estimate.marginal_gain)
                / (estimate.curvature * estimate.ls_slope * estimate.ls_slope);
            if (variance_term > relaxed_rhs + merge_tolerance) {
                return false;
            }
        }

        merged_segment = left;
        merged_segment.end_idx = right.end_idx;
        merged_segment.max_key = right.max_key;
        merged_segment.query_count = left.query_count + right.query_count;
        merged_segment.local_corr = Dmix;
        merged_segment.unique_prefix_count = estimate.unique_prefix_count;
        merged_segment.optimal_ls = estimate.optimal_ls;
        merged_segment.fp_len = estimate.fp_len;
        merged_segment.prefix_fpr = estimate.prefix_fpr;
        merged_segment.suffix_fpr = estimate.suffix_fpr;
        merged_segment.optimal_fpr = estimate.estimated_fpr;
        merged_segment.branch = estimate.branch;
        merged_segment.curvature = estimate.curvature;
        merged_segment.ls_slope = estimate.ls_slope;
        merged_segment.marginal_gain = estimate.marginal_gain;
        return true;
    };

    for (size_t pass = 0; pass < kMaxRepairPasses && calibrated.size() >= 2; ++pass) {
        bool changed = false;

        for (size_t i = 0; i < calibrated.size(); ++i) {
            if (EstimateCuckooTableLoad(calibrated[i]) >= kLoadGuard) {
                continue;
            }

            bool has_best = false;
            size_t best_left_idx = 0;
            size_t best_right_idx = 0;
            Segment best_segment;
            double best_load = -1.0;
            double best_corr_gap = std::numeric_limits<double>::infinity();

            auto consider_neighbor = [&](size_t left_idx, size_t right_idx) {
                const Segment& left = calibrated[left_idx];
                const Segment& right = calibrated[right_idx];
                const double corr_gap = std::abs(left.local_corr - right.local_corr);
                if (corr_gap > kCorrGuard) {
                    return;
                }

                Segment candidate;
                SegmentModelEstimate estimate;
                double variance_term = 0.0;
                double relaxed_rhs = 0.0;
                if (!build_merged_segment(left, right, candidate, estimate,
                                          variance_term, relaxed_rhs)) {
                    return;
                }

                const double candidate_load = EstimateCuckooTableLoad(candidate);
                if (!has_best
                    || candidate_load > best_load + 1e-12
                    || (std::abs(candidate_load - best_load) <= 1e-12
                        && corr_gap < best_corr_gap)) {
                    has_best = true;
                    best_left_idx = left_idx;
                    best_right_idx = right_idx;
                    best_segment = candidate;
                    best_load = candidate_load;
                    best_corr_gap = corr_gap;
                }
            };

            if (i > 0) {
                consider_neighbor(i - 1, i);
            }
            if (i + 1 < calibrated.size()) {
                consider_neighbor(i, i + 1);
            }

            if (!has_best) {
                continue;
            }

            calibrated[best_left_idx] = best_segment;
            calibrated.erase(calibrated.begin() + static_cast<std::ptrdiff_t>(best_right_idx));
            ++accepted_repairs;
            changed = true;
            break;
        }

        if (!changed) {
            break;
        }
    }

    return calibrated;
}

template <typename t_itr>
inline std::vector<Segment> TheoryAwareMergeSegments(
    const std::vector<Segment>& initial_segments,
    const t_itr begin,
    double tau,
    double metadata_bits_per_segment,
    double bpk,
    size_t max_range_size,
    size_t base_block_size,
    double rho,
    double sample_rate) {
    (void)base_block_size;
    (void)rho;
    (void)sample_rate;

    std::vector<Segment> merged;
    if (initial_segments.empty()) return merged;

    if (tau < 0.0) {
        merged = initial_segments;
        return merged;
    }

    const double merge_tolerance = std::max(0.0, tau);
    const double saved_metadata_bits =
        EstimateLogicalSegmentMetadataBits(metadata_bits_per_segment);
    size_t accepted_merges = 0;

    auto try_merge_pair = [&](const Segment& left, const Segment& right, Segment& merged_segment) {
        const double n1 = static_cast<double>(left.num_keys());
        const double n2 = static_cast<double>(right.num_keys());
        const double n = n1 + n2;
        if (n == 0.0) return false;

        const double Dmix = (n1 * left.local_corr + n2 * right.local_corr) / n;
        auto merged_begin = begin + left.start_idx;
        auto merged_end = begin + right.end_idx;
        const auto estimate = EstimateSegmentModel(
            merged_begin, merged_end, Dmix, max_range_size, bpk);

        const double variance_term =
            (n1 * n2 / n) * (right.local_corr - left.local_corr)
          * (right.local_corr - left.local_corr);

        const bool same_config = left.optimal_ls == right.optimal_ls && estimate.optimal_ls == left.optimal_ls;
        if (variance_term <= merge_tolerance || same_config) {
            merged_segment = left;
            merged_segment.end_idx = right.end_idx;
            merged_segment.max_key = right.max_key;
            merged_segment.query_count = left.query_count + right.query_count;
            merged_segment.local_corr = Dmix;
            merged_segment.unique_prefix_count = estimate.unique_prefix_count;
            merged_segment.optimal_ls = estimate.optimal_ls;
            merged_segment.fp_len = estimate.fp_len;
            merged_segment.prefix_fpr = estimate.prefix_fpr;
            merged_segment.suffix_fpr = estimate.suffix_fpr;
            merged_segment.optimal_fpr = estimate.estimated_fpr;
            merged_segment.branch = estimate.branch;
            merged_segment.curvature = estimate.curvature;
            merged_segment.ls_slope = estimate.ls_slope;
            merged_segment.marginal_gain = estimate.marginal_gain;
            return true;
        }

        if (estimate.branch == SegmentBranch::kUncorrelated
            || estimate.curvature <= 0.0
            || estimate.ls_slope <= 0.0
            || estimate.marginal_gain <= 0.0) {
            return false;
        }

        const double merge_rhs =
            (2.0 * saved_metadata_bits * estimate.marginal_gain)
            / (estimate.curvature * estimate.ls_slope * estimate.ls_slope);

        if (variance_term > merge_rhs + merge_tolerance) {
            return false;
        }

        merged_segment = left;
        merged_segment.end_idx = right.end_idx;
        merged_segment.max_key = right.max_key;
        merged_segment.query_count = left.query_count + right.query_count;
        merged_segment.local_corr = Dmix;
        merged_segment.unique_prefix_count = estimate.unique_prefix_count;
        merged_segment.optimal_ls = estimate.optimal_ls;
        merged_segment.fp_len = estimate.fp_len;
        merged_segment.prefix_fpr = estimate.prefix_fpr;
        merged_segment.suffix_fpr = estimate.suffix_fpr;
        merged_segment.optimal_fpr = estimate.estimated_fpr;
        merged_segment.branch = estimate.branch;
        merged_segment.curvature = estimate.curvature;
        merged_segment.ls_slope = estimate.ls_slope;
        merged_segment.marginal_gain = estimate.marginal_gain;
        return true;
    };

    for (const auto& seg : initial_segments) {
        merged.push_back(seg);
        while (merged.size() >= 2) {
            Segment merged_segment;
            const size_t right_idx = merged.size() - 1;
            const size_t left_idx = right_idx - 1;
            if (!try_merge_pair(merged[left_idx], merged[right_idx], merged_segment)) {
                break;
            }
            merged.pop_back();
            merged.back() = merged_segment;
            ++accepted_merges;
        }
    }

    return CapacityCalibrateSegments(merged, begin, merge_tolerance,
                                     metadata_bits_per_segment, bpk, max_range_size);
}

class PiecewiseHourglass {
public:
    std::vector<Segment> routing_segments_;
    std::vector<std::unique_ptr<Hourglass>> instances_;
    long long routing_index_build_ns_ = 0;
    long long local_filters_build_ns_ = 0;

    std::vector<std::pair<uint64_t, uint64_t>> sorted_samples_;

    template <typename t_itr>
    PiecewiseHourglass(const std::vector<Segment>& macro_segments,
                       const t_itr begin,
                       size_t max_range_size,
                       double bpk,
                       size_t base_block_size,
                       const std::vector<std::pair<uint64_t, uint64_t>>& sample_queries = {}) {
        sorted_samples_ = sample_queries;
        std::sort(sorted_samples_.begin(), sorted_samples_.end());
        if (base_block_size == 0) {
            throw std::runtime_error("base_block_size must be positive");
        }
        for (const auto& segment : macro_segments) {
            if (segment.num_keys() == 0) continue;
            CreateBinaryInstances(begin, segment, max_range_size, bpk, base_block_size);
        }
    }

    uint64_t SizeInBits() const {
        uint64_t total_bits = sizeof(*this) * 8
                            + routing_segments_.size() * sizeof(Segment) * 8
                            + instances_.size() * sizeof(std::unique_ptr<Hourglass>) * 8;
        for (const auto& instance : instances_) {
            total_bits += SafeHourglassSizeInBits(*instance);
        }
        return total_bits;
    }

    uint64_t CompactSizeInBits() const {
        uint64_t total_bits = sizeof(*this) * 8
                            + static_cast<uint64_t>(routing_segments_.size() * EstimateSegmentIndexBits())
                            + instances_.size() * sizeof(std::unique_ptr<Hourglass>) * 8;
        for (const auto& instance : instances_) {
            total_bits += SafeHourglassSizeInBits(*instance);
        }
        return total_bits;
    }

    uint64_t PayloadSizeInBits() const {
        uint64_t total_bits = 0;
        for (const auto& instance : instances_) {
            total_bits += SafeHourglassPayloadSizeInBits(*instance);
        }
        return total_bits;
    }

    bool Query(uint64_t left, uint64_t right) {
        auto it = std::upper_bound(
            routing_segments_.begin(), routing_segments_.end(), left,
            [](uint64_t value, const Segment& segment) { return value < segment.min_key; });

        size_t idx = 0;
        if (it != routing_segments_.begin()) {
            idx = static_cast<size_t>(std::distance(routing_segments_.begin(), std::prev(it)));
        }

        while (idx < routing_segments_.size() && routing_segments_[idx].min_key <= right) {
            const auto& segment = routing_segments_[idx];
            if (left <= segment.max_key) {
                const uint64_t clipped_left = std::max(left, segment.min_key);
                const uint64_t clipped_right = std::min(right, segment.max_key);
                if (QueryLocalInstance(*instances_[idx], clipped_left, clipped_right)) {
                    return true;
                }
            }
            ++idx;
        }
        return false;
    }

private:
    static size_t LargestPowerOfTwoMultiple(size_t remaining, size_t base_block_size) {
        if (remaining <= base_block_size) return remaining;
        const size_t blocks = remaining / base_block_size;
        const size_t power = 1ULL << static_cast<size_t>(std::floor(std::log2(blocks)));
        return power * base_block_size;
    }

    static bool QueryLocalInstance(Hourglass& hourglass, uint64_t query_left, uint64_t query_right) {
        if (!hourglass.correlated_) {
            return hourglass.Lookup_Adaptive(query_left, query_right);
        }

        const size_t prefix_bits = hourglass.num_layers_ + hourglass.pruned_layers_;
        if (prefix_bits == 0 || prefix_bits >= 64) {
            return hourglass.Lookup_Adaptive(query_left, query_right);
        }

        uint64_t current_left = query_left;
        while (current_left <= query_right) {
            const uint64_t prefix = current_left >> prefix_bits;
            const uint64_t prefix_right = ((prefix + 1ULL) << prefix_bits) - 1ULL;
            const uint64_t current_right = std::min(query_right, prefix_right);
            if (hourglass.Lookup_Adaptive(current_left, current_right)) {
                return true;
            }
            if (current_right == std::numeric_limits<uint64_t>::max()) {
                break;
            }
            current_left = current_right + 1ULL;
        }
        return false;
    }

    template <typename t_itr>
    void CreateBinaryInstances(const t_itr begin, const Segment& logical_segment,
                               size_t max_range_size, double bpk,
                               size_t base_block_size) {
        size_t cursor = logical_segment.start_idx;
        while (cursor < logical_segment.end_idx) {
            const size_t remaining = logical_segment.end_idx - cursor;
            const size_t chunk_keys = LargestPowerOfTwoMultiple(remaining, base_block_size);

            Segment physical_segment = logical_segment;
            physical_segment.start_idx = cursor;
            physical_segment.end_idx = cursor + chunk_keys;
            physical_segment.query_count = static_cast<size_t>(
                std::llround(static_cast<double>(logical_segment.query_count)
                             * static_cast<double>(chunk_keys)
                             / static_cast<double>(logical_segment.num_keys())));
            CreateInstance(begin, physical_segment, max_range_size, bpk);

            cursor += chunk_keys;
        }
    }

    template <typename t_itr>
    void CreateInstance(const t_itr begin, const Segment& logical_segment,
                        size_t max_range_size, double bpk) {
        Segment segment = logical_segment;
        segment.min_key = *(begin + segment.start_idx);
        segment.max_key = *(begin + segment.end_idx - 1);

        auto seg_begin = begin + segment.start_idx;
        auto seg_end = begin + segment.end_idx;
        const auto estimate =
            EstimateSegmentModel(seg_begin, seg_end, segment.local_corr, max_range_size, bpk);
        segment.unique_prefix_count = estimate.unique_prefix_count;
        segment.optimal_ls = estimate.optimal_ls;
        segment.fp_len = estimate.fp_len;
        segment.prefix_fpr = estimate.prefix_fpr;
        segment.suffix_fpr = estimate.suffix_fpr;
        segment.optimal_fpr = estimate.estimated_fpr;
        segment.branch = estimate.branch;
        segment.curvature = estimate.curvature;
        segment.ls_slope = estimate.ls_slope;
        segment.marginal_gain = estimate.marginal_gain;

        const auto index_start = std::chrono::high_resolution_clock::now();
        routing_segments_.push_back(segment);
        const auto index_end = std::chrono::high_resolution_clock::now();
        routing_index_build_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
            index_end - index_start).count();

        const auto filter_start = std::chrono::high_resolution_clock::now();

        int forced_ls = -1;
        if (!sorted_samples_.empty()) {
            const uint64_t span = 1ULL << max_range_size;
            const uint64_t lo = segment.min_key >= span ? segment.min_key - span : 0;
            auto first = std::lower_bound(sorted_samples_.begin(), sorted_samples_.end(),
                                          std::make_pair(lo, uint64_t{0}));
            std::vector<std::pair<uint64_t, uint64_t>> local;
            for (auto it = first; it != sorted_samples_.end() && it->first <= segment.max_key; ++it) {
                if (it->second >= segment.min_key) local.push_back(*it);
            }
            forced_ls = SelectSuffixLengthBySamples(seg_begin, seg_end, local, max_range_size, bpk);
        }

        constexpr int kMaxSameConfigRetries = 5;
        std::unique_ptr<Hourglass> hourglass;
        float corr = static_cast<float>(segment.local_corr);
        if (forced_ls > 0 && corr <= 0.0f) corr = 0.01f;
        int failures = 0;
        for (;;) {
            try {
                hourglass = std::make_unique<Hourglass>(
                    seg_begin, seg_end, max_range_size, corr, static_cast<size_t>(std::llround(bpk)), forced_ls);
                break;
            } catch (const std::runtime_error& e) {
                ++failures;
                if (failures < kMaxSameConfigRetries) {
                    continue;
                }
                failures = 0;
                if (forced_ls > 0) {
                    --forced_ls;
                } else {
                    if (corr <= 0.0f) throw;
                    corr = std::max(0.0f, corr - 0.05f);
                }
            }
        }
        routing_segments_.back().local_corr = corr;
        instances_.push_back(std::move(hourglass));
        const auto filter_end = std::chrono::high_resolution_clock::now();
        local_filters_build_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
            filter_end - filter_start).count();
    }
};

class HourglassPlus {
public:
    template <typename t_itr>
    HourglassPlus(const t_itr begin, const t_itr end, size_t max_range_size, double bpk,
                  const std::vector<std::pair<uint64_t, uint64_t>>& sample_queries,
                  size_t base_m = 12, double sample_rate = 0.15) {
        const size_t num_keys = static_cast<size_t>(end - begin);
        const size_t nominal_block_size = 4ULL * (1ULL << base_m);
        const size_t block_size = std::max<size_t>(
            1, static_cast<size_t>(std::floor(0.95 * static_cast<double>(nominal_block_size))));
        const double segment_bits = EstimateLogicalSegmentMetadataBits(0.0);
        auto initial = BuildInitialSegments_StrictM(begin, end, sample_queries, max_range_size, bpk, block_size);
        auto merged = TheoryAwareMergeSegments(initial, begin, 0.0, segment_bits, bpk,
                                               max_range_size, block_size, 0.0, sample_rate);
        num_logical_segments_ = merged.size();
        const double metadata_bpk = num_keys == 0 ? 0.0
            : static_cast<double>(std::llround(static_cast<double>(merged.size()) * segment_bits)) / num_keys;
        const double payload_bpk = std::max(1.0, bpk - metadata_bpk);
        filter_ = std::make_unique<PiecewiseHourglass>(merged, begin, max_range_size, payload_bpk,
                                                       block_size, sample_queries);
    }

    bool Query(uint64_t left, uint64_t right) { return filter_->Query(left, right); }

    uint64_t SizeInBits() const {
        return filter_->PayloadSizeInBits()
             + static_cast<uint64_t>(std::llround(filter_->instances_.size() * EstimateImplementationSegmentMetadataBits()));
    }

    size_t NumLogicalSegments() const { return num_logical_segments_; }

    size_t NumInstances() const { return filter_->instances_.size(); }

private:
    std::unique_ptr<PiecewiseHourglass> filter_;
    size_t num_logical_segments_ = 0;
};

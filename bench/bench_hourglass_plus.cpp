#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../hourglass_plus.hpp"

template <typename T>
std::vector<T> load_binary(const std::string& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + file);
    T n = 0;
    in.read(reinterpret_cast<char*>(&n), sizeof(T));
    std::vector<T> v(static_cast<size_t>(n));
    in.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(T));
    return v;
}

template <typename F>
void Evaluate(const std::string& name, double input_bpk, size_t n_keys, double build_ms,
              uint64_t size_bits, size_t segments, F&& query,
              const std::vector<std::pair<uint64_t, uint64_t>>& bounds, const std::vector<int>& truth) {
    size_t negatives = 0, fp = 0, fn = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < bounds.size(); ++i) {
        bool r = query(bounds[i].first, bounds[i].second);
        if (truth[i] == 0) {
            ++negatives;
            fp += r;
        } else if (!r) {
            ++fn;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    double query_ns = std::chrono::duration<double, std::nano>(end - start).count() / bounds.size();
    std::cout << name << ',' << input_bpk << ',' << static_cast<double>(size_bits) / n_keys << ','
              << static_cast<double>(fp) / negatives << ',' << fn << ',' << build_ms << ','
              << query_ns << ',' << segments << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 7) {
        std::cerr << "usage: " << argv[0]
                  << " <keys> <left> <right> <result> <bpk> <range_size> [base_m=12] [sample_rate=0.15] [seed=20260422]\n";
        return 1;
    }
    const double bpk = std::stod(argv[5]);
    const size_t range_size = std::stoull(argv[6]);
    const size_t base_m = argc > 7 ? std::stoull(argv[7]) : 12;
    const double sample_rate = argc > 8 ? std::stod(argv[8]) : 0.15;
    const uint32_t seed = argc > 9 ? static_cast<uint32_t>(std::stoul(argv[9])) : 20260422u;

    auto keys = load_binary<uint64_t>(argv[1]);
    auto left = load_binary<uint64_t>(argv[2]);
    auto right = load_binary<uint64_t>(argv[3]);
    auto truth = load_binary<int>(argv[4]);
    std::sort(keys.begin(), keys.end());

    std::vector<std::pair<uint64_t, uint64_t>> bounds;
    for (size_t i = 0; i < left.size(); ++i) bounds.emplace_back(left[i], right[i]);
    std::vector<std::pair<uint64_t, uint64_t>> samples;
    std::sample(bounds.begin(), bounds.end(), std::back_inserter(samples),
                static_cast<size_t>(bounds.size() * sample_rate), std::mt19937{seed});

    std::cout << "method,input_bpk,bpk,fpr,false_negatives,build_ms,query_ns,segments" << std::endl;
    {
        auto start = std::chrono::high_resolution_clock::now();
        const double corr = EstimateLocalCorrValue(keys.begin(), keys.end(), samples, range_size);
        Hourglass h(keys.begin(), keys.end(), range_size, corr, static_cast<size_t>(std::llround(bpk)));
        auto end = std::chrono::high_resolution_clock::now();
        double build_ms = std::chrono::duration<double, std::milli>(end - start).count();
        Evaluate("hourglass", bpk, keys.size(), build_ms, SafeHourglassPayloadSizeInBits(h) + static_cast<uint64_t>(EstimateImplementationSegmentMetadataBits()), 1,
                 [&](uint64_t l, uint64_t r) { return h.Lookup_Adaptive(l, r); }, bounds, truth);
    }
    {
        auto start = std::chrono::high_resolution_clock::now();
        HourglassPlus h(keys.begin(), keys.end(), range_size, bpk, samples, base_m, sample_rate);
        auto end = std::chrono::high_resolution_clock::now();
        double build_ms = std::chrono::duration<double, std::milli>(end - start).count();
        Evaluate("hourglass_plus", bpk, keys.size(), build_ms, h.SizeInBits(), h.NumLogicalSegments(),
                 [&](uint64_t l, uint64_t r) { return h.Query(l, r); }, bounds, truth);
    }
    return 0;
}

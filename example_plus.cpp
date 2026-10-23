#include "hourglass_plus.hpp"
#include <chrono>
#include <random>
#include <set>

int main()
{
    std::mt19937_64 gen(42);
    std::uniform_int_distribution<uint64_t> distr_value(0, UINT64_MAX - 1);

    size_t n_keys = 15'938'355;
    size_t n_queries = 1'000'000;
    size_t range_size = 10;
    double bpk = 16;

    std::set<uint64_t> key_set;
    while (key_set.size() < n_keys) {
        key_set.insert(distr_value(gen));
    }
    std::vector<uint64_t> keys(key_set.begin(), key_set.end());

    std::vector<std::pair<uint64_t, uint64_t>> queries;
    std::uniform_int_distribution<size_t> distr_key(0, n_keys - 2);
    std::uniform_int_distribution<uint64_t> distr_gap(1, 1ULL << 16);
    while (queries.size() < n_queries) {
        uint64_t left;
        if (queries.size() < n_queries / 2) {
            left = distr_value(gen);
        } else {
            size_t i = distr_key(gen);
            left = keys[i] + distr_gap(gen);
        }
        uint64_t right = left + (1ULL << range_size) - 1;
        if (right < left) continue;
        auto it = key_set.lower_bound(left);
        if (it != key_set.end() && *it <= right) continue;
        queries.emplace_back(left, right);
    }

    std::vector<std::pair<uint64_t, uint64_t>> samples;
    std::sample(queries.begin(), queries.end(), std::back_inserter(samples),
                static_cast<size_t>(queries.size() * 0.15), std::mt19937{20260422});

    auto start = std::chrono::high_resolution_clock::now();
    HourglassPlus h(keys.begin(), keys.end(), range_size, bpk, samples);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "[+] Construction Time: " << duration.count() << " ms" << std::endl;
    std::cout << "[+] Segments: " << h.NumLogicalSegments() << std::endl;

    size_t fp = 0;
    start = std::chrono::high_resolution_clock::now();
    for (const auto& q : queries) {
        if (h.Query(q.first, q.second)) fp++;
    }
    end = std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "[+] Query Time: " << ns / (double)n_queries << " ns/query" << std::endl;
    std::cout << "[+] FPR: " << fp / (double)n_queries << std::endl;
    std::cout << "[+] BPK: " << h.SizeInBits() / (double)n_keys << std::endl;

    return 0;
}

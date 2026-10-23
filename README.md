# Hourglass+

Hourglass+ is an adaptive range filter that defends against recurring false positives in skewed and adversarial range queries.
It extends Hourglass (SIGMOD 2026) with correlation-aware segmentation, which partitions the key space into segments and selects fingerprint and suffix lengths for each segment based on its key-query correlation estimated from sampled queries.

## Tested settings

- CMake 3.28.3
- GCC 13.2.0
- OpenSSL 3.0.13

## Build

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j16
```

Adaptivity is enabled by default. Use `-DADAPT=OFF` to disable it.
If CMake cannot find OpenSSL, set `-DOPENSSL_ROOT_DIR=<your_openssl_path>`.

## Usage

`hourglass.cpp` implements one Hourglass instance and `hourglass_plus.hpp` implements Hourglass+.

```cpp
#include "hourglass_plus.hpp"

// keys: sorted std::vector<uint64_t>
// samples: sampled range queries, std::vector<std::pair<uint64_t, uint64_t>>
// range_size: log2 of the maximum range size, bpk: bits per key
HourglassPlus filter(keys.begin(), keys.end(), range_size, bpk, samples);
bool maybe_non_empty = filter.Query(left, right);
```

Run the examples from the `build` directory.

```
./hourglass
./hourglass_plus
```

## Benchmark

See [bench/benchmark.md](bench/benchmark.md).

## Third-party code

- `roaring.c` and `roaring.h` are from CRoaring (Apache License 2.0).
- `sux/` is from the Sux library (GNU LGPL 3.0 or later).
- `hashutil.cc`, `hashutil.h`, and `permencoding.h` are adapted from the cuckoo filter library (Apache License 2.0).
- `bench/bench_utils.hpp`, `bench/workload_gen.cpp`, and `bench/scripts/download_datasets.sh` are from Grafite (GNU GPL 3.0 or later).

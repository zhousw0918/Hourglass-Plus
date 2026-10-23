# Build
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DADAPT=OFF
make -j16
```

# Generate workloads
Download the datasets and generate keys with the Grafite workload generator.
```
bash ./bench/scripts/download_datasets.sh
./build/bench/workload_gen --binary-keys ./real_datasets/books_200M_uint64 -n 15938355 -q 1000000 --range-size 0 5 10
```

Generate regional workloads whose key-query correlation varies across 32 regions from 0 to 0.9.
```
python3 ./bench/scripts/generate_regional_workload.py ./books/keys ./workloads/books
```

# Run benchmark
```
./build/bench/bench_hourglass_plus ./workloads/books/r10/{keys,left,right,result} <bpk> 10
```
The output reports the space, FPR, construction time, and query time of Hourglass and Hourglass+.

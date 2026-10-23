include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG        v3.0
)
FetchContent_MakeAvailable(argparse)

add_executable(workload_gen workload_gen.cpp)
target_compile_features(workload_gen PRIVATE cxx_std_20)
target_link_libraries(workload_gen argparse)

add_executable(bench_hourglass_plus bench_hourglass_plus.cpp)
target_link_libraries(bench_hourglass_plus hourglasslib)

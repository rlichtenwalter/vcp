# Counts-container study. Built only when both `VCP_BUILD_BENCHMARKS=ON`
# AND `VCP_BUILD_COUNTS_BENCH=ON` are set. The study is a one-time
# investigation into replacing the std::map return value of
# vcp<4, r, d>::generate_vector; its artifacts are not shipped.

option(VCP_BUILD_COUNTS_BENCH "Build the counts-container study harness" OFF)

if(VCP_BUILD_COUNTS_BENCH)
    add_executable(bench_counts_containers bench_counts_containers.cpp)
    target_link_libraries(bench_counts_containers
        PRIVATE vcp::vcp Catch2::Catch2WithMain)
    target_compile_definitions(bench_counts_containers PRIVATE NDEBUG)
    target_compile_options(bench_counts_containers PRIVATE
        -O3 -Wall -Wextra -pedantic)
endif()

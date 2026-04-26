# Small-map study microbenchmark. Built only when both
# `VCP_BUILD_BENCHMARKS=ON` AND `VCP_BUILD_SMALL_MAP_STUDY=ON` are set.
# The study is a one-time investigation that motivated
# `detail::dense_or_sparse_map`; its artifacts are not shipped with the
# library. Recorded results live under `results/`.

option(VCP_BUILD_SMALL_MAP_STUDY "Build the small-map study microbenchmark" OFF)

if(VCP_BUILD_SMALL_MAP_STUDY)
    add_executable(bench_container_mix bench_container_mix.cpp)
    target_link_libraries(bench_container_mix PRIVATE vcp::vcp Catch2::Catch2WithMain)
    target_compile_definitions(bench_container_mix PRIVATE NDEBUG)
    target_compile_options(bench_container_mix PRIVATE -O3 -Wall -Wextra -pedantic)
endif()

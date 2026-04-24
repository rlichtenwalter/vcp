# Small-map study instrumentation harness. Built only when both
# `VCP_BUILD_BENCHMARKS=ON` AND `VCP_BUILD_K_PROBE=ON` are set. The study
# is a one-time investigation that motivated `detail::dense_or_sparse_map`;
# its artifacts are not shipped with the library.

option(VCP_BUILD_K_PROBE "Build the k-probe study harness" OFF)

if(VCP_BUILD_K_PROBE)
    add_executable(k_probe_tool k_probe_tool.cpp)
    target_link_libraries(k_probe_tool PRIVATE vcp::vcp CLI11::CLI11)
    target_compile_definitions(k_probe_tool PRIVATE
        VCP_INSTRUMENT_K
        VCP_VERSION="${PROJECT_VERSION}"
        NDEBUG
    )
    target_include_directories(k_probe_tool PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_options(k_probe_tool PRIVATE -O3 -Wall -Wextra -pedantic)

    add_executable(bench_container_mix bench_container_mix.cpp)
    target_link_libraries(bench_container_mix PRIVATE vcp::vcp Catch2::Catch2WithMain)
    target_compile_definitions(bench_container_mix PRIVATE NDEBUG)
    target_compile_options(bench_container_mix PRIVATE -O3 -Wall -Wextra -pedantic)
endif()

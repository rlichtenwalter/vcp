# Benchmark target. Built only when VCP_BUILD_BENCHMARKS is ON.
#
# Reuses the Catch2 instance fetched by test/CMakeLists.txt when both are
# enabled; otherwise fetches its own. The benchmark executable is NOT
# registered with ctest (matches the kdtree/mRMR pattern) so it does not
# run during `make test`. Invoke it explicitly:
#
#   ./build/benchmark/bench_vcp "[!benchmark]"
#
# See benchmark/README.md.

include(FetchContent)

# Fetch Catch2 only if the test subdirectory did not already fetch it.
if(NOT TARGET Catch2::Catch2WithMain)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.13.0
    )
    FetchContent_MakeAvailable(Catch2)
endif()

add_executable(bench_vcp bench_vcp.cpp)
target_link_libraries(bench_vcp PRIVATE vcp::vcp Catch2::Catch2WithMain)

# Always optimize the benchmark binary, regardless of the parent build
# type - benchmarks only make sense at -O3 -DNDEBUG. Use per-target
# definitions/options rather than abusing -DNDEBUG as a compile flag, and
# rely on the LINK_OPTIONS form so the Release linker flags (notably
# -flto) do not get silently dropped when the parent is Debug.
target_compile_definitions(bench_vcp PRIVATE NDEBUG)
target_compile_options(bench_vcp PRIVATE
    -O3
    -Wall -Wextra -pedantic -Wno-unused-local-typedefs
)
# Propagate VCP_SANITIZE flags if enabled at the root. Sanitized
# benchmarks are mostly a smoke signal - the -O3 force here fights with
# -O0 defaults sanitizers want - but running the micro-benchmarks under
# ASan+UBSan at least proves the library headers don't UB under the
# hot-path instantiations the CLI tests do not cover.
target_compile_options(bench_vcp PRIVATE ${VCP_SANITIZE_FLAGS})
target_link_options(bench_vcp    PRIVATE ${VCP_SANITIZE_FLAGS})

# Optional one-off instrumentation study. Opt-in via VCP_BUILD_K_PROBE.
add_subdirectory(small_map_study)

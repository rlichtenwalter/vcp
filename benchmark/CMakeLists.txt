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
target_compile_options(bench_vcp PRIVATE
    -O3 -DNDEBUG
    -Wall -Wextra -pedantic -Wno-unused-local-typedefs
)

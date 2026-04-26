cmake_minimum_required(VERSION 3.24)

# Use Boost's own CMake config (BoostConfig.cmake) instead of the removed
# FindBoost module. Requires Boost >= 1.70 (Ubuntu 20.04+, Rocky 9+).
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VCP_VERSION_STRING)
string(STRIP "${VCP_VERSION_STRING}" VCP_VERSION_STRING)

project(vcp
    VERSION ${VCP_VERSION_STRING}
    LANGUAGES CXX
    DESCRIPTION "Header-only C++ library for Vertex Collocation Profile computation on graphs"
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
include(FetchContent)

# --- Dependencies ---
# Boost (headers-only; uses boost/multiprecision/cpp_int.hpp).
find_package(Boost REQUIRED)

# CLI11 (modern header-only command-line parser) fetched at configure time.
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.2
)
FetchContent_MakeAvailable(CLI11)

# --- Header-only library ---
add_library(vcp INTERFACE)
add_library(vcp::vcp ALIAS vcp)

target_compile_features(vcp INTERFACE cxx_std_20)

target_include_directories(vcp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(vcp INTERFACE Boost::headers)

set(VCP_PUBLIC_HEADERS
    include/vcp/directed_graph.hpp
    include/vcp/graph.hpp
    include/vcp/multirelational_directed_graph.hpp
    include/vcp/multirelational_graph.hpp
    include/vcp/square_matrix.hpp
    include/vcp/vcp.hpp
    include/vcp/vcp_3_1_0.hpp
    include/vcp/vcp_3_1_1.hpp
    include/vcp/vcp_3_r_0.hpp
    include/vcp/vcp_3_r_1.hpp
    include/vcp/vcp_4_1_0.hpp
    include/vcp/vcp_4_1_1.hpp
    include/vcp/vcp_4_r_0.hpp
    include/vcp/vcp_4_r_1.hpp
    include/vcp/vcp_dynamic_mapper.hpp
    include/vcp/vcp_static_mapper.hpp
    include/vcp/detail/dense_or_sparse_map.hpp
    include/vcp/detail/v3_buffer_bound.hpp
)

target_sources(vcp INTERFACE
    FILE_SET HEADERS
    BASE_DIRS include
    FILES ${VCP_PUBLIC_HEADERS}
)

# --- CLI tools ---
set(VCP_TOOLS
    vcp_generate
    vcp_map
    directed_to_undirected
    ell_2_pairs
    vcp_gen_er_fixture
)

# --- Warning flags ---
# Shared baseline applied to every target that consumes the library
# (tools, tests, benchmarks). Centralized so that a flag added here lands
# in every consumer build at once and the test/benchmark targets cannot
# silently slide past warnings the tools would error on.
#
# -Wno-unused-local-typedefs suppresses noise from CRTP/SFINAE-style
# helpers in the public headers. Removing it surfaces hundreds of
# warnings across template instantiations that have no actionable fix.
#
# Beyond -Wall -Wextra -pedantic, the expanded set catches:
#   -Wconversion           lossy implicit numeric conversions
#   -Wsign-conversion      signed<->unsigned implicit conversions
#   -Wshadow               variables shadowing outer-scope names
#   -Wnull-dereference     statically-provable null derefs
#   -Wdouble-promotion     implicit float->double (float literal to
#                          double via "3.14" instead of "3.14f", etc.)
#   -Wimplicit-fallthrough switch case missing [[fallthrough]] attr
# GCC adds two more that Clang does not implement:
#   -Wlogical-op           `a || a` / `a && a` typos
#   -Wduplicated-cond      `if(a) ... else if(a)` typos
set(VCP_WARNING_FLAGS
    -Wall -Wextra -Werror -pedantic -Wno-unused-local-typedefs
    -Wconversion -Wsign-conversion -Wshadow -Wnull-dereference
    -Wdouble-promotion -Wimplicit-fallthrough
    $<$<CXX_COMPILER_ID:GNU>:-Wlogical-op>
    $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-cond>
)

# Clang's -flto produces LLVM bitcode that the default GNU `ld` cannot
# read ("file format not recognized"). Route Clang link steps through
# `lld`, which speaks bitcode natively. CI installs `lld` alongside
# `clang`. GCC builds keep their default linker.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_link_options(-fuse-ld=lld)
endif()

foreach(tool IN LISTS VCP_TOOLS)
    add_executable(${tool} tools/${tool}.cpp)
    target_link_libraries(${tool} PRIVATE vcp::vcp CLI11::CLI11)
    target_compile_definitions(${tool} PRIVATE
        VCP_VERSION="${PROJECT_VERSION}"
    )
    target_compile_options(${tool} PRIVATE
        ${VCP_WARNING_FLAGS}
        $<$<CONFIG:Release>:-O3 -fomit-frame-pointer -DNDEBUG>
        $<$<CONFIG:Debug>:-Og -g -fno-omit-frame-pointer>
    )
    # Enable LTO via the CMake property rather than a raw -flto flag.
    # CMake handles the gcc-vs-clang toolchain difference (gcc-ar /
    # gcc-ranlib for gcc; llvm-ar / llvm-ranlib for clang) automatically,
    # selecting -flto / -flto=thin / -flto=auto as appropriate per
    # compiler.
    set_target_properties(${tool} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE
    )
endforeach()

# --- Sanitizers ---
# When -DVCP_SANITIZE=ON, every consumer of the library (tools, tests,
# benchmarks) is compiled and linked with AddressSanitizer +
# UndefinedBehaviorSanitizer. Guarded on Debug builds so Release artifacts
# are unaffected.
#
# -fno-sanitize-recover=all turns every sanitizer diagnostic into a hard
# error (abort + non-zero exit). Without this, UBSan reports and
# continues, which makes it too easy for CI to miss regressions.
#
# The flags are applied per-target rather than via the header-only
# INTERFACE library to avoid leaking build-time sanitizer requirements
# into the installed/exported targets that downstream consumers would
# pick up via find_package(vcp).
option(VCP_SANITIZE "Enable AddressSanitizer + UndefinedBehaviorSanitizer in Debug builds" OFF)
# Each flag is its own generator expression. Bundling them into a single
# genex with semicolon separators happens to expand correctly on Ninja +
# GNU Make with GCC today, but the behavior depends on how each generator
# splits the genex output into argv entries - a multi-element list form
# is specified and portable.
set(VCP_SANITIZE_FLAGS
    $<$<BOOL:${VCP_SANITIZE}>:$<$<CONFIG:Debug>:-fsanitize=address,undefined>>
    $<$<BOOL:${VCP_SANITIZE}>:$<$<CONFIG:Debug>:-fno-sanitize-recover=all>>
    $<$<BOOL:${VCP_SANITIZE}>:$<$<CONFIG:Debug>:-fno-omit-frame-pointer>>
)
foreach(tool IN LISTS VCP_TOOLS)
    target_compile_options(${tool} PRIVATE ${VCP_SANITIZE_FLAGS})
    target_link_options(${tool}    PRIVATE ${VCP_SANITIZE_FLAGS})
endforeach()

# --- Tests ---
option(VCP_BUILD_TESTS "Build tests" ${PROJECT_IS_TOP_LEVEL})
if(VCP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# --- Benchmarks ---
option(VCP_BUILD_BENCHMARKS "Build library-level micro-benchmarks" OFF)
if(VCP_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()

# --- Install ---
include(CMakePackageConfigHelpers)

install(TARGETS vcp
    EXPORT vcpTargets
    FILE_SET HEADERS
)

install(TARGETS ${VCP_TOOLS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT vcpTargets
    FILE vcpTargets.cmake
    NAMESPACE vcp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vcp
)

configure_package_config_file(
    cmake/vcpConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/vcpConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vcp
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/vcpConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/vcpConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/vcpConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vcp
)

configure_file(cmake/vcp.pc.in "${CMAKE_CURRENT_BINARY_DIR}/vcp.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/vcp.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

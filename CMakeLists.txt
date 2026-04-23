cmake_minimum_required(VERSION 3.21)

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

# --- User-configurable parameters ---
set(VCP_MAX_NEIGHBORS 16384 CACHE STRING
    "Maximum number of neighbors per vertex (compile-time upper bound for buffer allocation)")

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

target_compile_features(vcp INTERFACE cxx_std_14)

target_include_directories(vcp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(vcp INTERFACE Boost::headers)

target_compile_definitions(vcp INTERFACE MAX_NEIGHBORS=${VCP_MAX_NEIGHBORS})

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
)

foreach(tool IN LISTS VCP_TOOLS)
    add_executable(${tool} tools/${tool}.cpp)
    target_link_libraries(${tool} PRIVATE vcp::vcp CLI11::CLI11)
    target_compile_definitions(${tool} PRIVATE
        VCP_VERSION="${PROJECT_VERSION}"
    )
    target_compile_options(${tool} PRIVATE
        -Wall -Wextra -Werror -pedantic -Wno-unused-local-typedefs
        $<$<CONFIG:Release>:-O3 -flto -fomit-frame-pointer -DNDEBUG>
        $<$<CONFIG:Debug>:-Og -g -fno-omit-frame-pointer>
    )
endforeach()

option(VCP_SANITIZE "Enable AddressSanitizer in Debug builds" OFF)
if(VCP_SANITIZE)
    foreach(tool IN LISTS VCP_TOOLS)
        target_compile_options(${tool} PRIVATE $<$<CONFIG:Debug>:-fsanitize=address>)
        target_link_options(${tool} PRIVATE $<$<CONFIG:Debug>:-fsanitize=address>)
    endforeach()
endif()

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

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.13.0
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# --- Unit tests ---
# Each test executable is registered for CTest via catch_discover_tests.
# Sanitizer flags propagate from the root so tests exercise the same
# instrumentation as tools when VCP_SANITIZE=ON.

function(vcp_add_catch_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE vcp::vcp Catch2::Catch2WithMain)
    target_compile_options(${name} PRIVATE ${VCP_SANITIZE_FLAGS})
    target_link_options(${name}    PRIVATE ${VCP_SANITIZE_FLAGS})
    catch_discover_tests(${name})
endfunction()

vcp_add_catch_test(test_placeholder)
vcp_add_catch_test(test_square_matrix)
vcp_add_catch_test(test_graph)
vcp_add_catch_test(test_directed_graph)
vcp_add_catch_test(test_multirelational_graph)
vcp_add_catch_test(test_multirelational_directed_graph)
vcp_add_catch_test(test_vcp_dynamic_mapper)
vcp_add_catch_test(test_vcp_static_mapper)
vcp_add_catch_test(test_dense_or_sparse_map)
vcp_add_catch_test(test_vcp_3_1_0)
vcp_add_catch_test(test_vcp_3_1_1)
vcp_add_catch_test(test_vcp_4_1_0)
vcp_add_catch_test(test_vcp_4_1_1)
vcp_add_catch_test(test_vcp_4_r_0)
vcp_add_catch_test(test_vcp_4_r_1)

# --- CLI smoke tests ---
# Each tool responds to --version and --help. These verify CLI11 wiring and
# that each binary boots without segfault.
foreach(tool IN LISTS VCP_TOOLS)
    set(cli $<TARGET_FILE:${tool}>)

    add_test(NAME ${tool}_version COMMAND ${cli} --version)
    set_tests_properties(${tool}_version PROPERTIES
        PASS_REGULAR_EXPRESSION "v[0-9]+\\.[0-9]+\\.[0-9]+")

    add_test(NAME ${tool}_help COMMAND ${cli} --help)
    set_tests_properties(${tool}_help PROPERTIES
        PASS_REGULAR_EXPRESSION "[Uu]sage")
endforeach()

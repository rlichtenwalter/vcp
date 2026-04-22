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
# Placeholder for the tooling-modernization pass; real unit tests land in a
# follow-up release.
add_executable(test_placeholder test_placeholder.cpp)
target_link_libraries(test_placeholder PRIVATE vcp::vcp Catch2::Catch2WithMain)
catch_discover_tests(test_placeholder)

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

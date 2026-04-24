// Unit tests for vcp<3, 1, true> — the directed n=3, r=1 specialization.
//
// Returns std::array<unsigned long, 64>. Each slot encodes a 6-bit
// (V1V2, V1V3, V2V3) tuple where each pair uses 2 bits (OUT=1, IN=2,
// BOTH=3, unconnected=0). Sum invariant: sum(counts) == V - 2.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <numeric>
#include <sstream>
#include <string>

#include <vcp/directed_graph.hpp>
#include <vcp/vcp_3_1_1.hpp>

namespace {

vcp::directed_graph parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::directed_graph g;
  iss >> g;
  return g;
}

vcp::const_vertex_iterator vertex_at(vcp::directed_graph const &g, std::size_t id) {
  return g.vertices_begin() + id;
}

// Enum tags reproduced locally (private in the header). Per-pair stride
// is 4 in 2-bit encoding.
constexpr std::size_t OUT = 1;
constexpr std::size_t IN = 2;
constexpr std::size_t BOTH = 3;
constexpr std::size_t V1V2 = 1;
constexpr std::size_t V1V3 = 4;
constexpr std::size_t V2V3 = 16;

} // namespace

TEST_CASE("vcp<3,1,true>::element_count is 64", "[vcp_3_1_1]") {
  REQUIRE(vcp::vcp<3, 1, true>::element_count() == 64);
}

TEST_CASE("vcp<3,1,true> directed triangle routes to OUT*V1V2 + IN*V1V3 + OUT*V2V3 bucket",
          "[vcp_3_1_1]") {
  // 0 -> 1, 1 -> 2, 2 -> 0. From pivot (v1=0, v2=1):
  //   v1v2 = 0->1 = OUT (no reverse)
  //   v1v3 (v3=2) = 2->0 so v1's perspective is IN
  //   v2v3 (v3=2) = 1->2 so v2's perspective is OUT
  vcp::directed_graph g = parse("1\n2\n0\n");
  vcp::vcp<3, 1, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == g.vertex_count() - 2);
  // Expected bucket: OUT*V1V2 + IN*V1V3 + OUT*V2V3.
  std::size_t const expected = OUT * V1V2 + IN * V1V3 + OUT * V2V3;
  REQUIRE(counts[expected] == 1);
}

TEST_CASE("vcp<3,1,true> mutual triangle routes to BOTH,BOTH,BOTH bucket", "[vcp_3_1_1]") {
  // All three pairs mutual.
  vcp::directed_graph g = parse("1 2\n0 2\n0 1\n");
  vcp::vcp<3, 1, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  std::size_t const expected = BOTH * V1V2 + BOTH * V1V3 + BOTH * V2V3;
  REQUIRE(counts[expected] == 1);
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == 1);
}

TEST_CASE("vcp<3,1,true> sum invariant across larger graph", "[vcp_3_1_1][invariant]") {
  // Mixed directed 5-vertex graph.
  vcp::directed_graph g = parse("1 2\n0\n3\n1\n\n");
  REQUIRE(g.vertex_count() == 5);
  vcp::vcp<3, 1, true> v(g);
  for (std::size_t a(0); a < 5; ++a) {
    for (std::size_t b(0); b < 5; ++b) {
      if (a != b) {
        auto const counts = v.generate_vector(vertex_at(g, a), vertex_at(g, b));
        unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
        REQUIRE(sum == g.vertex_count() - 2);
      }
    }
  }
}

TEST_CASE("vcp<3,1,true> empty directed graph counts the fully-unconnected bucket", "[vcp_3_1_1]") {
  vcp::directed_graph g = parse("\n\n\n");
  vcp::vcp<3, 1, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  REQUIRE(counts[0] == 1);
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == 1);
}

// Unit tests for vcp<3, 1, false> — the undirected n=3, r=1 specialization.
//
// This specialization returns a std::array<unsigned long, 8> where the
// index is the direct (V1V2 + V1V3 + V2V3) bit-packed connectivity
// address:
//   bit 0 = V1V2 (edge between the pivot pair)
//   bit 1 = V1V3 (edge v1 -> v3 candidate)
//   bit 2 = V2V3 (edge v2 -> v3 candidate)
//
// There are 8 elements (one per 3-bit subset). For V vertices, the sum
// invariant is sum(counts) == V - 2 (each of the V-2 non-pivot vertices
// contributes exactly one bucket).

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <numeric>
#include <sstream>
#include <string>

#include <vcp/graph.hpp>
#include <vcp/vcp_3_1_0.hpp>

namespace {

vcp::graph parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::graph g;
  iss >> g;
  return g;
}

vcp::const_vertex_iterator vertex_at(vcp::graph const &g, std::size_t id) {
  return g.vertices_begin() + id;
}

// Enum tags reproduced locally from vcp_3_1_0.hpp (private in the header).
constexpr std::size_t V1V2 = 1;
constexpr std::size_t V1V3 = 2;
constexpr std::size_t V2V3 = 4;

} // namespace

TEST_CASE("vcp<3,1,false>::element_count is 8", "[vcp_3_1_0]") {
  REQUIRE(vcp::vcp<3, 1, false>::element_count() == 8);
}

TEST_CASE("vcp<3,1,false> triangle at pair (0,1) counts exactly one triangle bucket",
          "[vcp_3_1_0]") {
  // Triangle 0-1, 1-2, 0-2. V1V2=1, V1V3=1, V2V3=1 → bucket 7.
  vcp::graph g = parse("1 2\n0 2\n0 1\n");
  vcp::vcp<3, 1, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  REQUIRE(counts[V1V2 + V1V3 + V2V3] == 1);
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == g.vertex_count() - 2);
}

TEST_CASE("vcp<3,1,false> P3 path at non-adjacent endpoints routes to open-triad bucket",
          "[vcp_3_1_0]") {
  // Path 0-1-2. At pair (v1=0, v2=2): v1v2 absent, v3=vertex 1 has edges
  // to both → V1V3+V2V3 bucket. V=3, sum = 1.
  vcp::graph g = parse("1\n0 2\n1\n");
  vcp::vcp<3, 1, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 2));
  REQUIRE(counts[V1V3 + V2V3] == 1);
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == 1);
}

TEST_CASE("vcp<3,1,false> sum invariant across larger graph", "[vcp_3_1_0][invariant]") {
  // 5-cycle. Any pair sees 3 non-pivot vertices.
  vcp::graph g = parse("1 4\n0 2\n1 3\n2 4\n0 3\n");
  REQUIRE(g.vertex_count() == 5);
  vcp::vcp<3, 1, false> v(g);
  for (std::size_t a(0); a < 5; ++a) {
    for (std::size_t b(a + 1); b < 5; ++b) {
      auto const counts = v.generate_vector(vertex_at(g, a), vertex_at(g, b));
      unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
      REQUIRE(sum == g.vertex_count() - 2);
    }
  }
}

TEST_CASE("vcp<3,1,false> isolated-only graph counts empty bucket", "[vcp_3_1_0]") {
  // 3 isolated vertices. V1V2=0, V1V3=0, V2V3=0 → bucket 0.
  vcp::graph g = parse("\n\n\n");
  vcp::vcp<3, 1, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  REQUIRE(counts[0] == 1);
  unsigned long const sum = std::accumulate(counts.begin(), counts.end(), 0UL);
  REQUIRE(sum == 1);
}

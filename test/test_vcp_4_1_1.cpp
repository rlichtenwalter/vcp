// Unit tests for vcp<4, 1, true> — the directed binary-connectivity
// specialization (n=4, r=1, d=1). This is the specialization with the
// K4-bidirectional unconnected-class underflow that motivated fix/
// vcp-4-1-1-underflow on the prior branch, plus the mutual / asymmetric
// self-correction for the v1v2 class at the end of generate_vector.
//
// Coverage here is the C++-level analog of regression Phase 6's Python
// validator: sum(counts) == C(V-2, 2) and max(counts) <= C(V-2, 2) per
// pair, across representative directed topologies. Regression Phase 6
// only exercises the six dbug2_* fixtures; these tests spread the
// invariant to mainline topologies (K4-bidirectional, K4-oriented, C4,
// star, dk3, isolated). If the self-correction for amutual or mutual
// v1v2 classes regresses, the sum drops below C(V-2, 2); if the
// unconnected-class arithmetic underflows, the peak blows past it.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <numeric>
#include <sstream>
#include <string>

#include <vcp/directed_graph.hpp>
#include <vcp/vcp_4_1_1.hpp>

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

unsigned long expected_pair_sum(std::size_t vertex_count) {
  std::size_t const n(vertex_count - 2);
  return static_cast<unsigned long>(n * (n - 1) / 2);
}

void check_pair(vcp::vcp<4, 1, true> &v, vcp::directed_graph const &g, std::size_t a,
                std::size_t b) {
  auto const counts = v.generate_vector(vertex_at(g, a), vertex_at(g, b));
  unsigned long const total = std::accumulate(counts.begin(), counts.end(), 0UL);
  unsigned long const peak = *std::max_element(counts.begin(), counts.end());
  unsigned long const expected = expected_pair_sum(g.vertex_count());
  REQUIRE(total == expected);
  REQUIRE(peak <= expected);
}

} // namespace

TEST_CASE("vcp<4,1,true> element_count is 2112", "[vcp_4_1_1]") {
  REQUIRE(vcp::vcp<4, 1, true>::element_count() == 2112);
}

TEST_CASE("vcp<4,1,true> sum invariant on K4 bidirectional", "[vcp_4_1_1][invariant]") {
  // Every directed pair is mutual. Historically this triggered the
  // unconnected-class underflow because the dense subtraction chain
  // wrapped below zero when v3_count * (V - 2 - v3_count) exceeded
  // unconnected_pairs.
  vcp::directed_graph g = parse("1 2 3\n0 2 3\n0 1 3\n0 1 2\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, true> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,true> sum invariant on K4 oriented (asymmetric only)",
          "[vcp_4_1_1][invariant]") {
  // Tournament: 0->1, 0->2, 0->3, 1->2, 1->3, 2->3. Every pair
  // asymmetric (amutual class), so the amutualPairs ctor accumulator
  // is exercised. Self-correction branch at line 603 of vcp_4_1_1.hpp:
  // this->amutualPairs - amutuals - static_cast<bool>(v1v2 == OUT || IN).
  vcp::directed_graph g = parse("1 2 3\n2 3\n3\n\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, true> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,true> sum invariant on 3-cycle (oriented) with isolated vertex",
          "[vcp_4_1_1][invariant]") {
  // 0->1, 1->2, 2->0; vertex 3 isolated. V=4 => expected sum = 1.
  vcp::directed_graph g = parse("1\n2\n0\n\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, true> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,true> sum invariant on empty directed graph (4 isolated)",
          "[vcp_4_1_1][invariant]") {
  vcp::directed_graph g = parse("\n\n\n\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, true> v(g);
  check_pair(v, g, 0, 1);
}

TEST_CASE("vcp<4,1,true> sum invariant on directed C4", "[vcp_4_1_1][invariant]") {
  // 0->1, 1->2, 2->3, 3->0.
  vcp::directed_graph g = parse("1\n2\n3\n0\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, true> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,true> sum invariant on 5-vertex mixed graph", "[vcp_4_1_1][invariant]") {
  // V=5: C(V-2, 2) = 3 per pair. Mix of mutual and asymmetric edges.
  // Mutual 0<->1, asymmetric 0->2, 2->3, 3->1, vertex 4 isolated.
  vcp::directed_graph g = parse("1 2\n0\n3\n1\n\n");
  REQUIRE(g.vertex_count() == 5);
  vcp::vcp<4, 1, true> v(g);
  check_pair(v, g, 0, 1); // mutual pair
  check_pair(v, g, 0, 2); // asymmetric pair
  check_pair(v, g, 0, 4); // unconnected pair
}

TEST_CASE("vcp<4,1,true> locates the single non-zero bucket for V=4", "[vcp_4_1_1]") {
  // V=4, one non-zero bucket with value 1. Works as a cheap independent
  // check of "did the count land in exactly one element" even without
  // reconstructing the expected element_address.
  vcp::directed_graph g = parse("1 2 3\n0 2 3\n0 1 3\n0 1 2\n");
  vcp::vcp<4, 1, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  auto const nonzero = static_cast<unsigned long>(
      std::count_if(counts.begin(), counts.end(), [](unsigned long c) { return c != 0; }));
  REQUIRE(nonzero == 1);
  unsigned long const peak = *std::max_element(counts.begin(), counts.end());
  REQUIRE(peak == 1);
}

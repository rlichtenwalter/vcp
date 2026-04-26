// Unit tests for vcp<4, 1, false> — the undirected binary-connectivity
// specialization (n=4, r=1, d=0).
//
// This specialization had the K4 unconnected-class underflow hazard; the
// NDEBUG assertions are live only in debug builds, so we cover the same
// invariants at release-capable unit-test level via sum + peak checks.
//
// Ground-truth approach: for any vertex pair (v1, v2) and V-vertex input,
// the sum of counts across all 40 VCP elements must equal C(V-2, 2) (the
// number of unordered {v3, v4} pairs in the remaining V-2 vertices).
// Peak over the 40 slots must not exceed C(V-2, 2) — an explicit
// diagnostic for unsigned underflow (would produce a ~2**64 value).
//
// These invariants are "absolute": they hold regardless of topology, so
// they extend cleanly to richer fixtures without the Python regression
// harness's dependency on legacy-build binaries.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <numeric>
#include <sstream>
#include <string>

#include <vcp/graph.hpp>
#include <vcp/vcp_4_1_0.hpp>

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

unsigned long expected_pair_sum(std::size_t vertex_count) {
  // C(V-2, 2) = (V-2)(V-3)/2, the number of unordered {v3, v4} pairs.
  std::size_t const n(vertex_count - 2);
  return static_cast<unsigned long>(n * (n - 1) / 2);
}

void check_pair(vcp::vcp<4, 1, false> &v, vcp::graph const &g, std::size_t a, std::size_t b) {
  auto const counts = v.generate_vector(vertex_at(g, a), vertex_at(g, b));
  unsigned long const total = std::accumulate(counts.begin(), counts.end(), 0UL);
  unsigned long const peak = *std::ranges::max_element(counts);
  unsigned long const expected = expected_pair_sum(g.vertex_count());
  REQUIRE(total == expected);
  REQUIRE(peak <= expected);
}

} // namespace

TEST_CASE("vcp<4,1,false> element_count is 40", "[vcp_4_1_0]") {
  REQUIRE(vcp::vcp<4, 1, false>::element_count() == 40);
}

TEST_CASE("vcp<4,1,false> sum invariant on K4 (all pairs)", "[vcp_4_1_0][invariant]") {
  vcp::graph g = parse("1 2 3\n0 2 3\n0 1 3\n0 1 2\n");
  REQUIRE(g.vertex_count() == 4);
  REQUIRE(g.edge_count() == 6);
  vcp::vcp<4, 1, false> v(g);
  // Every ordered pair; V=4 means C(V-2, 2) = 1, so each call returns a
  // vector summing to 1. K4 is the maximally-symmetric undirected case
  // and a common regression trigger.
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,false> sum invariant on P4 path", "[vcp_4_1_0][invariant]") {
  // 0-1-2-3 path.
  vcp::graph g = parse("1\n0 2\n1 3\n2\n");
  REQUIRE(g.vertex_count() == 4);
  REQUIRE(g.edge_count() == 3);
  vcp::vcp<4, 1, false> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,false> sum invariant on star4 (center=0)", "[vcp_4_1_0][invariant]") {
  // Vertex 0 connected to 1, 2, 3; others isolated besides 0.
  vcp::graph g = parse("1 2 3\n0\n0\n0\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, false> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,false> sum invariant on C4 cycle", "[vcp_4_1_0][invariant]") {
  // 0-1-2-3-0 cycle.
  vcp::graph g = parse("1 3\n0 2\n1 3\n0 2\n");
  REQUIRE(g.vertex_count() == 4);
  vcp::vcp<4, 1, false> v(g);
  for (std::size_t a(0); a < 4; ++a) {
    for (std::size_t b(0); b < 4; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,false> sum invariant on empty graph (4 isolated vertices)",
          "[vcp_4_1_0][invariant]") {
  // Pathological case: no edges at all. C(V-2, 2) = 1; the non-zero
  // bucket should be the fully-unconnected one. Catches underflow in the
  // unconnected_pairs subtraction chain (most-negative case — gaps = 0,
  // v3_count = 0, v4_count = 0, unconnected_pairs = 6).
  vcp::graph g = parse("\n\n\n\n");
  REQUIRE(g.vertex_count() == 4);
  REQUIRE(g.edge_count() == 0);
  vcp::vcp<4, 1, false> v(g);
  check_pair(v, g, 0, 1);
}

TEST_CASE("vcp<4,1,false> sum invariant on 5-vertex graphs", "[vcp_4_1_0][invariant]") {
  // V=5: C(V-2, 2) = 3 per pair. Pentagon cycle. Exercise every pair —
  // the invariant is per-pair and a topology-dependent bug could pass
  // on adjacent pairs while failing on non-adjacent ones.
  vcp::graph g = parse("1 4\n0 2\n1 3\n2 4\n0 3\n");
  REQUIRE(g.vertex_count() == 5);
  vcp::vcp<4, 1, false> v(g);
  for (std::size_t a(0); a < 5; ++a) {
    for (std::size_t b(0); b < 5; ++b) {
      if (a != b) {
        check_pair(v, g, a, b);
      }
    }
  }
}

TEST_CASE("vcp<4,1,false> single-bucket concentration on K5", "[vcp_4_1_0]") {
  // In K5 every non-pivot pair {v3, v4} is fully connected both to each
  // other and to the pivot pair. All C(V-2, 2) = 3 counts therefore
  // land in the same bucket — a strictly stronger invariant than sum
  // alone, and cheap.
  vcp::graph g = parse("1 2 3 4\n0 2 3 4\n0 1 3 4\n0 1 2 4\n0 1 2 3\n");
  REQUIRE(g.vertex_count() == 5);
  vcp::vcp<4, 1, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  auto const nonzero = static_cast<unsigned long>(
      std::count_if(counts.begin(), counts.end(), [](unsigned long c) { return c != 0; }));
  REQUIRE(nonzero == 1);
  unsigned long const peak = *std::ranges::max_element(counts);
  REQUIRE(peak == 3);
}

TEST_CASE("vcp<4,1,false> locates the single non-zero bucket for V=4", "[vcp_4_1_0]") {
  // V=4 means C(V-2, 2) = 1, so exactly one count slot is 1 and the
  // other 39 are zero. This is the finest-grained invariant we can
  // assert without re-deriving the element_address lookup table.
  vcp::graph g = parse("1 2 3\n0 2 3\n0 1 3\n0 1 2\n");
  vcp::vcp<4, 1, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));
  auto const nonzero = static_cast<unsigned long>(
      std::count_if(counts.begin(), counts.end(), [](unsigned long c) { return c != 0; }));
  REQUIRE(nonzero == 1);
  unsigned long const peak = *std::ranges::max_element(counts);
  REQUIRE(peak == 1);
}

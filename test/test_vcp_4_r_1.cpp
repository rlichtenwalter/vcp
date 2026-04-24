// Regression tests for vcp<4, r, true> (multirelational directed).
//
// Two bugs under investigation:
//
//   1. vcp_4_r_1.hpp:337 writes v4's in-direction edge value into
//      it1->second(3, 3) — the diagonal — instead of (3, 2). The
//      diagonal is not part of the canonical subgraph address, so the
//      value is silently discarded and the address is computed with the
//      stale (3, 2) cell (typically 0) instead of the actual v4→v3
//      edge value. Every v3-v4 "it2 role promotion with edge" case in
//      the `it1 < it2` branch is affected.
//
//   2. The post-enumeration accounting loop applies the
//      "−1 iff v1v2 is this directionality class" self-correction only
//      to the (0, 0) unconnected class (lines 391–394). The analogous
//      correction in vcp_4_1_1 is present for all three classes
//      (asymmetric / mutual / unconnected) — lines 587–590 of that file.
//      If edge_types counts v1v2, then non-unconnected classes lose a
//      count of 1 when v1v2 belongs to that class. These tests verify
//      empirically whether this hypothesis produces wrong output.

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>

#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_4_r_1.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace {

vcp::multirelational_directed_graph<2> parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::multirelational_directed_graph<2> g;
  iss >> g;
  return g;
}

vcp::const_vertex_iterator vertex_at(vcp::multirelational_directed_graph<2> const &g,
                                     std::size_t id) {
  return g.vertices_begin() + id;
}

using conn_t = typename vcp::multirelational_directed_graph<2>::connectivity_address_type;
using matrix_t = vcp::square_matrix<conn_t, 4>;
using address_t = typename vcp::vcp_dynamic_mapper<4, 2, true>::subgraph_address_type;

// Directed matrix has 12 off-diagonal cells; the diagonal is never part
// of the canonical address. Arguments follow (row, col): out-direction
// from row to col. `ij` = edge i→j value.
address_t canonical_of(conn_t c01, conn_t c10, conn_t c02, conn_t c20, conn_t c03, conn_t c30,
                       conn_t c12, conn_t c21, conn_t c13, conn_t c31, conn_t c23, conn_t c32) {
  matrix_t m;
  m(0, 1) = c01;
  m(1, 0) = c10;
  m(0, 2) = c02;
  m(2, 0) = c20;
  m(0, 3) = c03;
  m(3, 0) = c30;
  m(1, 2) = c12;
  m(2, 1) = c21;
  m(1, 3) = c13;
  m(3, 1) = c31;
  m(2, 3) = c23;
  m(3, 2) = c32;
  vcp::vcp_dynamic_mapper<4, 2, true> mapper;
  return mapper.canonical_subgraph_address(m);
}

} // namespace

TEST_CASE("vcp<4,2,true> it2 role promotion with v3-v4 edge writes (3,2), not (3,3)",
          "[vcp_4_r_1][bug][ground-truth]") {
  // Four-vertex directed multirelational graph:
  //   0 → 2 rel 1        (v1 out to v3 candidate 2)
  //   0 → 3 rel 1        (v1 out to v3 candidate 3)
  //   2 → 3 rel 1        (v3 = 2 out to v4 = 3)
  //   3 → 2 rel 2        (v4 out back to v3 — populates (3, 2))
  // v1 = 0, v2 = 1. Vertex 1 has no edges.
  //
  // The algorithm enters the it1/it2 loop with it1 = v3=2, it2 = v3=3,
  // hits the "v3-v4 edge + role promotion" branch at line 325, and
  // should write (2, 3) = 1 and (3, 2) = 2. The bug writes (3, 3) = 2
  // (into the unused diagonal) and leaves (3, 2) = 0, routing the
  // enumerated subgraph to the wrong canonical bucket.
  auto const g = parse("2,1 3,1\n"
                       "\n"
                       "3,1\n"
                       "2,2\n");
  REQUIRE(g.vertex_count() == 4);

  vcp::vcp<4, 2, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));

  // Expected induced subgraph on (v1=0, v2=1, v3=2, v4=3):
  //   0→1 = 0, 1→0 = 0, 0→2 = 1, 2→0 = 0, 0→3 = 1, 3→0 = 0,
  //   1→2 = 0, 2→1 = 0, 1→3 = 0, 3→1 = 0, 2→3 = 1, 3→2 = 2
  address_t const expected = canonical_of(0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 2);
  auto it = counts.find(expected);
  REQUIRE(it != counts.end());
  REQUIRE(it->second == 1);
}

// Investigation of the suspected missing self-correction (task #68).
// Construct a graph where v1v2 is mutual (directionality class BOTH,BOTH)
// and the graph has additional mutual pairs elsewhere. If the correction
// is missing, the mutual-class bucket will be inflated by 1.
//
// Sum invariant: sum(counts) must equal C(V-2, 2) = 1 for V=4, one pair.
// A missing −1 correction would break the sum invariant by leaving it at 2.
TEST_CASE("vcp<4,2,true> post-enum self-correction preserves sum invariant "
          "when v1v2 is non-unconnected",
          "[vcp_4_r_1][investigate]") {
  // Four-vertex graph. v1=0, v2=1 mutual via rel 1:
  //   0 → 1 rel 1
  //   1 → 0 rel 1
  // Additional isolated structure so no v3 candidates exist, letting the
  // post-enum loop drive the result alone.
  //   Vertex 2, 3 isolated.
  auto const g = parse("1,1\n"
                       "0,1\n"
                       "\n"
                       "\n");
  REQUIRE(g.vertex_count() == 4);

  vcp::vcp<4, 2, true> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));

  unsigned long total = 0;
  for (auto const &kv : counts) {
    total += kv.second;
  }
  // C(4-2, 2) = 1. If the mutual-class self-correction is missing, this
  // reports 2 instead of 1.
  REQUIRE(total == 1);
}

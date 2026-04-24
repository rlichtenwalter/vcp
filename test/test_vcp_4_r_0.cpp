// Regression tests for vcp<4, r, false> (multirelational undirected).
//
// The v1/v2 neighbor merge in vcp_4_r_0.hpp::generate_vector had a
// copy-paste defect: the v2-exclusive branch (line 109) and the v2-drain
// tail loop (line 143) wrote v3Vertices(1,2) — the (v2,v3) edge value —
// using edge_value(v1_neighbors_it) instead of edge_value(v2_neighbors_it).
// That corrupts the connectivity matrix for every v3 candidate that is a
// neighbor of v2 but not of v1, which in turn routes the enumerated
// subgraph to a wrong canonical VCP bucket.
//
// Why uncaught before: the existing Phase 3 regression is a byte-for-byte
// legacy-vs-current diff. Both builds carry the same bug since 2012, so
// the diff reports PASS. These tests compute the expected canonical
// bucket via the same dynamic mapper the algorithm uses and assert the
// enumerated subgraph lands there — a ground-truth check independent of
// any reference build.

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>

#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_4_r_0.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace {

vcp::multirelational_graph<2> parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::multirelational_graph<2> g;
  iss >> g;
  return g;
}

vcp::const_vertex_iterator vertex_at(vcp::multirelational_graph<2> const &g, std::size_t id) {
  return g.vertices_begin() + id;
}

using conn_t = typename vcp::multirelational_graph<2>::connectivity_address_type;
using matrix_t = vcp::square_matrix<conn_t, 4>;
using address_t = typename vcp::vcp_dynamic_mapper<4, 2, false>::subgraph_address_type;

// Convenience: build the (undirected) connectivity matrix from the six
// unordered pair values, then ask the mapper for its canonical address.
address_t canonical_of(conn_t v01, conn_t v02, conn_t v03, conn_t v12, conn_t v13, conn_t v23) {
  matrix_t m;
  m(0, 1) = v01;
  m(0, 2) = v02;
  m(0, 3) = v03;
  m(1, 2) = v12;
  m(1, 3) = v13;
  m(2, 3) = v23;
  vcp::vcp_dynamic_mapper<4, 2, false> mapper;
  return mapper.canonical_subgraph_address(m);
}

} // namespace

TEST_CASE("vcp<4,2,false> merge v2-exclusive branch routes to correct bucket",
          "[vcp_4_r_0][bug][ground-truth]") {
  // Four-vertex graph with edges
  //   (0,1) relation 1
  //   (0,3) relation 1
  //   (1,2) relation 2
  // At pair (v1=0, v2=1) the merge loop reaches the v2-exclusive branch
  // with v3 candidate = vertex 2. The correct (1,2) slot value is 2
  // (the relation-2 edge between v2 and v3). The bug wrote the value of
  // edge (0, 3) = 1 instead. The induced subgraph on {0,1,2,3} has:
  //   (0,1)=1, (0,2)=0, (0,3)=1, (1,2)=2, (1,3)=0, (2,3)=0
  // which must land in its canonical bucket with count == 1.
  auto const g = parse("1,1 3,1\n"
                       "0,1 2,2\n"
                       "1,2\n"
                       "0,1\n");
  REQUIRE(g.vertex_count() == 4);

  vcp::vcp<4, 2, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));

  address_t const expected = canonical_of(1, 0, 1, 2, 0, 0);
  auto it = counts.find(expected);
  REQUIRE(it != counts.end());
  REQUIRE(it->second == 1);
}

TEST_CASE("vcp<4,2,false> merge v2-drain tail loop routes to correct bucket",
          "[vcp_4_r_0][bug][ground-truth]") {
  // Same shape but arranged so the bug trips in the v2-drain tail loop
  // (line 143) instead of the main merge's v2-exclusive branch (line 109).
  // After the merge loop the v1 iterator is at end(); the buggy code
  // dereferenced it via edge_value(), yielding an unrelated edge's value
  // or 0 (via the sentinel guard in multirelational_graph::edge_value).
  //   (0,1) relation 1     — v1 connected only to v2
  //   (1,2) relation 2     — v2 connected to an exclusive v3
  //   vertex 3 isolated
  // Expected matrix: (0,1)=1, (1,2)=2, everything else 0.
  auto const g = parse("1,1\n"
                       "0,1 2,2\n"
                       "1,2\n"
                       "\n");
  REQUIRE(g.vertex_count() == 4);

  vcp::vcp<4, 2, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));

  address_t const expected = canonical_of(1, 0, 0, 2, 0, 0);
  auto it = counts.find(expected);
  REQUIRE(it != counts.end());
  REQUIRE(it->second == 1);
}

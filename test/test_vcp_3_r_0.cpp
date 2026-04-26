// Regression tests for vcp<3, r, false> (multirelational undirected, n=3).
//
// vcp_3_r_0.hpp::generate_vector had the same family of copy-paste defect
// fixed earlier in vcp_4_r_0.hpp: the V2V3 branch of the main merge loop
// (line 117) wrote the (v2, v3) edge bitmask using edge_value(v1_it)
// instead of edge_value(v2_it). For r == 1 the edge value is always 1 so
// the bug is invisible; for r >= 2 it routes the enumerated subgraph to
// the wrong canonical bucket.
//
// Why uncaught before: the byte-for-byte legacy regression diff carries
// the same bug on both sides. These tests compute the expected canonical
// bucket via the dynamic mapper and assert the enumerated subgraph lands
// there — independent of any reference build.

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>

#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_3_r_0.hpp>
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
using matrix_t = vcp::square_matrix<conn_t, 3>;
using address_t = typename vcp::vcp_dynamic_mapper<3, 2, false>::subgraph_address_type;

// Build the (undirected) 3-vertex connectivity matrix from the three
// unordered-pair bitmasks and return its canonical address.
address_t canonical_of(conn_t v01, conn_t v02, conn_t v12) {
  matrix_t m;
  m(0, 1) = v01;
  m(0, 2) = v02;
  m(1, 2) = v12;
  vcp::vcp_dynamic_mapper<3, 2, false> mapper;
  return mapper.canonical_subgraph_address(m);
}

} // namespace

TEST_CASE("vcp<3,2,false> merge V2V3 branch routes to correct bucket",
          "[vcp_3_r_0][bug][ground-truth]") {
  // Four-vertex graph with edges
  //   (0,1) bitmask 3 — pivot edge (sets V1V2)
  //   (0,3) bitmask 1 — v1-only neighbor at target=3
  //   (1,2) bitmask 2 — v2-only neighbor at target=2
  //
  // At pair (v1=0, v2=1) the merge loop reaches the V2V3 branch
  // (target_of(v1_it)=3 > target_of(v2_it)=2) for v3=2. Correct V2V3
  // value is edge_value(v2_it) = 2; the bug substituted
  // edge_value(v1_it) = 1, routing the subgraph to the wrong bucket.
  // Distinct edge bitmasks (3, 1, 2) are required: at r=1 every edge
  // value is 1 and the bug is masked.
  auto const g = parse("1,3 3,1\n"
                       "0,3 2,2\n"
                       "1,2\n"
                       "0,1\n");
  REQUIRE(g.vertex_count() == 4);

  vcp::vcp<3, 2, false> v(g);
  auto const counts = v.generate_vector(vertex_at(g, 0), vertex_at(g, 1));

  address_t const v3_is_2 = canonical_of(3, 0, 2);
  auto it = counts.find(v3_is_2);
  REQUIRE(it != counts.end());
  REQUIRE(it->second == 1);

  // Sanity: the v3=3 subgraph (handled by the v1-tail loop, which is
  // correct in both buggy and fixed versions) lands in its own bucket.
  address_t const v3_is_3 = canonical_of(3, 1, 0);
  auto it2 = counts.find(v3_is_3);
  REQUIRE(it2 != counts.end());
  REQUIRE(it2->second == 1);
}

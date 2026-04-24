// Unit tests for vcp::multirelational_directed_graph — the directed
// multi-relation CSR graph used by the vcp<n, r, true> specializations
// for r >= 2.
//
// Coverage parallels test_multirelational_graph.cpp: construction,
// istream/ostream round-trip, out/in-neighbor iteration, edge_value
// lookup, and relation_count. relation_count() had the same two latent
// bugs as the undirected counterpart (function-local `static` caching
// across instances, and std::max_element UB on empty edge_values).
// Fixes and test rationale are symmetric.

#include <catch2/catch_test_macros.hpp>

#include <sstream>

#include <vcp/multirelational_directed_graph.hpp>

TEST_CASE("multirelational_directed_graph default-constructs empty",
          "[multirelational_directed_graph]") {
  vcp::multirelational_directed_graph<2> g;
  REQUIRE(g.vertex_count() == 0);
  REQUIRE(g.out_edge_count() == 0);
  REQUIRE(g.in_edge_count() == 0);
  REQUIRE(g.relation_count() == 0);
}

TEST_CASE("multirelational_directed_graph loads a 3-cycle with a single active relation",
          "[multirelational_directed_graph]") {
  // 0 -> 1, 1 -> 2, 2 -> 0, every edge value 1.
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss("1,1\n2,1\n0,1\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.out_edge_count() == 3);
  REQUIRE(g.in_edge_count() == 3);
  REQUIRE(g.relation_count() == 1);
}

TEST_CASE("multirelational_directed_graph loads a 3-cycle with two active relations",
          "[multirelational_directed_graph]") {
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss("1,3\n2,3\n0,3\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.out_edge_count() == 3);
  REQUIRE(g.relation_count() == 2);
}

TEST_CASE("multirelational_directed_graph::relation_count is per-instance, not cached",
          "[multirelational_directed_graph][regression]") {
  // Same bug as the undirected counterpart: function-local `static`
  // cached across every instance of the same template instantiation.
  vcp::multirelational_directed_graph<2> g_low;
  {
    std::istringstream iss("1,1\n2,1\n0,1\n");
    iss >> g_low;
  }
  vcp::multirelational_directed_graph<2> g_high;
  {
    std::istringstream iss("1,3\n2,3\n0,3\n");
    iss >> g_high;
  }
  REQUIRE(g_low.relation_count() == 1);
  REQUIRE(g_high.relation_count() == 2);
  REQUIRE(g_high.relation_count() == 2);
  REQUIRE(g_low.relation_count() == 1);
}

TEST_CASE("multirelational_directed_graph out and in neighbor iteration",
          "[multirelational_directed_graph]") {
  // 0 -> 1 (value 2), 1 -> 2 (value 3), 2 -> 0 (value 1).
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss("1,2\n2,3\n0,1\n");
  iss >> g;
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v1 = v0 + 1;

  // v0 has a single out-neighbor (v1, value 2) and a single in-neighbor
  // (v2, value 1 — from the edge 2 -> 0).
  vcp::const_edge_iterator const out0 = g.out_neighbors_begin(v0);
  REQUIRE(g.out_neighbors_end(v0) - out0 == 1);
  REQUIRE(g.vertex_id(g.target_of(out0)) == 1);
  REQUIRE(g.edge_value(out0) == 2);

  vcp::const_edge_iterator const in0 = g.in_neighbors_begin(v0);
  REQUIRE(g.in_neighbors_end(v0) - in0 == 1);
  REQUIRE(g.vertex_id(g.target_of(in0)) == 2);
  REQUIRE(g.edge_value(in0) == 1);

  // v1's in-neighbor is v0 with value 2.
  vcp::const_edge_iterator const in1 = g.in_neighbors_begin(v1);
  REQUIRE(g.in_neighbors_end(v1) - in1 == 1);
  REQUIRE(g.vertex_id(g.target_of(in1)) == 0);
  REQUIRE(g.edge_value(in1) == 2);
}

TEST_CASE("multirelational_directed_graph::out_edge returns sentinel for non-adjacent",
          "[multirelational_directed_graph]") {
  // Single edge 0 -> 1 with value 1; other two vertices have no out-edges.
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss("1,1\n\n\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v2 = v0 + 2;
  auto const no_edge = g.out_edge(v0, v2);
  REQUIRE(no_edge == g.in_edges_end());
  // edge_value on the sentinel must return 0 (unconnected), not UB.
  REQUIRE(g.edge_value(no_edge) == 0);
}

TEST_CASE("multirelational_directed_graph round-trips through operator<< / operator>>",
          "[multirelational_directed_graph]") {
  std::string const input("1,2\n2,3\n0,1\n");
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss(input);
  iss >> g;
  std::ostringstream oss;
  oss << g;
  REQUIRE(oss.str() == input);
}

TEST_CASE("multirelational_directed_graph copy constructor deep-copies",
          "[multirelational_directed_graph]") {
  vcp::multirelational_directed_graph<2> g;
  std::istringstream iss("1,2\n2,3\n0,1\n");
  iss >> g;
  vcp::multirelational_directed_graph<2> g_copy(g);
  REQUIRE(g_copy.vertex_count() == g.vertex_count());
  REQUIRE(g_copy.out_edge_count() == g.out_edge_count());
  REQUIRE(g_copy.in_edge_count() == g.in_edge_count());
  REQUIRE(g_copy.relation_count() == g.relation_count());
  std::ostringstream o1, o2;
  o1 << g;
  o2 << g_copy;
  REQUIRE(o1.str() == o2.str());
}

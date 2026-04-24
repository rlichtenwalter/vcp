// Unit tests for vcp::multirelational_graph — the undirected multi-relation
// CSR graph used by the vcp<n, r, false> specializations for r >= 2.
//
// Coverage focuses on (a) construction and size accessors, (b) istream/
// ostream round-trip, (c) neighbor iteration + edge_value lookup, and
// (d) relation_count, which had two latent bugs fixed alongside this suite:
//   - function-local `static` cached the first call's result across every
//     instance of the same template instantiation (cross-instance bleed).
//   - std::max_element on an empty edge_values range produced UB on any
//     default-constructed or empty-input graph.
// The class has no internal callers of relation_count(), which is why the
// bugs survived: it's a public API surface with zero consumers in this
// codebase. These tests lock in the correct per-instance, empty-safe
// semantics for future callers.

#include <catch2/catch_test_macros.hpp>

#include <sstream>

#include <vcp/multirelational_graph.hpp>

TEST_CASE("multirelational_graph default-constructs empty", "[multirelational_graph]") {
  vcp::multirelational_graph<2> g;
  REQUIRE(g.vertex_count() == 0);
  REQUIRE(g.edge_count() == 0);
  // Empty graph has no relations to count. Previously, the static-local
  // cache would read past-the-end of a zero-size edge_values allocation.
  REQUIRE(g.relation_count() == 0);
}

TEST_CASE("multirelational_graph loads a triangle with a single active relation",
          "[multirelational_graph]") {
  vcp::multirelational_graph<2> g;
  std::istringstream iss("1,1 2,1\n0,1 2,1\n0,1 1,1\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.edge_count() == 3);
  // Every edge value is 1 = 0b01: only relation 0 is active.
  REQUIRE(g.relation_count() == 1);
}

TEST_CASE("multirelational_graph loads a triangle with two active relations",
          "[multirelational_graph]") {
  vcp::multirelational_graph<2> g;
  std::istringstream iss("1,3 2,3\n0,3 2,3\n0,3 1,3\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.edge_count() == 3);
  // Every edge value is 3 = 0b11: both relations 0 and 1 are active.
  REQUIRE(g.relation_count() == 2);
}

TEST_CASE("multirelational_graph::relation_count is per-instance, not cached",
          "[multirelational_graph][regression]") {
  // Regression for: relation_count() used a function-local `static`,
  // computing its value once and caching it across every instance of
  // the same template instantiation. Under the buggy implementation,
  // whichever graph called relation_count() first would imprint its
  // answer on every subsequent graph of the same <r>.
  vcp::multirelational_graph<2> g_low;
  {
    std::istringstream iss("1,1 2,1\n0,1 2,1\n0,1 1,1\n");
    iss >> g_low;
  }
  vcp::multirelational_graph<2> g_high;
  {
    std::istringstream iss("1,3 2,3\n0,3 2,3\n0,3 1,3\n");
    iss >> g_high;
  }
  REQUIRE(g_low.relation_count() == 1);
  REQUIRE(g_high.relation_count() == 2);
  // Query again in reverse order: the answers must still be per-instance.
  REQUIRE(g_high.relation_count() == 2);
  REQUIRE(g_low.relation_count() == 1);
}

TEST_CASE("multirelational_graph round-trips through operator<< / operator>>",
          "[multirelational_graph]") {
  std::string const input("1,2 2,3\n0,2 2,1\n0,3 1,1\n");
  vcp::multirelational_graph<2> g;
  std::istringstream iss(input);
  iss >> g;
  std::ostringstream oss;
  oss << g;
  REQUIRE(oss.str() == input);
}

TEST_CASE("multirelational_graph neighbor iteration yields correct targets and values",
          "[multirelational_graph]") {
  // Vertex 0 adjacencies: (v1, value 2), (v2, value 3).
  vcp::multirelational_graph<2> g;
  std::istringstream iss("1,2 2,3\n0,2 2,1\n0,3 1,1\n");
  iss >> g;
  vcp::const_vertex_iterator const vIt = g.vertices_begin();
  vcp::const_edge_iterator const nBegin = g.neighbors_begin(vIt);
  vcp::const_edge_iterator const nEnd = g.neighbors_end(vIt);
  REQUIRE(nEnd - nBegin == 2);
  REQUIRE(g.vertex_id(g.target_of(nBegin)) == 1);
  REQUIRE(g.edge_value(nBegin) == 2);
  REQUIRE(g.vertex_id(g.target_of(nBegin + 1)) == 2);
  REQUIRE(g.edge_value(nBegin + 1) == 3);
}

TEST_CASE("multirelational_graph::edge returns edges_end for non-adjacent vertices",
          "[multirelational_graph]") {
  // Only edge is 0 <-> 1 (value 1). Vertex 2 is isolated.
  vcp::multirelational_graph<2> g;
  std::istringstream iss("1,1\n0,1\n\n");
  iss >> g;
  REQUIRE(g.vertex_count() == 3);
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v2 = v0 + 2;
  auto const it = g.edge(v0, v2);
  REQUIRE(it == g.edges_end());
  // edge_value on the sentinel must return 0 (unconnected), not UB.
  REQUIRE(g.edge_value(it) == 0);
}

TEST_CASE("multirelational_graph copy constructor deep-copies", "[multirelational_graph]") {
  vcp::multirelational_graph<2> g;
  std::istringstream iss("1,1 2,2\n0,1 2,3\n0,2 1,3\n");
  iss >> g;
  vcp::multirelational_graph<2> g_copy(g);
  REQUIRE(g_copy.vertex_count() == g.vertex_count());
  REQUIRE(g_copy.edge_count() == g.edge_count());
  REQUIRE(g_copy.relation_count() == g.relation_count());
  std::ostringstream o1, o2;
  o1 << g;
  o2 << g_copy;
  REQUIRE(o1.str() == o2.str());
}

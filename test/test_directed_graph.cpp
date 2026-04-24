// Unit tests for vcp::directed_graph — the foundational directed CSR
// graph used by the vcp<n, 1, true> specializations.
//
// Coverage parallels test_graph.cpp, with additional checks for the
// in-neighbor derivation that operator>> performs after parsing
// out-neighbors (unique to the directed class):
//   - out/in-neighbor iteration is symmetric: if u → v is in u's out
//     list, v → u is in v's in list with value = out-side value.
//   - out_edge() / in_edge() both return in_edges_end() for non-adjacent
//     pairs (shared sentinel convention).
//   - Copy constructor preserves the 2*vertex_count + 1 slot layout
//     (out starts / in starts / end sentinel).

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

#include <vcp/directed_graph.hpp>

namespace {

vcp::directed_graph parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::directed_graph g;
  iss >> g;
  return g;
}

} // namespace

TEST_CASE("directed_graph default-constructs empty", "[directed_graph]") {
  vcp::directed_graph g;
  REQUIRE(g.vertex_count() == 0);
  REQUIRE(g.out_edge_count() == 0);
  REQUIRE(g.in_edge_count() == 0);
}

TEST_CASE("directed_graph parses an oriented 3-cycle", "[directed_graph]") {
  // 0 -> 1, 1 -> 2, 2 -> 0. Each vertex has exactly one out- and one
  // in-neighbor.
  vcp::directed_graph g = parse("1\n2\n0\n");
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.out_edge_count() == 3);
  REQUIRE(g.in_edge_count() == 3);
}

TEST_CASE("directed_graph derives in-neighbors from the parsed out list", "[directed_graph]") {
  // 0 -> 1, 1 -> 2, 2 -> 0. Vertex 1's in-neighbor should be vertex 0.
  vcp::directed_graph g = parse("1\n2\n0\n");
  vcp::const_vertex_iterator const v1 = g.vertices_begin() + 1;
  vcp::const_edge_iterator const in_begin = g.in_neighbors_begin(v1);
  vcp::const_edge_iterator const in_end = g.in_neighbors_end(v1);
  REQUIRE(in_end - in_begin == 1);
  REQUIRE(g.vertex_id(g.target_of(in_begin)) == 0);
}

TEST_CASE("directed_graph out_edge and in_edge sentinel share in_edges_end()", "[directed_graph]") {
  // Sole edge 0 -> 1. Vertex 2 is isolated. out_edge(0, 2) must return
  // the shared sentinel (==in_edges_end), and in_edge(0, 2) likewise.
  // The sentinel convention is: both out_edge and in_edge return
  // in_edges_end() on miss.
  vcp::directed_graph g = parse("1\n\n\n");
  REQUIRE(g.vertex_count() == 3);
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v2 = v0 + 2;
  auto const out_miss = g.out_edge(v0, v2);
  auto const in_miss = g.in_edge(v0, v2);
  REQUIRE(out_miss == g.in_edges_end());
  REQUIRE(in_miss == g.in_edges_end());
  REQUIRE_FALSE(g.out_edge_exists(v0, v2));
  REQUIRE_FALSE(g.in_edge_exists(v0, v2));
}

TEST_CASE("directed_graph out/in asymmetry: a one-way edge is seen from one side only",
          "[directed_graph]") {
  // 0 -> 1 only. out_edge_exists(0, 1) is true; in_edge_exists(0, 1) is
  // false. The reverse (0 <- 1) is not present: in_edge_exists(1, 0) is
  // true while out_edge_exists(1, 0) is false.
  vcp::directed_graph g = parse("1\n\n");
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v1 = v0 + 1;
  REQUIRE(g.out_edge_exists(v0, v1));
  REQUIRE_FALSE(g.in_edge_exists(v0, v1));
  REQUIRE(g.in_edge_exists(v1, v0));
  REQUIRE_FALSE(g.out_edge_exists(v1, v0));
}

TEST_CASE("directed_graph mutual edge: both out and in present on both endpoints",
          "[directed_graph]") {
  vcp::directed_graph g = parse("1\n0\n");
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v1 = v0 + 1;
  REQUIRE(g.out_edge_exists(v0, v1));
  REQUIRE(g.in_edge_exists(v0, v1));
  REQUIRE(g.out_edge_exists(v1, v0));
  REQUIRE(g.in_edge_exists(v1, v0));
}

TEST_CASE("directed_graph round-trips through operator<< / operator>>", "[directed_graph]") {
  // operator<< emits only out-neighbors (in-neighbors are derived on
  // parse). Round-trip must equal the original string.
  std::string const input("1 2\n2\n0\n");
  vcp::directed_graph g = parse(input);
  std::ostringstream oss;
  oss << g;
  REQUIRE(oss.str() == input);
}

TEST_CASE("directed_graph operator>> rejects malformed vertex id tokens",
          "[directed_graph][parse-error]") {
  std::istringstream iss("1 bogus\n\n");
  vcp::directed_graph g;
  REQUIRE_THROWS_AS(iss >> g, std::invalid_argument);
}

TEST_CASE("directed_graph operator>> rejects partial-prefix vertex ids",
          "[directed_graph][parse-error]") {
  std::istringstream iss("1 2x\n\n");
  vcp::directed_graph g;
  REQUIRE_THROWS_AS(iss >> g, std::invalid_argument);
}

TEST_CASE("directed_graph copy constructor deep-copies out and in layouts", "[directed_graph]") {
  vcp::directed_graph g = parse("1 2\n2\n0\n");
  // Exercising the copy ctor is the point of this test; clang-tidy's
  // performance-unnecessary-copy-initialization can't see that.
  // NOLINTNEXTLINE(performance-unnecessary-copy-initialization)
  vcp::directed_graph g_copy(g);
  REQUIRE(g_copy.vertex_count() == g.vertex_count());
  REQUIRE(g_copy.out_edge_count() == g.out_edge_count());
  REQUIRE(g_copy.in_edge_count() == g.in_edge_count());
  std::ostringstream a, b;
  a << g;
  b << g_copy;
  REQUIRE(a.str() == b.str());
  // In-neighbor correctness on the copy: vertex 2's in-neighbors should
  // be vertex 0 (from 0 -> 2) and vertex 1 (from 1 -> 2). Validates that
  // the copy preserves the second half of the vertices buffer (the
  // in-neighbor start pointers) and not just the out half.
  vcp::const_vertex_iterator const v2 = g_copy.vertices_begin() + 2;
  vcp::const_edge_iterator const in_begin = g_copy.in_neighbors_begin(v2);
  vcp::const_edge_iterator const in_end = g_copy.in_neighbors_end(v2);
  REQUIRE(in_end - in_begin == 2);
}

TEST_CASE("directed_graph copy assignment deep-copies", "[directed_graph]") {
  vcp::directed_graph g = parse("1\n2\n0\n");
  vcp::directed_graph dst;
  dst = g;
  REQUIRE(dst.vertex_count() == 3);
  REQUIRE(dst.out_edge_count() == 3);
  std::ostringstream a, b;
  a << g;
  b << dst;
  REQUIRE(a.str() == b.str());
}

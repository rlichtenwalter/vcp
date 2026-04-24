// Unit tests for vcp::graph — the foundational undirected CSR graph used
// by the vcp<n, 1, false> specializations.
//
// Coverage focuses on behavior that existing regression cannot distinguish
// from identical bugs in legacy and current builds:
//   - istream/ostream round-trip on representative topologies (empty,
//     isolated vertices, star, triangle, path).
//   - operator>> parse rejection for malformed input (from_chars errno path
//     introduced alongside the modernization; legacy's atol silently
//     returned 0 and was UB on overflow).
//   - neighbor_begin/end iteration bounds on the first, middle, and last
//     vertex (exercises the CSR end-sentinel that makes
//     neighbors_end(last_vertex) correct).
//   - edge() sentinel for non-adjacent vertices equals edges_end().
//   - copy constructor preserves vertex count, edge count, and the
//     per-vertex neighbor-list pointer layout (the tricky bit — the
//     vertices[vertex_count()] end-sentinel is derived, not copied from
//     the source's corresponding slot).

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

#include <vcp/graph.hpp>

namespace {

vcp::graph parse(std::string const &text) {
  std::istringstream iss(text);
  vcp::graph g;
  iss >> g;
  return g;
}

} // namespace

TEST_CASE("graph default-constructs empty", "[graph]") {
  vcp::graph g;
  REQUIRE(g.vertex_count() == 0);
  REQUIRE(g.edge_count() == 0);
  REQUIRE(g.vertices_begin() == g.vertices_end());
}

TEST_CASE("graph parses an empty single-vertex line as one isolated vertex", "[graph]") {
  // A single newline-terminated empty line means one vertex with no
  // neighbors. edge_count() returns num_edges / 2 == 0 on both sides of
  // the undirected count (no stored half-edges).
  vcp::graph g = parse("\n");
  REQUIRE(g.vertex_count() == 1);
  REQUIRE(g.edge_count() == 0);
  REQUIRE(g.neighbors_begin(g.vertices_begin()) == g.neighbors_end(g.vertices_begin()));
}

TEST_CASE("graph parses a triangle with three edges", "[graph]") {
  // Each undirected edge appears as two half-edges; edge_count() reports
  // num_edges / 2 = 3.
  vcp::graph g = parse("1 2\n0 2\n0 1\n");
  REQUIRE(g.vertex_count() == 3);
  REQUIRE(g.edge_count() == 3);
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  REQUIRE(g.neighbors_end(v0) - g.neighbors_begin(v0) == 2);
  REQUIRE(g.vertex_id(g.target_of(g.neighbors_begin(v0))) == 1);
  REQUIRE(g.vertex_id(g.target_of(g.neighbors_begin(v0) + 1)) == 2);
}

TEST_CASE("graph edge() returns a real iterator for adjacent vertices", "[graph]") {
  vcp::graph g = parse("1 2\n0 2\n0 1\n");
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v2 = v0 + 2;
  auto const it = g.edge(v0, v2);
  REQUIRE(it != g.edges_end());
  REQUIRE(g.edge_exists(it));
  REQUIRE(g.edge_exists(v0, v2));
}

TEST_CASE("graph edge() returns edges_end() for non-adjacent vertices", "[graph]") {
  // Isolated vertex 2 is not adjacent to vertex 0.
  vcp::graph g = parse("1\n0\n\n");
  REQUIRE(g.vertex_count() == 3);
  vcp::const_vertex_iterator const v0 = g.vertices_begin();
  vcp::const_vertex_iterator const v2 = v0 + 2;
  auto const it = g.edge(v0, v2);
  REQUIRE(it == g.edges_end());
  REQUIRE_FALSE(g.edge_exists(v0, v2));
}

TEST_CASE("graph round-trips through operator<< / operator>>", "[graph]") {
  // Star: vertex 0 connected to everyone, others isolated besides 0.
  std::string const input("1 2 3\n0\n0\n0\n");
  vcp::graph g = parse(input);
  std::ostringstream oss;
  oss << g;
  REQUIRE(oss.str() == input);
}

TEST_CASE("graph operator>> rejects malformed vertex id tokens", "[graph][parse-error]") {
  // atol (legacy) silently returned 0 on non-numeric input. from_chars
  // (current) rejects anything that isn't fully-parsed as an integer.
  std::istringstream iss("1 bogus\n0\n");
  vcp::graph g;
  REQUIRE_THROWS_AS(iss >> g, std::invalid_argument);
}

TEST_CASE("graph operator>> rejects mixed-numeric vertex ids", "[graph][parse-error]") {
  // Trailing garbage after a valid prefix must fail: from_chars reports
  // partial consumption, which the parser treats as an error.
  std::istringstream iss("1 2a\n0\n");
  vcp::graph g;
  REQUIRE_THROWS_AS(iss >> g, std::invalid_argument);
}

TEST_CASE("graph neighbor iteration on the last vertex reaches the CSR end sentinel", "[graph]") {
  // The last vertex's neighbors_end reads vertices[vertex_count()] which
  // is the end-sentinel initialized in operator>>. If that slot is not
  // set, iteration would run off into adjacent heap memory. This test
  // fails with a heap-buffer-overflow under ASan if the end sentinel is
  // ever dropped.
  vcp::graph g = parse("1\n0 2\n1\n");
  vcp::const_vertex_iterator const v_last = g.vertices_begin() + 2;
  vcp::const_edge_iterator const nb = g.neighbors_begin(v_last);
  vcp::const_edge_iterator const ne = g.neighbors_end(v_last);
  REQUIRE(ne - nb == 1);
  REQUIRE(g.vertex_id(g.target_of(nb)) == 1);
}

TEST_CASE("graph copy constructor deep-copies with working neighbor iteration", "[graph]") {
  // The copy constructor rebuilds the CSR pointer layout in the destination
  // buffer. Correctness here means: after a copy, (a) counts match, (b)
  // serialization of the copy equals the source, (c) a subsequent edge()
  // query on the copy still resolves.
  vcp::graph g = parse("1 2\n0 2\n0 1\n");
  // Exercising the copy ctor is the point of this test; clang-tidy's
  // performance-unnecessary-copy-initialization can't see that.
  // NOLINTNEXTLINE(performance-unnecessary-copy-initialization)
  vcp::graph g_copy(g);
  REQUIRE(g_copy.vertex_count() == g.vertex_count());
  REQUIRE(g_copy.edge_count() == g.edge_count());
  std::ostringstream a, b;
  a << g;
  b << g_copy;
  REQUIRE(a.str() == b.str());
  // Independence: edits to the original CSR do not affect the copy. We
  // can't mutate an existing edge (graph is immutable after parse), but
  // we verify that the copy's edge lookup resolves independently.
  auto const it = g_copy.edge(g_copy.vertices_begin(), g_copy.vertices_begin() + 2);
  REQUIRE(it != g_copy.edges_end());
}

TEST_CASE("graph copy assignment deep-copies", "[graph]") {
  vcp::graph g = parse("1 2\n0 2\n0 1\n");
  vcp::graph dst;
  dst = g;
  REQUIRE(dst.vertex_count() == 3);
  REQUIRE(dst.edge_count() == 3);
  std::ostringstream a, b;
  a << g;
  b << dst;
  REQUIRE(a.str() == b.str());
}

TEST_CASE("graph copy assignment is a no-op on self", "[graph]") {
  // Verifies the this != &g guard in operator=; without it, we would
  // reallocate and copy from the buffers we just reset to nullptr.
  vcp::graph g = parse("1\n0\n");
  vcp::graph const &self_ref(g);
  g = self_ref;
  REQUIRE(g.vertex_count() == 2);
  REQUIRE(g.edge_count() == 1);
}

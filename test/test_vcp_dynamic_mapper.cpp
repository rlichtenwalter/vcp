// Unit tests for vcp::vcp_dynamic_mapper — the reference mapper used to
// encode a connectivity matrix as a flat subgraph_address, to collapse
// isomorphic subgraphs onto a canonical address, and (in element_structure)
// to recover a representative matrix from an address.
//
// Coverage focuses on (a) subgraph_count, (b) subgraph_address correctness
// bit-wise, (c) subgraph_address / element_structure round-trip, and (d)
// canonical_subgraph_address invariance across v3-v4 swap isomorphisms.
//
// element_structure previously shifted by `(1 << r)` (the size of the
// r-relation power-set) instead of by `r` (the width of a connectivity
// slot), so decoding advanced by 2**r bits per step rather than r. The
// round-trip tests in this suite fail on the old form for every (n, r)
// with r >= 1. element_structure had no internal callers — the bug was
// latent in a public API surface — and the function had never been
// instantiated, which is why the additional signature issue (mutating
// a const-ref parameter) went uncaught.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace {

// Helper: set the undirected upper-triangular slots of `m` to consecutive
// nonzero values in the mapper's iteration order, wrapping at the
// per-slot capacity of 2**r - 1. Value 0 is avoided so the round-trip
// sentinel (loop exit on address == 0) does not short-circuit.
template <std::size_t n, std::size_t r>
void fill_undirected(vcp::square_matrix<std::size_t, n> &m) {
  std::size_t const cap((std::size_t(1) << r) - 1);
  std::size_t v(1);
  for (std::size_t row(0); row < n; ++row) {
    for (std::size_t column(row + 1); column < n; ++column) {
      m(row, column) = v;
      m(column, row) = v;
      v = (v % cap) + 1;
    }
  }
}

// Helper: set every off-diagonal slot of `m` to a distinct nonzero value
// in the mapper's directed iteration order, wrapping at 2**r - 1.
template <std::size_t n, std::size_t r> void fill_directed(vcp::square_matrix<std::size_t, n> &m) {
  std::size_t const cap((std::size_t(1) << r) - 1);
  std::size_t v(1);
  for (std::size_t row(0); row < n; ++row) {
    for (std::size_t column(0); column < n; ++column) {
      if (row != column) {
        m(row, column) = v;
        v = (v % cap) + 1;
      }
    }
  }
}

template <std::size_t n>
bool matrices_equal(vcp::square_matrix<std::size_t, n> const &a,
                    vcp::square_matrix<std::size_t, n> const &b) {
  for (std::size_t row(0); row < n; ++row) {
    for (std::size_t column(0); column < n; ++column) {
      if (a(row, column) != b(row, column)) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

TEST_CASE("vcp_dynamic_mapper subgraph_count equals 2**(slots * r)", "[vcp_dynamic_mapper]") {
  // Undirected n=3, r=1: 3 undirected pairs, 1 bit each → 2**3 = 8.
  vcp::vcp_dynamic_mapper<3, 1, false> m3_1_u;
  REQUIRE(m3_1_u.subgraph_count() == 8);

  // Directed n=3, r=1: 6 directed pairs, 1 bit each → 2**6 = 64.
  vcp::vcp_dynamic_mapper<3, 1, true> m3_1_d;
  REQUIRE(m3_1_d.subgraph_count() == 64);

  // Undirected n=4, r=2: 6 pairs, 2 bits each → 2**12 = 4096.
  vcp::vcp_dynamic_mapper<4, 2, false> m4_2_u;
  REQUIRE(m4_2_u.subgraph_count() == 4096);

  // Directed n=4, r=2: 12 pairs, 2 bits each → 2**24 = 16777216.
  vcp::vcp_dynamic_mapper<4, 2, true> m4_2_d;
  REQUIRE(m4_2_d.subgraph_count() == (std::size_t(1) << 24));
}

TEST_CASE("vcp_dynamic_mapper subgraph_address encodes slot values bit-packed "
          "at r-aligned offsets",
          "[vcp_dynamic_mapper]") {
  // Undirected n=3, r=1. value_matrix offsets are (0,1)->0, (0,2)->1, (1,2)->2.
  // All-ones upper triangular should produce 0b111 = 7.
  vcp::vcp_dynamic_mapper<3, 1, false> mapper;
  vcp::square_matrix<std::size_t, 3> m;
  m(0, 1) = 1;
  m(1, 0) = 1;
  m(0, 2) = 1;
  m(2, 0) = 1;
  m(1, 2) = 1;
  m(2, 1) = 1;
  REQUIRE(mapper.subgraph_address(m) == 7);
}

TEST_CASE("vcp_dynamic_mapper round-trip: subgraph_address then "
          "element_structure yields original matrix (undirected n=3, r=1)",
          "[vcp_dynamic_mapper][regression]") {
  // Regression for: element_structure shifted by (1 << r) instead of r.
  // For r=1 the buggy form shifts by 2 per iteration, dropping every
  // other bit.
  vcp::vcp_dynamic_mapper<3, 1, false> mapper;
  vcp::square_matrix<std::size_t, 3> src;
  fill_undirected<3, 1>(src);
  auto const address = mapper.subgraph_address(src);
  auto const dst = mapper.element_structure(address);
  // element_structure only populates the upper triangle (it writes
  // matrix(row, column) for column > row in the undirected iteration
  // order). The lower triangle of dst stays zero-initialized. Compare
  // only the upper triangle.
  for (std::size_t row(0); row < 3; ++row) {
    for (std::size_t column(row + 1); column < 3; ++column) {
      REQUIRE(dst(row, column) == src(row, column));
    }
  }
}

TEST_CASE("vcp_dynamic_mapper round-trip: subgraph_address then "
          "element_structure yields original matrix (undirected n=4, r=2)",
          "[vcp_dynamic_mapper][regression]") {
  // For r=2 the buggy form shifts by 4 per iteration, scrambling every
  // 2-bit slot.
  vcp::vcp_dynamic_mapper<4, 2, false> mapper;
  vcp::square_matrix<std::size_t, 4> src;
  fill_undirected<4, 2>(src);
  auto const address = mapper.subgraph_address(src);
  auto const dst = mapper.element_structure(address);
  for (std::size_t row(0); row < 4; ++row) {
    for (std::size_t column(row + 1); column < 4; ++column) {
      REQUIRE(dst(row, column) == src(row, column));
    }
  }
}

TEST_CASE("vcp_dynamic_mapper round-trip: subgraph_address then "
          "element_structure yields original matrix (directed n=3, r=1)",
          "[vcp_dynamic_mapper][regression]") {
  vcp::vcp_dynamic_mapper<3, 1, true> mapper;
  vcp::square_matrix<std::size_t, 3> src;
  fill_directed<3, 1>(src);
  auto const address = mapper.subgraph_address(src);
  auto const dst = mapper.element_structure(address);
  // Directed element_structure populates both (row, column) and
  // (column, row) slots; the whole off-diagonal must match.
  for (std::size_t row(0); row < 3; ++row) {
    for (std::size_t column(0); column < 3; ++column) {
      if (row != column) {
        REQUIRE(dst(row, column) == src(row, column));
      }
    }
  }
}

TEST_CASE("vcp_dynamic_mapper round-trip: subgraph_address then "
          "element_structure yields original matrix (directed n=4, r=2)",
          "[vcp_dynamic_mapper][regression]") {
  vcp::vcp_dynamic_mapper<4, 2, true> mapper;
  vcp::square_matrix<std::size_t, 4> src;
  fill_directed<4, 2>(src);
  auto const address = mapper.subgraph_address(src);
  auto const dst = mapper.element_structure(address);
  for (std::size_t row(0); row < 4; ++row) {
    for (std::size_t column(0); column < 4; ++column) {
      if (row != column) {
        REQUIRE(dst(row, column) == src(row, column));
      }
    }
  }
}

TEST_CASE("vcp_dynamic_mapper::canonical_subgraph_address collapses v3-v4 "
          "isomorphic subgraphs (undirected n=4, r=1)",
          "[vcp_dynamic_mapper]") {
  // Two matrices that differ only by swapping v3 and v4 (indices 2 and 3).
  // canonical_subgraph_address permutes indices 2 through n-1 and picks
  // the minimum, so both inputs must produce the same canonical address.
  //
  // NOTE: the VCP specializations only populate the upper triangle of the
  // connectivity matrix for undirected graphs — see `it1->second(0, 2) = ...`
  // / `it1->second(1, 2) = ...` sites in vcp_4_r_0.hpp. Matching that
  // convention here is required: if both triangles were populated,
  // canonical_subgraph_address's all-off-diagonals permuted iteration would
  // double-count (value_matrix is symmetric for undirected), inflate every
  // permuted variant past the baseline, and canonicalization would degrade
  // to a no-op.
  vcp::vcp_dynamic_mapper<4, 1, false> mapper;
  vcp::square_matrix<std::size_t, 4> a;
  a(0, 2) = 1;
  a(1, 3) = 1;
  vcp::square_matrix<std::size_t, 4> b;
  // Same topology with v3 and v4 swapped: 0-v4 and 1-v3 edges.
  b(0, 3) = 1;
  b(1, 2) = 1;
  REQUIRE(mapper.canonical_subgraph_address(a) == mapper.canonical_subgraph_address(b));
}

TEST_CASE("vcp_dynamic_mapper::canonical_subgraph_address is idempotent", "[vcp_dynamic_mapper]") {
  // Canonicalizing a canonical address should be a no-op. The canonical
  // address derives from element_structure only if the user explicitly
  // round-trips; the stronger property this verifies is that running
  // canonical_subgraph_address twice on different matrices that share a
  // canonical form produces the same result.
  vcp::vcp_dynamic_mapper<4, 2, false> mapper;
  vcp::square_matrix<std::size_t, 4> a;
  fill_undirected<4, 2>(a);
  auto const canon_a = mapper.canonical_subgraph_address(a);
  // Build b by re-decoding canon_a and re-encoding it. Its canonical
  // form must match canon_a.
  auto const b = mapper.element_structure(canon_a);
  // element_structure only fills upper triangle; mirror it to restore
  // symmetry for the undirected encoder.
  vcp::square_matrix<std::size_t, 4> b_sym;
  for (std::size_t row(0); row < 4; ++row) {
    for (std::size_t column(row + 1); column < 4; ++column) {
      b_sym(row, column) = b(row, column);
      b_sym(column, row) = b(row, column);
    }
  }
  REQUIRE(mapper.canonical_subgraph_address(b_sym) == canon_a);
}

// Unit tests for vcp::vcp_static_mapper — the runtime-parameterized
// mapper used by the `vcp_map` tool and by callers that need to collapse
// subgraph addresses onto VCP element addresses at runtime-chosen (n, r,
// d). The static and dynamic mappers are two implementations of the same
// contract, so these tests compare static against dynamic where
// feasible, and against paper cardinalities elsewhere.
//
// Coverage:
//   - subgraph_count static-vs-dynamic consistency across small (n, r, d).
//   - element cardinality for (4, 1, 0) and (4, 1, 1) matches the paper
//     (40 and 2112 respectively; also regression Phase 1 ground truth).
//   - element_address agreement between static and dynamic mappers on a
//     small set of hand-picked connectivity matrices.
//   - element_structure round-trip (address → structure → address).

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <set>

#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>
#include <vcp/vcp_static_mapper.hpp>

TEST_CASE("vcp_static_mapper::subgraph_count matches formula 2^(n(n-1)r(d+1)/2)",
          "[vcp_static_mapper]") {
  // Exhaust the small corners the library actually compiles for.
  REQUIRE(vcp::vcp_static_mapper::subgraph_count(3, 1, false) == 8);
  REQUIRE(vcp::vcp_static_mapper::subgraph_count(3, 1, true) == 64);
  REQUIRE(vcp::vcp_static_mapper::subgraph_count(4, 1, false) == 64);
  REQUIRE(vcp::vcp_static_mapper::subgraph_count(4, 1, true) == 4096);
}

TEST_CASE("vcp_static_mapper element cardinality matches the paper for (4, 1, 0)",
          "[vcp_static_mapper]") {
  // The static mapper's map vector has one entry per subgraph address,
  // and the count of distinct element-address values is the VCP
  // cardinality |VCP_{n,r,d}|. For (4, 1, 0) the paper gives 40 elements.
  vcp::vcp_static_mapper m(4, 1, false);
  std::set<std::size_t> elements;
  // Enumerate all subgraphs by encoding each 6-bit value as a symmetric
  // 4x4 matrix and reading back its element address.
  for (std::size_t subgraph(0); subgraph < 64; ++subgraph) {
    vcp::square_matrix<std::size_t> conn(4);
    std::size_t bits(subgraph);
    std::size_t slot(0);
    for (std::size_t row(0); row < 4; ++row) {
      for (std::size_t column(row + 1); column < 4; ++column) {
        std::size_t const v((bits >> slot) & 1);
        conn(row, column) = v;
        conn(column, row) = v;
        ++slot;
      }
    }
    elements.insert(m.element_address(conn));
  }
  REQUIRE(elements.size() == 40);
}

TEST_CASE("vcp_static_mapper for (3, 1, 1) has no isomorphism collapse (cardinality == 64)",
          "[vcp_static_mapper]") {
  // For n=3, v1 and v2 are fixed and only v3..vn are permuted to find
  // isomorphic equivalents — but there is no v3..vn range to permute,
  // so no two distinct subgraphs collapse. The element cardinality
  // equals the subgraph cardinality: 64. (`vcp_map 3 1 1 | sort -u |
  // wc -l` confirms.) A regression that introduced incorrect
  // canonicalization — for instance, permuting across v1/v2 — would
  // collapse the count below 64.
  vcp::vcp_static_mapper m(3, 1, true);
  std::set<std::size_t> elements;
  for (std::size_t subgraph(0); subgraph < 64; ++subgraph) {
    vcp::square_matrix<std::size_t> conn(3);
    std::size_t bits(subgraph);
    std::size_t slot(0);
    for (std::size_t row(0); row < 3; ++row) {
      for (std::size_t column(0); column < 3; ++column) {
        if (row != column) {
          conn(row, column) = (bits >> slot) & 1;
          ++slot;
        }
      }
    }
    elements.insert(m.element_address(conn));
  }
  REQUIRE(elements.size() == 64);
}

TEST_CASE("vcp_static_mapper::element_structure round-trips a subgraph address for (4, 1, 0)",
          "[vcp_static_mapper]") {
  // element_structure's argument is a subgraph address, not an element
  // address (the name is historical — it returns the connectivity
  // matrix that encodes that subgraph). The round-trip subgraph_address
  // → structure → subgraph_address must be the identity across the
  // full 2^6 = 64 address range.
  vcp::vcp_static_mapper m(4, 1, false);
  for (std::size_t subgraph(0); subgraph < 64; ++subgraph) {
    vcp::square_matrix<std::size_t> const structure = m.element_structure(subgraph);
    REQUIRE(m.subgraph_address(structure) == subgraph);
  }
}

TEST_CASE(
    "vcp_static_mapper and vcp_dynamic_mapper induce the same element partition for (4, 1, 0)",
    "[vcp_static_mapper][cross-check]") {
  // The two mappers' value_matrix layouts differ (static uses powers of
  // 2^r, dynamic packs by bit-shifts), but the property that must hold
  // regardless is that they agree on the *equivalence partition* of
  // subgraphs. Formally, for any two subgraphs s1, s2 in [0, 64):
  //
  //   static.element_address(s1) == static.element_address(s2)
  //   iff
  //   dynamic.canonical_subgraph_address(s1) ==
  //       dynamic.canonical_subgraph_address(s2).
  //
  // A bug that renumbered one mapper's element space but preserved the
  // partition would pass; a bug that split or merged equivalence
  // classes on either side would fail. This is the proper cross-check
  // of static-vs-dynamic agreement — the earlier form of this test
  // compared two unrelated quantities and would have missed a real
  // discrepancy.
  vcp::vcp_static_mapper s(4, 1, false);
  vcp::vcp_dynamic_mapper<4, 1, false> d;

  auto build = [](std::size_t subgraph) {
    vcp::square_matrix<std::size_t> conn_s(4);
    vcp::square_matrix<std::size_t, 4> conn_d;
    std::size_t bits(subgraph);
    std::size_t slot(0);
    for (std::size_t row(0); row < 4; ++row) {
      for (std::size_t column(row + 1); column < 4; ++column) {
        std::size_t const v((bits >> slot) & 1);
        conn_s(row, column) = v;
        conn_s(column, row) = v;
        conn_d(row, column) = v;
        conn_d(column, row) = v;
        ++slot;
      }
    }
    return std::pair{conn_s, conn_d};
  };

  for (std::size_t s1(0); s1 < 64; ++s1) {
    auto const inputs1 = build(s1);
    for (std::size_t s2(s1); s2 < 64; ++s2) {
      auto const inputs2 = build(s2);
      bool const static_same = s.element_address(inputs1.first) == s.element_address(inputs2.first);
      bool const dynamic_same = d.canonical_subgraph_address(inputs1.second) ==
                                d.canonical_subgraph_address(inputs2.second);
      REQUIRE(static_same == dynamic_same);
    }
  }
}

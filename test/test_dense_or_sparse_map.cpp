// Unit tests for vcp::detail::dense_or_sparse_map — the two-tier
// std::map replacement used internally by the vcp<4, r, d> hot path
// for the temp_edge_types and edge_types maps.
//
// Coverage focuses on the behaviors the hot path actually depends on:
//   - insert_or_zero correctly creates a zeroed entry on miss and
//     returns a stable reference on hit (so `++m.insert_or_zero(k)`
//     is a proper insert-or-increment).
//   - find distinguishes legitimately-zero entries from never-inserted
//     ones (the footgun the previous iteration of dense_impl hit: a
//     value==0 sentinel would confuse gaps' zero value with absence,
//     silently dropping the unconnected-class correction in the tail
//     loop of vcp_4_r_0/1.hpp).
//   - clear() resets presence without disturbing the backing array,
//     so a reused class-member instance behaves like a fresh one.
//   - for_each walks only the populated entries and does not emit
//     stale data from before the last clear().
//   - The compile-time tier choice is the expected one for each
//     KeyCount threshold (dense up to DENSE_BYTE_BUDGET / sizeof(V),
//     sparse above).
//   - The directed-pair packer (pack_pair_key_by_bits) collapses
//     (lo, hi) into a single index and re-inflates correctly for
//     find/for_each.

#include <catch2/catch_test_macros.hpp>

#include <climits>
#include <cstddef>
#include <utility>

#include <vcp/detail/dense_or_sparse_map.hpp>

using vcp::detail::dense_or_sparse_map;
using vcp::detail::pack_pair_key_by_bits;

TEST_CASE("dense_or_sparse_map selects the dense tier under the byte budget",
          "[dense_or_sparse_map]") {
  // Value = unsigned long (8 bytes), KeyCount = 4 ⇒ 32 bytes ≪ 8 MB ⇒ dense.
  static_assert((dense_or_sparse_map<std::size_t, unsigned long, 4>::is_dense));

  // Value = unsigned long (8 bytes), KeyCount = 1 << 20 (1M) ⇒ 8 MB == budget ⇒ dense.
  static_assert(
      (dense_or_sparse_map<std::size_t, unsigned long, (std::size_t(1) << 20)>::is_dense));

  // Value = unsigned long (8 bytes), KeyCount = (1 << 20) + 1 ⇒ just over budget ⇒ sparse.
  static_assert(
      !(dense_or_sparse_map<std::size_t, unsigned long, (std::size_t(1) << 20) + 1>::is_dense));

  // KeyCount = 0 is the "cannot statically size" sentinel ⇒ sparse.
  static_assert(!(dense_or_sparse_map<std::size_t, unsigned long, 0>::is_dense));
}

TEST_CASE("dense_or_sparse_map::insert_or_zero creates a zeroed entry on miss",
          "[dense_or_sparse_map][dense]") {
  dense_or_sparse_map<std::size_t, unsigned long, 4> m;
  unsigned long &slot = m.insert_or_zero(2);
  REQUIRE(slot == 0);
  REQUIRE(m.size() == 1);
  // Write through the reference and verify a subsequent find sees it.
  slot = 42;
  unsigned long const *p = m.find(2);
  REQUIRE(p != nullptr);
  REQUIRE(*p == 42);
}

TEST_CASE("dense_or_sparse_map::insert_or_zero returns a stable reference on re-insert",
          "[dense_or_sparse_map][dense][regression]") {
  // Regression for: the ++m.insert_or_zero(k) pattern requires that
  // re-inserting the same key returns the EXISTING slot's reference,
  // not a freshly-zeroed one. An earlier iteration of this file
  // returned slot-by-value for the on-hit branch and would have
  // silently dropped accumulated increments.
  dense_or_sparse_map<std::size_t, unsigned long, 4> m;
  ++m.insert_or_zero(1);
  ++m.insert_or_zero(1);
  ++m.insert_or_zero(1);
  unsigned long const *p = m.find(1);
  REQUIRE(p != nullptr);
  REQUIRE(*p == 3);
  REQUIRE(m.size() == 1);
}

TEST_CASE("dense_or_sparse_map::find distinguishes zero-valued entries from absent ones",
          "[dense_or_sparse_map][dense][regression]") {
  // Regression for the original dense_impl iteration where a value of 0
  // was used as "not inserted" sentinel. The `gaps` ref in
  // vcp_4_r_0/1.hpp is inserted at key 0 with initial value 0 and can
  // remain at 0 for the entire call (no disconnected v3 candidates).
  // The tail loop's find(0) must return a non-null pointer even when
  // the stored value is 0 — otherwise the unconnected-class correction
  // is silently skipped.
  dense_or_sparse_map<std::size_t, unsigned long, 4> m;
  unsigned long &slot = m.insert_or_zero(0);
  REQUIRE(slot == 0);
  unsigned long const *p = m.find(0);
  REQUIRE(p != nullptr);
  REQUIRE(*p == 0);
  // A never-inserted key correctly reports absent.
  REQUIRE(m.find(3) == nullptr);
}

TEST_CASE("dense_or_sparse_map::clear resets presence but reuses storage",
          "[dense_or_sparse_map][dense]") {
  dense_or_sparse_map<std::size_t, unsigned long, 4> m;
  m.insert_or_zero(1) = 10;
  m.insert_or_zero(3) = 30;
  REQUIRE(m.size() == 2);
  m.clear();
  REQUIRE(m.size() == 0);
  REQUIRE(m.find(1) == nullptr);
  REQUIRE(m.find(3) == nullptr);
  // Re-insert works on the cleared instance and starts from zero.
  unsigned long &slot = m.insert_or_zero(1);
  REQUIRE(slot == 0);
  ++slot;
  REQUIRE(*m.find(1) == 1);
}

TEST_CASE("dense_or_sparse_map::for_each walks populated entries", "[dense_or_sparse_map][dense]") {
  dense_or_sparse_map<std::size_t, unsigned long, 8> m;
  m.insert_or_zero(2) = 20;
  m.insert_or_zero(5) = 50;
  m.insert_or_zero(7) = 70;

  std::size_t key_sum = 0;
  unsigned long value_sum = 0;
  std::size_t count = 0;
  m.for_each([&](std::size_t k, unsigned long v) {
    key_sum += k;
    value_sum += v;
    ++count;
  });
  REQUIRE(count == 3);
  REQUIRE(key_sum == 2 + 5 + 7);
  REQUIRE(value_sum == 20 + 50 + 70);
}

TEST_CASE("dense_or_sparse_map::for_each does not emit entries from before clear",
          "[dense_or_sparse_map][dense][regression]") {
  dense_or_sparse_map<std::size_t, unsigned long, 8> m;
  m.insert_or_zero(2) = 20;
  m.insert_or_zero(5) = 50;
  m.clear();
  m.insert_or_zero(7) = 70;

  std::size_t count = 0;
  m.for_each([&](std::size_t /*k*/, unsigned long /*v*/) { ++count; });
  REQUIRE(count == 1);
}

TEST_CASE("dense_or_sparse_map exercises the sparse tier for key spaces over the budget",
          "[dense_or_sparse_map][sparse]") {
  // (1 << 21) keys at 8 bytes = 16 MB ⇒ sparse tier.
  constexpr std::size_t KEY_COUNT = std::size_t(1) << 21;
  dense_or_sparse_map<std::size_t, unsigned long, KEY_COUNT> m;
  static_assert(!(decltype(m)::is_dense));

  m.insert_or_zero(0) = 100;
  m.insert_or_zero(1 << 20) = 200;
  REQUIRE(m.size() == 2);
  REQUIRE(*m.find(0) == 100);
  REQUIRE(*m.find(1 << 20) == 200);
  REQUIRE(m.find(1 << 19) == nullptr);

  unsigned long sum = 0;
  m.for_each([&](std::size_t /*k*/, unsigned long v) { sum += v; });
  REQUIRE(sum == 300);
}

TEST_CASE("dense_or_sparse_map with pair_hash handles std::pair keys at the sparse tier",
          "[dense_or_sparse_map][sparse][pair]") {
  // KeyCount = 0 sentinel forces sparse tier; Key = std::pair.
  using pair_t = std::pair<std::size_t, std::size_t>;
  dense_or_sparse_map<pair_t, unsigned long, 0> m;
  static_assert(!(decltype(m)::is_dense));

  m.insert_or_zero(std::make_pair<std::size_t, std::size_t>(1, 2)) = 12;
  m.insert_or_zero(std::make_pair<std::size_t, std::size_t>(3, 4)) = 34;
  REQUIRE(*m.find(std::make_pair<std::size_t, std::size_t>(1, 2)) == 12);
  REQUIRE(*m.find(std::make_pair<std::size_t, std::size_t>(3, 4)) == 34);
  REQUIRE(m.find(std::make_pair<std::size_t, std::size_t>(2, 1)) == nullptr);
}

TEST_CASE("dense_or_sparse_map with pack_pair_key_by_bits handles std::pair keys at the "
          "dense tier",
          "[dense_or_sparse_map][dense][pair]") {
  // r = 3 ⇒ key space 2^6 = 64, well under the byte budget ⇒ dense.
  // pack_pair_key_by_bits<3> collapses (lo, hi) to (lo | (hi << 3)).
  using pair_t = std::pair<std::size_t, std::size_t>;
  constexpr std::size_t KEY_COUNT = std::size_t(1) << 6;
  dense_or_sparse_map<pair_t, unsigned long, KEY_COUNT, pack_pair_key_by_bits<3>> m;
  static_assert(decltype(m)::is_dense);

  m.insert_or_zero(std::make_pair<std::size_t, std::size_t>(1, 2)) = 12;
  m.insert_or_zero(std::make_pair<std::size_t, std::size_t>(3, 4)) = 34;
  // (1, 2) packs to 1 | (2 << 3) = 17.
  // (3, 4) packs to 3 | (4 << 3) = 35.
  REQUIRE(*m.find(std::make_pair<std::size_t, std::size_t>(1, 2)) == 12);
  REQUIRE(*m.find(std::make_pair<std::size_t, std::size_t>(3, 4)) == 34);
  // Verify for_each yields the original pair, not the packed index.
  bool saw_12 = false;
  bool saw_34 = false;
  m.for_each([&](pair_t const &key, unsigned long v) {
    if (key.first == 1 && key.second == 2) {
      REQUIRE(v == 12);
      saw_12 = true;
    }
    if (key.first == 3 && key.second == 4) {
      REQUIRE(v == 34);
      saw_34 = true;
    }
  });
  REQUIRE(saw_12);
  REQUIRE(saw_34);
}

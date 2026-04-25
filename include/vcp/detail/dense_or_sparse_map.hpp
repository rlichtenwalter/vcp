/*
Copyright (C) 2026 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with
the Vertex Collocation Profiles code base. If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef VCP_DETAIL_DENSE_OR_SPARSE_MAP_HPP
#define VCP_DETAIL_DENSE_OR_SPARSE_MAP_HPP

// `detail::dense_or_sparse_map` — a two-tier std::map-like container used
// internally by the vcp<4, r, d> hot path in place of the previous
// std::map usage for the two `edge_types`-family maps. The tier is
// chosen at compile time from KeyCount and sizeof(Value):
//
//   tier 1 (dense):   KeyCount * sizeof(Value) <= DENSE_BYTE_BUDGET
//                     backing = std::array<Value, KeyCount> (heap)
//                             + std::vector<std::size_t> present_keys
//
//   tier 2 (sparse):  otherwise
//                     backing = std::unordered_map<Key, Value>
//
// The single-tier `std::vector<std::pair>` linear-scan option in between
// was evaluated in benchmark/small_map_study/ and rejected: it only
// outperforms std::unordered_map at k ≤ ~8, and the dense tier already
// dominates that regime for any key space where dense is feasible. The
// middle tier had no operating point where it beat both neighbors, so
// it was dropped from the design. See
// benchmark/small_map_study/results/summary_phase_C.md for the
// microbenchmark numbers.
//
// Dense tier design:
//
//   - Value storage uses `new std::array<Value, KeyCount>` (default-
//     initialized, not value-initialized). For POD Value types this
//     leaves the array UNINITIALIZED at construction; a value is
//     zeroed lazily on its first `insert_or_zero`. Saves the O(N)
//     zero-init write that would otherwise dominate per-construction
//     cost at large KeyCount. The caller must never read a slot for
//     which `insert_or_zero` has not returned first (enforced by the
//     API — `find` returns nullptr for absent keys, and iteration
//     walks only the `present_keys` tracker).
//
//   - `present_keys` is a companion vector tracking which indices
//     have been populated in the current logical "session." This
//     distinguishes legitimately-zero entries from never-inserted
//     ones: the `gaps` ref used in vcp_4_r_*.hpp starts at 0 and can
//     stay 0 when no disconnected v3 candidates exist, and the tail
//     loop's `find(0)` must find it (for the unconnected-class
//     correction) not report absence.
//
//   - `clear()` clears present_keys in O(1) (vector::clear is O(k)
//     but the destructors are trivial for std::size_t). The values
//     array is NOT reset — its stale contents are hidden by the
//     presence tracker until re-populated.
//
//   - Since clear is cheap, a single `dense_or_sparse_map` instance
//     can be reused across many generate_vector calls as a class
//     member, avoiding the per-call malloc-path overhead that a
//     function-local std::array<V, N> at N ≈ 10^5 would incur.
//
// Sparse tier design:
//
//   - Uses std::unordered_map, not std::map. Phase C microbench
//     showed hash is 2-17× faster than tree across every tested
//     k ≥ 16, and the hot path has no iteration-order dependency
//     (the public-API counts map preserves ordering via a
//     conversion at the return boundary, not here).
//
//   - For `std::pair` keys (directed edge_types), a non-commutative
//     hash-combine is provided by `pair_hash`; stdlib does not
//     specialize std::hash for std::pair.
//
// Thread-safety: not thread-safe. Consistent with the pre-existing
// vcp<4, r, d> contract (v3Vertices, edge_types). See README.

#include <array>
#include <climits>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace vcp {
namespace detail {

// Byte budget for the dense tier. 8 MB fits comfortably in typical L2
// caches while covering r ≤ 20 (undirected, 8-byte Value) and r ≤ 10
// (directed, where the effective bit count is 2r). Typed as std::size_t
// at the literal boundary so downstream multiplications avoid an int
// overflow on 32-bit-int platforms.
inline constexpr std::size_t DENSE_BYTE_BUDGET = std::size_t(8) * 1024 * 1024;

// Compile-time predicate selecting the dense tier. KeyCount == 0 is
// treated as sparse (sentinel for "cannot statically size the key
// space").
template <typename Value, std::size_t KeyCount>
inline constexpr bool dense_fits = KeyCount > 0 && KeyCount <= DENSE_BYTE_BUDGET / sizeof(Value);

// --------------------------------------------------------------------
// Hash support for std::pair (directed edge_types key).
// --------------------------------------------------------------------

template <typename A, typename B> struct pair_hash {
  std::size_t operator()(std::pair<A, B> const &p) const noexcept {
    std::size_t const h1 = std::hash<A>()(p.first);
    std::size_t const h2 = std::hash<B>()(p.second);
    // Boost's hash_combine constant, width-selected: the 32-bit form
    // 0x9e3779b9 is what you see in most reference snippets, but on a
    // 64-bit size_t it leaves the upper 32 bits of h2 unmixed — a real
    // collision risk when the sparse tier fires (r > 20 undirected /
    // r > 10 directed), which is exactly the regime where key halves
    // span beyond 32 bits. The 64-bit constant 0x9e3779b97f4a7c15 is
    // the mixed-width analog used by Boost on 64-bit platforms.
    if constexpr (sizeof(std::size_t) == 8) {
      return h1 ^ (h2 + std::size_t(0x9e3779b97f4a7c15ULL) + (h1 << 6) + (h1 >> 2));
    } else {
      return h1 ^ (h2 + std::size_t(0x9e3779b9UL) + (h1 << 6) + (h1 >> 2));
    }
  }
};

template <typename Key> struct hash_for {
  using type = std::hash<Key>;
};

template <typename A, typename B> struct hash_for<std::pair<A, B>> {
  using type = pair_hash<A, B>;
};

// --------------------------------------------------------------------
// Default Key → dense-index projection.
// --------------------------------------------------------------------

// For integral Key types, the dense-tier index is just the key cast to
// std::size_t. The directed edge-types map uses std::pair<...> as its
// key type, which cannot be cast directly; that call site supplies its
// own KeyToIndex functor (`pack_pair_key_by_bits`) to collapse the
// pair into a single index via `lo | (hi << r)`.
template <typename Key> struct default_key_to_index {
  std::size_t operator()(Key const &k) const { return static_cast<std::size_t>(k); }
};

// Pack a std::pair<std::size_t, std::size_t> into a single std::size_t
// index by shifting the high half by `Bits`. Valid when 2*Bits fits
// in std::size_t. Used only by the dense tier; the sparse tier takes
// the raw pair as its hash key.
template <std::size_t Bits> struct pack_pair_key_by_bits {
  template <typename A, typename B> std::size_t operator()(std::pair<A, B> const &p) const {
    return static_cast<std::size_t>(p.first) | (static_cast<std::size_t>(p.second) << Bits);
  }
};

// --------------------------------------------------------------------
// Tier 1: dense.
// --------------------------------------------------------------------

template <typename Key, typename Value, std::size_t KeyCount> class dense_impl {
  static_assert(KeyCount > 0, "dense tier requires KeyCount > 0");

public:
  dense_impl()
      : values_(std::unique_ptr<std::array<Value, KeyCount>>(new std::array<Value, KeyCount>)) {
    // NB: `new std::array<Value, KeyCount>` (no parens) is
    // default-initialization; for POD Value types the elements are
    // UNINITIALIZED. insert_or_zero lazily zeroes on first touch.
  }

  // Clear the logical contents without touching the values array.
  // After clear(), find() returns nullptr for every key and for_each
  // yields nothing, until insert_or_zero repopulates.
  void clear() { present_keys_.clear(); }

  // Insert a zero if (index, original_key) is not present, returning a
  // reference to the slot either way. The `index` is the packed
  // std::size_t that the dense array uses as its direct index;
  // `original_key` is the caller's Key type (possibly std::pair<...>),
  // saved for for_each callbacks so the caller does not need to
  // re-derive it from the index.
  Value &insert_or_zero(std::size_t index, Key const &original_key) {
    for (std::size_t i = 0; i < present_keys_.size(); ++i) {
      if (present_keys_[i].first == index) {
        return (*values_)[index];
      }
    }
    present_keys_.emplace_back(index, original_key);
    (*values_)[index] = Value{}; // lazy zero-init on first touch
    return (*values_)[index];
  }

  // Return a pointer to the stored value, or nullptr if `index` was not
  // previously inserted. Read-only; does not populate.
  Value const *find(std::size_t index) const {
    for (auto const &p : present_keys_) {
      if (p.first == index) {
        return &(*values_)[index];
      }
    }
    return nullptr;
  }

  std::size_t size() const { return present_keys_.size(); }

  // Invoke fn(key, value) for each populated entry. Iteration order is
  // insertion order — NOT the sorted order std::map would provide.
  // The hot path has no iteration-order dependency (verified by the
  // regression byte-diff suite that must stay green after the
  // std::map → dense_or_sparse_map swap).
  template <typename Fn> void for_each(Fn fn) const {
    for (auto const &p : present_keys_) {
      fn(p.second, (*values_)[p.first]);
    }
  }

private:
  std::unique_ptr<std::array<Value, KeyCount>> values_;
  // Each entry is (packed_index_into_values, original_key). The
  // original key is stored so for_each yields the caller's Key type
  // without needing an inverse of the KeyToIndex projection. For
  // integral Key (undirected edge types), this is a small amount of
  // duplicate storage (size_t + size_t); for std::pair<size_t, size_t>
  // (directed), the "index" collapses the pair under a packing
  // function while the full pair stays available for iteration
  // callbacks.
  std::vector<std::pair<std::size_t, Key>> present_keys_;
};

// --------------------------------------------------------------------
// Tier 2: sparse (std::unordered_map).
// --------------------------------------------------------------------

template <typename Key, typename Value> class sparse_impl {
public:
  using hash_type = typename hash_for<Key>::type;

  void clear() { map_.clear(); }

  Value &insert_or_zero(Key const &key) {
    return map_.insert(std::make_pair(key, Value{})).first->second;
  }

  Value const *find(Key const &key) const {
    auto it = map_.find(key);
    return it == map_.end() ? nullptr : &it->second;
  }

  std::size_t size() const { return map_.size(); }

  template <typename Fn> void for_each(Fn fn) const {
    for (auto const &kv : map_) {
      fn(kv.first, kv.second);
    }
  }

private:
  std::unordered_map<Key, Value, hash_type> map_;
};

// --------------------------------------------------------------------
// Facade.
// --------------------------------------------------------------------

template <typename Key, typename Value, std::size_t KeyCount,
          typename KeyToIndex = default_key_to_index<Key>>
class dense_or_sparse_map {
public:
  // Exposed so tests can assert which tier was selected without
  // inspecting sizeof()-derived types.
  static constexpr bool is_dense = dense_fits<Value, KeyCount>;

  void clear() { impl_.clear(); }

  Value &insert_or_zero(Key const &key) {
    if constexpr (is_dense) {
      return impl_.insert_or_zero(KeyToIndex{}(key), key);
    } else {
      return impl_.insert_or_zero(key);
    }
  }

  Value const *find(Key const &key) const {
    if constexpr (is_dense) {
      return impl_.find(KeyToIndex{}(key));
    } else {
      return impl_.find(key);
    }
  }

  std::size_t size() const { return impl_.size(); }

  template <typename Fn> void for_each(Fn fn) const { impl_.for_each(fn); }

private:
  using impl_type =
      std::conditional_t<is_dense, dense_impl<Key, Value, KeyCount>, sparse_impl<Key, Value>>;
  impl_type impl_;
};

// Compile-time key-count helpers. These drive the KeyCount template
// parameter for the per-specialization instantiations.
//
// Undirected edge key: r bits → 2^r distinct values.
// Directed edge key (as a packed std::size_t index into the dense
//   tier): 2r bits → 2^(2r) distinct values.
//
// When the bit count exceeds std::size_t's width, the key count is
// sentinel'd to 0, which forces the sparse tier via dense_fits's
// KeyCount==0 check. This keeps the template well-formed at r=30
// directed (60 bits, larger than size_t on 32-bit platforms).
template <std::size_t bits> inline constexpr std::size_t key_count_from_bits() {
  return bits >= sizeof(std::size_t) * CHAR_BIT ? std::size_t(0) : (std::size_t(1) << bits);
}

} // namespace detail
} // namespace vcp

#endif

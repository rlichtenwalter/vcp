// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_DYNAMIC_MAPPER_HPP
#define VCP_VCP_DYNAMIC_MAPPER_HPP

#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_int.hpp>
#include <climits>
#include <cstddef>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>

namespace vcp {

template <std::size_t base, std::size_t exponent> struct TMP_power {
  enum { value_matrix = base * TMP_power<base, exponent - 1>::value_matrix };
};

template <std::size_t base> struct TMP_power<base, 0> {
  enum { value_matrix = 1 };
};

template <std::size_t value1, std::size_t value2> struct TMP_min {
  enum { value_matrix = value1 < value2 ? value1 : value2 };
};

template <std::size_t value1, std::size_t value2> struct TMP_max {
  enum { value_matrix = value1 > value2 ? value1 : value2 };
};

/**
 * @brief Compute-on-demand mapper from connectivity matrices to canonical element addresses.
 *
 * Maps n-vertex subgraph connectivity matrices to their canonical (isomorphism-reduced)
 * element address by iterating over all permutations of the n-2 non-pivot vertices
 * and taking the minimum encoded address. This avoids the O(2^k) precomputation and
 * memory cost of `vcp_static_mapper`, at the expense of per-call computation.
 *
 * The subgraph address is a packed integer: each undirected (or directed) pair of
 * vertices occupies `r` bits, laid out in the order assigned by the internal
 * `value_matrix`. The canonical address is the minimum across all isomorphic
 * representations (i.e., permutations of vertices 2..n-1; vertices 0 and 1 are
 * the pivot pair and are never permuted).
 *
 * This class is the backbone of the generic `vcp<n, r, d>` template and the
 * multirelational partial specializations `vcp<3, r, d>` and `vcp<4, r, d>`.
 *
 * @tparam n Number of vertices in each subgraph.
 * @tparam r Number of edge relations (bits per edge pair).
 * @tparam d True for directed graphs; false for undirected.
 */
template <std::size_t n, std::size_t r, bool d> class vcp_dynamic_mapper {
public:
  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   *
   * Alias of `multirelational_graph<r>::connectivity_address_type`.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete n-vertex subgraph.
   *
   * Wide enough to hold all n*(n-1)*r*(d+1)/2 bits of the subgraph encoding.
   * Resolves to `std::size_t` when the bit count fits within a machine word,
   * or to a fixed-width Boost.Multiprecision unsigned integer otherwise.
   */
  using subgraph_address_type = typename std::conditional<
      n *(n - 1) * r *(d + 1) / 2 <= CHAR_BIT * sizeof(std::size_t), std::size_t,
      boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
          n *(n - 1) * r *(d + 1) / 2, n *(n - 1) * r *(d + 1) / 2,
          boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>>::
      type;

  /** @brief Construct the mapper and initialize the internal bit-position table. */
  vcp_dynamic_mapper();

  /**
   * @brief Return the total number of distinct subgraph addresses.
   *
   * Equal to 2^(n*(n-1)*r*(d+1)/2): the number of ways to assign r-bit
   * bitmasks to all vertex pairs.
   *
   * @return The number of distinct subgraph addresses.
   */
  constexpr subgraph_address_type subgraph_count() const;

  /**
   * @brief Compute the subgraph address for the given connectivity matrix.
   *
   * Encodes the connectivity matrix into a packed integer by summing
   * `connectivity(i, j) << value_matrix(i, j)` over all relevant pairs.
   * The result is not canonical; isomorphic subgraphs may yield different
   * addresses.
   *
   * @param connectivity n × n matrix of r-bit edge bitmasks.
   * @return Packed integer encoding of the subgraph.
   */
  subgraph_address_type
  subgraph_address(square_matrix<connectivity_address_type, n> const &connectivity) const;

  /**
   * @brief Compute the canonical element address for the given connectivity matrix.
   *
   * Returns the minimum subgraph address across all permutations of the
   * non-pivot vertices (vertices 2 through n-1). Isomorphic subgraphs
   * always yield the same canonical address, which corresponds to the VCP
   * element index as defined in Lichtenwalter & Chawla (2012, 2014).
   *
   * @param connectivity n × n matrix of r-bit edge bitmasks.
   * @return Canonical element address (minimum over isomorphic permutations).
   */
  subgraph_address_type
  canonical_subgraph_address(square_matrix<connectivity_address_type, n> const &connectivity) const;

  /**
   * @brief Reconstruct the connectivity matrix for the given element address.
   *
   * Decodes the packed address back into an n × n connectivity matrix. The
   * returned matrix corresponds to one representative of the isomorphism class;
   * it is the same representative that was used to compute the canonical address.
   *
   * @param address A canonical element address in [0, subgraph_count()).
   * @return n × n connectivity matrix reconstructed from the address.
   */
  square_matrix<connectivity_address_type, n>
  element_structure(subgraph_address_type const &address) const;

private:
  square_matrix<std::size_t, n> value_matrix;
  subgraph_address_type
  canonical_subgraph_address(square_matrix<connectivity_address_type, n> const &connectivity,
                             subgraph_address_type subgraph_address) const;
};

template <std::size_t n, std::size_t r, bool d> vcp_dynamic_mapper<n, r, d>::vcp_dynamic_mapper() {
  std::size_t index(0);
  for (std::size_t row(0); row < n; ++row) {
    value_matrix(row, row) = 0;
    for (std::size_t column(row + 1); column < n; ++column) {
      value_matrix(row, column) = r * index++;
      value_matrix(column, row) = d ? (r * index++) : value_matrix(row, column);
    }
  }
}

template <std::size_t n, std::size_t r, bool d>
constexpr typename vcp_dynamic_mapper<n, r, d>::subgraph_address_type
vcp_dynamic_mapper<n, r, d>::subgraph_count() const {
  return subgraph_address_type(1) << n * (n - 1) * r * (d + 1) / 2;
}

template <std::size_t n, std::size_t r, bool d>
typename vcp_dynamic_mapper<n, r, d>::subgraph_address_type
vcp_dynamic_mapper<n, r, d>::subgraph_address(
    square_matrix<connectivity_address_type, n> const &connectivity) const {
  subgraph_address_type subgraph_address(0);
  for (std::size_t row(0); row < n; ++row) {
    for (std::size_t column(d ? 0 : row + 1); column < n; ++column) {
      if (row != column) {
        subgraph_address += subgraph_address_type(connectivity(row, column))
                            << value_matrix(row, column);
      }
    }
  }
  return subgraph_address;
}

template <std::size_t n, std::size_t r, bool d>
typename vcp_dynamic_mapper<n, r, d>::subgraph_address_type
vcp_dynamic_mapper<n, r, d>::canonical_subgraph_address(
    square_matrix<connectivity_address_type, n> const &connectivity) const {
  return canonical_subgraph_address(connectivity, subgraph_address(connectivity));
}

template <std::size_t n, std::size_t r, bool d>
typename vcp_dynamic_mapper<n, r, d>::subgraph_address_type
vcp_dynamic_mapper<n, r, d>::canonical_subgraph_address(
    square_matrix<connectivity_address_type, n> const &connectivity,
    subgraph_address_type subgraph_address) const {
  std::array<std::size_t, n> permuter;
  for (std::size_t row(0); row < n; ++row) {
    permuter[row] = row;
  }
  // The inner-loop bounds match those of `subgraph_address` above so the
  // baseline and every permuted variant iterate the same set of slots.
  // In particular, for undirected (d == false) we must iterate the upper
  // triangle only — otherwise the permuted term double-counts each pair
  // (value_matrix is symmetric for undirected, so it adds the same slot
  // contribution twice), the permuted address always exceeds the baseline,
  // `std::min` always picks the baseline, and canonicalization silently
  // degrades to a no-op for any caller that populates both triangles of
  // the connectivity matrix.
  while (std::next_permutation(permuter.begin() + 2, permuter.end())) {
    subgraph_address_type isomorphism_address(0);
    for (std::size_t row(0); row < n; ++row) {
      for (std::size_t column(d ? 0 : row + 1); column < n; ++column) {
        if (row != column) {
          isomorphism_address += subgraph_address_type(connectivity(row, column))
                                 << value_matrix(permuter[row], permuter[column]);
        }
      }
    }
    subgraph_address = std::min(subgraph_address, isomorphism_address);
  }
  return subgraph_address;
}

template <std::size_t n, std::size_t r, bool d>
square_matrix<typename vcp_dynamic_mapper<n, r, d>::connectivity_address_type, n>
vcp_dynamic_mapper<n, r, d>::element_structure(subgraph_address_type const &address) const {
  // Each (row, col) slot occupies r bits of `address`, so the slot-mask is
  // (1 << r) - 1 and advancing to the next slot shifts by r. The previous
  // form shifted by (1 << r) — the size of the r-relation power-set —
  // which scrambled decoding for every r (and was all-zeros UB for r with
  // 2**r >= the address bit-width).
  subgraph_address_type const slot_mask((vcp_dynamic_mapper::subgraph_address_type(1) << r) - 1);
  subgraph_address_type remaining(address);
  square_matrix<connectivity_address_type, n> matrix;
  std::size_t row(0);
  std::size_t column(1);
  while (remaining > 0) {
    matrix(row, column) = remaining & slot_mask;
    remaining >>= r;
    if (d) {
      std::swap(row, column);
      if (row < column) {
        ++column;
        if (column >= n) {
          ++row;
          column = row + 1;
        }
      }
    } else {
      ++column;
      if (column >= n) {
        ++row;
        column = row + 1;
      }
    }
  }
  return matrix;
}

} // namespace vcp

#endif

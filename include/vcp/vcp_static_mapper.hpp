// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_STATIC_MAPPER_HPP
#define VCP_VCP_STATIC_MAPPER_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ostream>
#include <ranges>
#include <stdexcept>
#include <vcp/square_matrix.hpp>
#include <vector>

namespace vcp {

/**
 * @brief Precomputed lookup-table mapper from subgraph addresses to canonical element addresses.
 *
 * At construction time, all 2^(n*(n-1)*r*(d+1)/2) subgraph encodings are enumerated
 * and grouped by isomorphism class. Each class is assigned a dense element address
 * starting at 0. The result is a flat `map` vector that converts any subgraph address
 * to its canonical element address in O(1).
 *
 * This is faster than `vcp_dynamic_mapper` for repeated lookups, but requires
 * O(subgraph_count) memory at construction time. For large `n`, `r`, or `d`
 * the table size may be prohibitive; `vcp_dynamic_mapper` is the alternative.
 *
 * Parameters `n`, `r`, `d` are runtime values, which is why this class uses
 * dynamic (`square_matrix<std::size_t>`) rather than fixed-size matrices.
 */
class vcp_static_mapper {
public:
  /**
   * @brief Return the total number of distinct subgraph addresses for the given parameters.
   *
   * Equal to 2^(n*(n-1)*r*(d+1)/2). This is the size of the internal lookup table.
   *
   * @param n Number of vertices in each subgraph.
   * @param r Number of edge relations.
   * @param d True for directed; false for undirected.
   * @return Total number of subgraph addresses.
   */
  [[nodiscard]] static std::size_t subgraph_count(std::size_t n, std::size_t r, bool d);

  /**
   * @brief Construct the mapper and build the full isomorphism lookup table.
   *
   * Enumerates all subgraph encodings and assigns canonical element addresses.
   * Construction time is O(subgraph_count * n!) due to the permutation search.
   *
   * @param n Number of vertices in each subgraph.
   * @param r Number of edge relations.
   * @param d True for directed; false for undirected.
   */
  vcp_static_mapper(std::size_t n, std::size_t r, bool d);

  /** @brief Return the subgraph vertex count. */
  [[nodiscard]] std::size_t n() const noexcept;

  /** @brief Return the number of edge relations. */
  [[nodiscard]] std::size_t r() const noexcept;

  /** @brief Return true if the mapper was constructed for directed graphs. */
  [[nodiscard]] bool d() const noexcept;

  /**
   * @brief Compute the subgraph address for the given fixed-size connectivity matrix.
   *
   * @tparam n Side length of the connectivity matrix.
   * @param connectivity n × n matrix of edge relation values.
   * @return Packed subgraph address integer.
   */
  template <std::size_t n>
  [[nodiscard]] std::size_t
  subgraph_address(square_matrix<std::size_t, n> const &connectivity) const;

  /**
   * @brief Compute the canonical element address for the given fixed-size connectivity matrix.
   *
   * @tparam n Side length of the connectivity matrix.
   * @param connectivity n × n matrix of edge relation values.
   * @return Canonical element address in [0, element_count).
   */
  template <std::size_t n>
  [[nodiscard]] std::size_t
  element_address(square_matrix<std::size_t, n> const &connectivity) const;

  /**
   * @brief Compute the subgraph address for the given dynamic connectivity matrix.
   *
   * @param connectivity Dynamic square matrix of edge relation values.
   * @return Packed subgraph address integer.
   */
  [[nodiscard]] std::size_t subgraph_address(square_matrix<std::size_t> const &connectivity) const;

  /**
   * @brief Compute the canonical element address for the given dynamic connectivity matrix.
   *
   * @param connectivity Dynamic square matrix of edge relation values.
   * @return Canonical element address in [0, element_count).
   */
  [[nodiscard]] std::size_t element_address(square_matrix<std::size_t> const &connectivity) const;

  /**
   * @brief Look up the canonical element address from a precomputed subgraph address.
   *
   * O(1) table lookup; the index must be in [0, subgraph_count()).
   *
   * @param subgraph_address Packed subgraph address previously computed by subgraph_address().
   * @return Canonical element address.
   */
  [[nodiscard]] std::size_t element_address(std::size_t subgraph_address) const;

  /**
   * @brief Reconstruct the connectivity matrix for the given element address.
   *
   * Decodes the element address back into a dynamic square matrix. The
   * returned matrix is one canonical representative of the isomorphism class.
   *
   * @param element_address Canonical element address in [0, element_count).
   * @return Connectivity matrix reconstructed from the element address.
   */
  [[nodiscard]] square_matrix<std::size_t> element_structure(std::size_t element_address) const;

  friend std::ostream &operator<<(std::ostream &os, vcp_static_mapper const &mapper);

private:
  std::size_t n_;
  std::size_t r_;
  bool d_;
  std::vector<std::size_t> map;
  square_matrix<std::size_t> value_matrix;
};

inline std::size_t vcp_static_mapper::subgraph_count(std::size_t n, std::size_t r, bool d) {
  std::size_t const bits = n * (n - 1) * r * (d + 1) / 2;
  if (bits >= std::numeric_limits<std::size_t>::digits) {
    throw std::length_error("vcp_static_mapper: parameter combination overflows std::size_t");
  }
  return std::size_t(1) << bits;
}

inline vcp_static_mapper::vcp_static_mapper(std::size_t n, std::size_t r, bool d)
    : n_(n), r_(r), d_(d), map(vcp_static_mapper::subgraph_count(n, r, d)), value_matrix(n_) {
  std::size_t const r_pset = std::size_t(1) << r_;
  std::size_t index(0);
  for (std::size_t row(0); row < n_; ++row) {
    value_matrix(row, row) = 0;
    for (std::size_t column(row + 1); column < n_; ++column) {
      value_matrix(row, column) = std::size_t(1) << (r_ * index++);
      value_matrix(column, row) =
          d_ ? (std::size_t(1) << (r_ * index++)) : value_matrix(row, column);
    }
  }

  square_matrix<std::size_t> connectivity(n_);
  std::vector<std::size_t> permuter(n_);

  std::size_t element_address(0);
  while (!connectivity(0, 0)) {
    std::size_t subgraph_address(0);
    for (std::size_t row(0); row < n_; ++row) {
      for (std::size_t column(d_ ? 0 : row + 1); column < n_; ++column) {
        subgraph_address += connectivity(row, column) * value_matrix(row, column);
      }
    }
    if (map[subgraph_address] ==
        0) { // we have not encountered this subgraph or its isomorphic equivalents before
      map[subgraph_address] = element_address;
      for (std::size_t i(0); i < n_; ++i) { // reset the permuter vector
        permuter[i] = i;
      }
      while (std::next_permutation(permuter.begin() + 2, permuter.end())) { // find all isomorphisms
        std::size_t isomorphism_address(0);
        for (std::size_t row(0); row < n_; ++row) {
          for (std::size_t column(0); column < n_; ++column) {
            isomorphism_address +=
                connectivity(row, column) * value_matrix(permuter[row], permuter[column]);
          }
        }
        map[isomorphism_address] = element_address;
      }
      ++element_address;
    }

    // move to the next subgraph
    std::size_t row(0);
    std::size_t column(1);
    std::size_t *cell(&connectivity(row, column));
    ++(*cell);
    while (*cell == r_pset) {
      *cell = 0;
      if (d_) {
        std::swap(row, column);
        if (row < column) {
          ++column;
          if (column >= n_) {
            ++row;
            column = row + 1;
            if (row == n_ - 1) {
              row = 0;
              column = 0;
            }
          }
        }
      } else {
        ++column;
        if (column >= n_) {
          ++row;
          column = row + 1;
          if (row == n_ - 1) {
            row = 0;
            column = 0;
          }
        }
      }
      cell = &connectivity(row, column);
      ++(*cell);
    }
  }
}

template <std::size_t nn>
std::size_t
vcp_static_mapper::subgraph_address(square_matrix<std::size_t, nn> const &connectivity) const {
  std::size_t address(0);
  for (std::size_t row(0); row < nn; ++row) {
    for (std::size_t column(d_ ? 0 : row + 1); column < nn; ++column) {
      address += connectivity(row, column) * value_matrix(row, column);
    }
  }
  return address;
}

template <std::size_t nn>
std::size_t
vcp_static_mapper::element_address(square_matrix<std::size_t, nn> const &connectivity) const {
  return map[subgraph_address(connectivity)];
}

std::size_t inline vcp_static_mapper::subgraph_address(
    square_matrix<std::size_t> const &connectivity) const {
  std::size_t address(0);
  for (std::size_t row(0); row < n_; ++row) {
    for (std::size_t column(d_ ? 0 : row + 1); column < n_; ++column) {
      address += connectivity(row, column) * value_matrix(row, column);
    }
  }
  return address;
}

std::size_t inline vcp_static_mapper::element_address(
    square_matrix<std::size_t> const &connectivity) const {
  return map[subgraph_address(connectivity)];
}

inline square_matrix<std::size_t> vcp_static_mapper::element_structure(std::size_t address) const {
  std::size_t const r_pset = std::size_t(1) << r_;
  square_matrix<std::size_t> matrix(n_);
  std::size_t i = 0;
  std::size_t j = 1;
  while (address > 0) {
    matrix(i, j) = address % r_pset;
    address /= r_pset;
    if (d_) {
      std::swap(i, j);
      if (i < j) {
        ++j;
        if (j >= n_) {
          ++i;
          j = i + 1;
        }
      }
    } else {
      ++j;
      if (j >= n_) {
        ++i;
        j = i + 1;
      }
    }
  }
  return matrix;
}

/**
 * @brief Write the full element-address lookup table to an output stream.
 *
 * Emits one element address per line in subgraph-address order, from 0 to
 * subgraph_count() - 1. Useful for debugging and for cross-checking against
 * independently computed VCP element mappings.
 *
 * @param os     Output stream.
 * @param mapper Mapper whose table is to be written.
 * @return Reference to @p os.
 */
inline std::ostream &operator<<(std::ostream &os, vcp_static_mapper const &mapper) {
  std::ranges::copy(mapper.map, std::ostream_iterator<std::size_t>(os, "\n"));
  return os;
}

} // namespace vcp

#endif

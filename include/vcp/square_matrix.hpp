// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_SQUARE_MATRIX_HPP
#define VCP_SQUARE_MATRIX_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

namespace vcp {

/**
 * @brief Fixed-size square matrix stored in row-major order on the stack.
 *
 * Wraps a `std::array` of `n * n` elements value-initialized to zero.
 * The compile-time dimension `n` must be non-zero; use the `n = 0`
 * partial specialization for runtime-sized matrices. VCP algorithms use
 * this type to represent the connectivity between the `n` pivot vertices
 * of a candidate subgraph.
 *
 * @tparam value_type Element type stored in each cell.
 * @tparam n          Side length of the matrix; `n * n` elements are stored.
 */
template <typename value_type, std::size_t n = 0> class square_matrix {
public:
  /** @brief Return the side length of the matrix. */
  std::size_t size() const;

  /**
   * @brief Return a mutable reference to the element at (row, column).
   *
   * @param row    Zero-based row index in [0, size()).
   * @param column Zero-based column index in [0, size()).
   * @return Mutable reference to the element.
   */
  value_type &operator()(std::size_t row, std::size_t column);

  /**
   * @brief Return a const reference to the element at (row, column).
   *
   * @param row    Zero-based row index in [0, size()).
   * @param column Zero-based column index in [0, size()).
   * @return Const reference to the element.
   */
  value_type const &operator()(std::size_t row, std::size_t column) const;

  template <typename value_type_, std::size_t n_>
  friend std::ostream &operator<<(std::ostream &os, square_matrix<value_type_, n_> const &matrix);

private:
  std::array<value_type, n * n> data = {{0}};
};

/**
 * @brief Write a fixed-size square matrix to an output stream.
 *
 * Outputs rows separated by newlines; elements within each row are
 * separated by commas.
 *
 * @tparam value_type Element type.
 * @tparam n          Side length.
 * @param os     Destination output stream.
 * @param matrix Matrix to write.
 * @return Reference to @p os.
 */
template <typename value_type, size_t n>
std::ostream &operator<<(std::ostream &os, square_matrix<value_type, n> const &matrix);

template <typename value_type, std::size_t n>
std::size_t square_matrix<value_type, n>::size() const {
  return std::sqrt(data.size());
}

template <typename value_type, std::size_t n>
value_type &square_matrix<value_type, n>::operator()(std::size_t row, std::size_t column) {
  return const_cast<value_type &>(
      static_cast<square_matrix<value_type, n> const &>(*this)(row, column));
}

template <typename value_type, std::size_t n>
value_type const &square_matrix<value_type, n>::operator()(std::size_t row,
                                                           std::size_t column) const {
  return data[size() * row + column];
}

template <typename value_type, std::size_t n>
std::ostream &operator<<(std::ostream &os, square_matrix<value_type, n> const &matrix) {
  for (std::size_t row(0); row < matrix.size(); ++row) {
    for (std::size_t column(0); column < matrix.size(); ++column) {
      os << matrix(row, column);
      if (column < matrix.size() - 1) {
        os << ',';
      }
    }
    os << '\n';
  }
  return os;
}

template <typename value_type>
std::ostream &operator<<(std::ostream &os, square_matrix<value_type, 0> const &matrix);

/**
 * @brief Dynamic square matrix whose side length is set at construction time.
 *
 * This partial specialization (n = 0) stores elements in a `std::vector`
 * rather than a `std::array`, allowing the size to be chosen at runtime.
 * The interface is intentionally identical to the fixed-size primary template
 * so that code can switch between the two without modification. The VCP
 * static-mapper and tool code use this variant when the graph parameters
 * `n`, `r`, `d` are not known at compile time.
 *
 * @tparam value_type Element type stored in each cell.
 */
template <typename value_type> class square_matrix<value_type, 0> {
public:
  /** @brief Return the side length of the matrix. */
  std::size_t size() const;

  /**
   * @brief Resize the matrix to @p n × @p n, zero-initializing new elements.
   *
   * Existing element values are not preserved across a resize.
   *
   * @param n New side length.
   */
  void resize(std::size_t n);

  /** @brief Construct an empty (0 × 0) matrix. */
  square_matrix() = default;

  /**
   * @brief Construct an n × n matrix with all elements zero-initialized.
   *
   * @param n Side length of the matrix.
   */
  square_matrix(std::size_t n);

  /**
   * @brief Construct a dynamic matrix from a fixed-size matrix by copying its elements.
   *
   * @tparam n Side length of the source fixed-size matrix.
   * @param matrix Source matrix to copy.
   */
  template <std::size_t n> square_matrix(square_matrix<value_type, n> const &matrix);

  /**
   * @brief Assign from a fixed-size matrix, resizing and copying elements.
   *
   * @tparam n Side length of the source fixed-size matrix.
   * @param matrix Source matrix to copy.
   * @return Reference to this matrix.
   */
  template <std::size_t n>
  square_matrix<value_type, 0> &operator=(square_matrix<value_type, n> const &matrix);

  /**
   * @brief Return a mutable reference to the element at (row, column).
   *
   * @param row    Zero-based row index in [0, size()).
   * @param column Zero-based column index in [0, size()).
   * @return Mutable reference to the element.
   */
  value_type &operator()(std::size_t row, std::size_t column);

  /**
   * @brief Return a const reference to the element at (row, column).
   *
   * @param row    Zero-based row index in [0, size()).
   * @param column Zero-based column index in [0, size()).
   * @return Const reference to the element.
   */
  value_type const &operator()(std::size_t row, std::size_t column) const;

private:
  std::vector<value_type> data;
};

template <typename value_type> std::size_t square_matrix<value_type, 0>::size() const {
  return std::sqrt(data.size());
}

template <typename value_type> void square_matrix<value_type, 0>::resize(std::size_t n) {
  data.resize(n * n);
}

template <typename value_type>
square_matrix<value_type, 0>::square_matrix(std::size_t n) : data(n * n, 0) {}

template <typename value_type>
template <std::size_t n>
square_matrix<value_type, 0>::square_matrix(square_matrix<value_type, n> const &matrix)
    : data(n * n) {
  std::copy(&matrix(0, 0), &matrix(0, 0) + n * n, &data[0]);
}

template <typename value_type>
template <std::size_t n>
square_matrix<value_type, 0> &
square_matrix<value_type, 0>::operator=(square_matrix<value_type, n> const &matrix) {
  if (static_cast<void const *>(this) != static_cast<void const *>(&matrix)) {
    data.resize(n * n);
    std::copy(&matrix(0, 0), &matrix(0, 0) + n * n, &data[0]);
  }
  return *this;
}

template <typename value_type>
value_type &square_matrix<value_type, 0>::operator()(std::size_t row, std::size_t column) {
  return const_cast<value_type &>(
      static_cast<square_matrix<value_type, 0> const &>(*this)(row, column));
}

template <typename value_type>
value_type const &square_matrix<value_type, 0>::operator()(std::size_t row,
                                                           std::size_t column) const {
  return data[size() * row + column];
}

/**
 * @brief Write a dynamic square matrix to an output stream.
 *
 * Outputs rows separated by newlines; elements within each row are
 * separated by commas.
 *
 * @tparam value_type Element type.
 * @param os     Destination output stream.
 * @param matrix Matrix to write.
 * @return Reference to @p os.
 */
template <typename value_type>
std::ostream &operator<<(std::ostream &os, square_matrix<value_type, 0> const &matrix) {
  for (std::size_t row(0); row < matrix.size(); ++row) {
    for (std::size_t column(0); column < matrix.size(); ++column) {
      os << matrix(row, column);
      if (column < matrix.size() - 1) {
        os << ',';
      }
    }
    os << '\n';
  }
  return os;
}

} // namespace vcp

#endif

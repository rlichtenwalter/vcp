// Unit tests for vcp::square_matrix — the flat-storage square matrix
// template used as a lookup table by the VCP mappers and as a connectivity
// matrix by the vcp<4, r, *> specializations.
//
// Coverage focuses on (a) value round-trip through operator()(row, col),
// (b) size() after construction and resize, and (c) the fixed-to-dynamic
// conversion paths (copy constructor and operator=), which had two
// latent bugs fixed alongside this suite:
//   - copy ctor std::copy used matrix(n-1, n-1) (the last element's
//     address) as the end iterator, dropping the last element.
//   - operator= resize(n) instead of resize(n*n), producing a heap
//     buffer overflow for n > 1.
// Both conversion paths were latent (unreferenced by current callers),
// but the bugs were real. These tests lock in the correct behavior.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

#include <vcp/square_matrix.hpp>

TEST_CASE("fixed-size square_matrix value-initializes to zero", "[square_matrix]") {
  vcp::square_matrix<int, 3> m;
  REQUIRE(m.size() == 3);
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      REQUIRE(m(row, col) == 0);
    }
  }
}

TEST_CASE("fixed-size square_matrix write/read round-trip", "[square_matrix]") {
  vcp::square_matrix<int, 4> m;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      m(row, col) = static_cast<int>(10 * row + col);
    }
  }
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      REQUIRE(m(row, col) == static_cast<int>(10 * row + col));
    }
  }
}

TEST_CASE("fixed-size square_matrix corner cells reachable", "[square_matrix]") {
  // Paranoia check: the four corners, including the (n-1, n-1) cell
  // that the buggy conversion ctor dropped.
  vcp::square_matrix<int, 4> m;
  m(0, 0) = 1;
  m(0, 3) = 2;
  m(3, 0) = 3;
  m(3, 3) = 4;
  REQUIRE(m(0, 0) == 1);
  REQUIRE(m(0, 3) == 2);
  REQUIRE(m(3, 0) == 3);
  REQUIRE(m(3, 3) == 4);
}

TEST_CASE("dynamic square_matrix default-constructs empty", "[square_matrix]") {
  vcp::square_matrix<int, 0> m;
  REQUIRE(m.size() == 0);
}

TEST_CASE("dynamic square_matrix sized constructor zero-initializes", "[square_matrix]") {
  vcp::square_matrix<int, 0> m(3);
  REQUIRE(m.size() == 3);
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      REQUIRE(m(row, col) == 0);
    }
  }
}

TEST_CASE("dynamic square_matrix resize changes size", "[square_matrix]") {
  // resize(n) takes the DIMENSION and re-sizes internal storage to n*n,
  // matching the sized constructor's semantics.
  vcp::square_matrix<int, 0> m;
  REQUIRE(m.size() == 0);
  m.resize(5);
  REQUIRE(m.size() == 5);
  m(4, 4) = 42;
  REQUIRE(m(4, 4) == 42);
  m.resize(2);
  REQUIRE(m.size() == 2);
}

TEST_CASE("dynamic square_matrix write/read round-trip", "[square_matrix]") {
  vcp::square_matrix<int, 0> m(4);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      m(row, col) = static_cast<int>(100 * row + col);
    }
  }
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      REQUIRE(m(row, col) == static_cast<int>(100 * row + col));
    }
  }
}

TEST_CASE("fixed-to-dynamic conversion ctor copies every cell", "[square_matrix][conversion]") {
  // Regression for: std::copy(&m(0,0), &m(n-1,n-1), &data[0]) — the end
  // iterator pointed at the last element instead of one-past, dropping
  // the (n-1, n-1) cell. This test writes distinct nonzero values to
  // every cell (including the diagonal) and requires that every one
  // survive the conversion.
  vcp::square_matrix<int, 4> src;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      src(row, col) = static_cast<int>(100 * row + col + 1);
    }
  }
  vcp::square_matrix<int, 0> dst(src);
  REQUIRE(dst.size() == 4);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      REQUIRE(dst(row, col) == src(row, col));
    }
  }
  // The (n-1, n-1) cell specifically: this is what the old end-iterator
  // form silently dropped.
  REQUIRE(dst(3, 3) == src(3, 3));
  REQUIRE(dst(3, 3) != 0);
}

TEST_CASE("fixed-to-dynamic operator= copies every cell without overflow",
          "[square_matrix][conversion]") {
  // Regression for: data.resize(n) followed by std::copy of n*n elements
  // — a heap buffer overflow. ASan/UBSan builds would trap immediately;
  // non-sanitized builds would silently smash the heap. This test
  // asserts that the post-assignment size is n (dimension) and that
  // every source cell survives.
  vcp::square_matrix<int, 4> src;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      src(row, col) = static_cast<int>(10 * row + col + 1);
    }
  }
  vcp::square_matrix<int, 0> dst;
  dst = src;
  REQUIRE(dst.size() == 4);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      REQUIRE(dst(row, col) == src(row, col));
    }
  }
}

TEST_CASE("fixed-to-dynamic operator= re-sizes an existing dynamic matrix",
          "[square_matrix][conversion]") {
  // Assigning into a dynamic that previously held a different size.
  vcp::square_matrix<int, 0> dst(2);
  dst(0, 0) = 99;
  vcp::square_matrix<int, 3> src;
  src(2, 2) = 7;
  dst = src;
  REQUIRE(dst.size() == 3);
  REQUIRE(dst(2, 2) == 7);
}

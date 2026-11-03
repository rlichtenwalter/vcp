// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

// Standalone consumer compiled against an *installed* vcp tree only.
//
// The rest of the test suite consumes the headers through the library's
// BUILD_INTERFACE, i.e. straight out of the source tree, so it stays green
// even when the FILE_SET HEADERS list omits a header and the installation is
// therefore uncompilable. This translation unit is built by a separate CMake
// project that locates vcp via find_package(), with no path into the
// source tree at all. It includes every public header.

#include <iostream>
#include <sstream>

#include <vcp/detail/dense_or_sparse_map.hpp>
#include <vcp/detail/v3_buffer_bound.hpp>
#include <vcp/vcp.hpp>

int main() {
  // A triangle: the pair (0, 1) has exactly one third vertex, adjacent to both.
  std::istringstream input{"1 2\n0 2\n0 1\n"};
  vcp::graph g;
  input >> g;

  vcp::vcp<3, 1, false> calculator{g};
  auto const counts = calculator.generate_vector(g.vertices_begin(), g.vertices_begin() + 1);

  if (counts[7] != 1) {
    std::cerr << "unexpected triangle count\n";
    return 1;
  }

  std::cout << "installed-tree consumer counted the triangle\n";
  return 0;
}

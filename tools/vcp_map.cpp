// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

#include <CLI/CLI.hpp>

#include <vcp/vcp_static_mapper.hpp>

int main(int argc, char *argv[]) {
  try {
    CLI::App app{"Output the subgraph-to-element mapping for a particular VCP."};
    app.set_version_flag("--version",
                         std::string("vcp_map by Ryan N. Lichtenwalter v") + VCP_VERSION);

    std::size_t n = 0;
    std::size_t r = 0;
    std::size_t d = 0;
    std::uint64_t max_bytes = 0;

    app.add_option("n", n, "Number of vertices in the VCP")
        ->required()
        ->check(CLI::IsMember({3, 4, 5, 6, 7, 8}));
    app.add_option("r", r, "Number of relations in the VCP")
        ->required()
        ->check(CLI::IsMember({1, 2, 3, 4, 5, 6, 7, 8}));
    app.add_option("d", d, "Whether the VCP considers directedness (0 or 1)")
        ->required()
        ->check(CLI::IsMember({0, 1}));
    app.add_option(
           "-m,--mmax", max_bytes,
           "If more than this much memory would be required for computation, do not proceed "
           "with computation. Print the number of subgraphs in the VCP and exit. Accepts numeric "
           "suffixes (e.g. 1024, 500k, 2M, 4G).")
        ->transform(CLI::AsSizeValue(true));

    CLI11_PARSE(app, argc, argv);

    std::size_t const max_bytes_value = (max_bytes == 0) ? std::numeric_limits<std::size_t>::max()
                                                         : static_cast<std::size_t>(max_bytes);

    bool const directed = (d != 0);
    std::size_t const subgraph_count = vcp::vcp_static_mapper::subgraph_count(n, r, directed);
    if (subgraph_count * sizeof(std::size_t) > max_bytes_value) {
      std::cout << subgraph_count << '\n';
      return 0;
    }
    vcp::vcp_static_mapper mapper(n, r, directed);
    std::cout << mapper;

    return 0;
  } catch (std::exception const &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}

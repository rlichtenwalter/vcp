// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#include <iostream>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include <vcp/directed_graph.hpp>

int main(int argc, char *argv[]) {
  try {
    CLI::App app{
        "Convert a directed graph to an undirected graph by considering all unidirectional "
        "edges bidirectional."};
    app.set_version_flag("--version",
                         std::string("directed_to_undirected by Ryan N. Lichtenwalter v") +
                             VCP_VERSION);

    bool bidirectional = false;
    app.add_flag("-b,--bidirectional", bidirectional,
                 "Only bidirectional edges become undirected edges in the output. Unidirectional "
                 "edges are removed.");

    CLI11_PARSE(app, argc, argv);

    vcp::directed_graph g;
    std::cin >> g;
    std::vector<vcp::const_vertex_iterator> neighbors;
    for (vcp::const_vertex_iterator vIt(g.vertices_begin()); vIt != g.vertices_end(); ++vIt) {
      vcp::const_edge_iterator outIt = g.out_neighbors_begin(vIt);
      vcp::const_edge_iterator outEnd = g.out_neighbors_end(vIt);
      vcp::const_edge_iterator inIt = g.in_neighbors_begin(vIt);
      vcp::const_edge_iterator inEnd = g.in_neighbors_end(vIt);
      while (outIt != outEnd && inIt != inEnd) {
        if (g.target_of(outIt) < g.target_of(inIt)) {
          if (!bidirectional) {
            neighbors.push_back(g.target_of(outIt));
          }
          ++outIt;
        } else if (g.target_of(outIt) > g.target_of(inIt)) {
          if (!bidirectional) {
            neighbors.push_back(g.target_of(inIt));
          }
          ++inIt;
        } else {
          neighbors.push_back(g.target_of(outIt));
          ++outIt;
          ++inIt;
        }
      }
      if (!bidirectional) {
        while (outIt != outEnd) {
          neighbors.push_back(g.target_of(outIt));
          ++outIt;
        }
        while (inIt != inEnd) {
          neighbors.push_back(g.target_of(inIt));
          ++inIt;
        }
      }
      // Guard against isolated vertices: neighbors.end() - 1 is ill-formed
      // when the vector is empty (pointer arithmetic below begin()).
      if (!neighbors.empty()) {
        std::vector<vcp::const_vertex_iterator>::const_iterator neighbors_it(neighbors.begin());
        while (neighbors_it < neighbors.end() - 1) {
          std::cout << g.vertex_id(*neighbors_it++) << ' ';
        }
        std::cout << g.vertex_id(*neighbors_it);
      }
      std::cout << '\n';
      neighbors.clear();
    }

    return 0;
  } catch (std::exception const &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}

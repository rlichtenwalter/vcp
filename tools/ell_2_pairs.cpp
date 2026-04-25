// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#include <iostream>
#include <set>
#include <string>

#include <CLI/CLI.hpp>

#include <vcp/graph.hpp>

int main(int argc, char *argv[]) {
  try {
    CLI::App app{"Print, in lexicographical order, all pairs of nodes that are two hops distant "
                 "from each other."};
    app.set_version_flag("--version",
                         std::string("ell_2_pairs by Ryan N. Lichtenwalter v") + VCP_VERSION);

    CLI11_PARSE(app, argc, argv);

    vcp::graph g;
    std::cin >> g;
    for (vcp::const_vertex_iterator vIt(g.vertices_begin()); vIt != g.vertices_end(); ++vIt) {
      std::set<vcp::const_vertex_iterator> neighbors;
      for (vcp::const_edge_iterator e1It(g.neighbors_begin(vIt)); e1It != g.neighbors_end(vIt);
           ++e1It) {
        for (vcp::const_edge_iterator e2It(g.neighbors_begin(g.target_of(e1It)));
             e2It != g.neighbors_end(g.target_of(e1It)); ++e2It) {
          if (vIt < g.target_of(e2It) && !g.edge_exists(vIt, g.target_of(e2It))) {
            neighbors.insert(g.target_of(e2It));
          }
        }
      }
      for (auto neighbor : neighbors) {
        std::cout << g.vertex_id(vIt) << ' ' << g.vertex_id(neighbor) << '\n';
      }
    }

    return 0;
  } catch (std::exception const &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}

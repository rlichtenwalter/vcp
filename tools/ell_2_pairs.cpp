/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

You should have received a copy of the GNU General Public License along with
the Vertex Collocation Profiles code base. If not, see
<http://www.gnu.org/licenses/>.
*/

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

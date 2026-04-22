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

#include <array>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <CLI/CLI.hpp>

#include <vcp/directed_graph.hpp>
#include <vcp/graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/vcp.hpp>

template <std::size_t n>
std::ostream &operator<<(std::ostream &os, std::array<unsigned long, n> const &array) {
  if (array.empty()) {
    return os;
  }
  os << array[0];
  for (std::size_t i(1); i < n; ++i) {
    os << ' ' << array[i];
  }
  os << '\n';
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::map<T, unsigned long> const &map) {
  if (map.empty()) {
    return os;
  }
  auto it(map.begin());
  os << it->first << ',' << it->second;
  ++it;
  for (; it != map.end(); ++it) {
    os << ' ' << it->first << ',' << it->second;
  }
  os << '\n';
  return os;
}

int main(int argc, char *argv[]) {
  try {
    CLI::App app{"Output VCP vectors for pairs read from standard input."};
    app.set_version_flag("--version",
                         std::string("vcp_generate by Ryan N. Lichtenwalter v") + VCP_VERSION);

    std::size_t n = 0;
    std::size_t r = 0;
    std::size_t d = 0;
    std::string filename;

    app.add_option("n", n, "Number of vertices in the VCP")
        ->required()
        ->check(CLI::IsMember({3, 4, 5, 6, 7, 8}));
    app.add_option("r", r, "Number of relations in the VCP (1 or more)")
        ->required()
        ->check(CLI::PositiveNumber);
    app.add_option("d", d, "Whether the VCP considers directedness (0 or 1)")
        ->required()
        ->check(CLI::IsMember({0, 1}));
    app.add_option("graph_filename", filename, "Path to the file containing the graph")
        ->required()
        ->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);

    std::ifstream file(filename, std::ifstream::in);
    if (!file) {
      std::cerr << "error opening file: " << filename << '\n';
      return 1;
    }

    vcp::vertex_id_t v1;
    vcp::vertex_id_t v2;
    bool const directed = (d != 0);

    if (directed) {
      if (n == 3) {
        if (r == 1) {
          vcp::directed_graph g;
          file >> g;
          vcp::vcp<3, 1, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 2) {
          vcp::multirelational_directed_graph<2> g;
          file >> g;
          vcp::vcp<3, 2, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 30) {
          vcp::multirelational_directed_graph<30> g;
          file >> g;
          vcp::vcp<3, 30, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        }
      } else if (n == 4) {
        if (r == 1) {
          vcp::directed_graph g;
          file >> g;
          vcp::vcp<4, 1, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 2) {
          vcp::multirelational_directed_graph<2> g;
          file >> g;
          vcp::vcp<4, 2, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 30) {
          vcp::multirelational_directed_graph<30> g;
          file >> g;
          vcp::vcp<4, 30, true> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        }
      }
    } else {
      if (n == 3) {
        if (r == 1) {
          vcp::graph g;
          file >> g;
          vcp::vcp<3, 1, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 2) {
          vcp::multirelational_graph<2> g;
          file >> g;
          vcp::vcp<3, 2, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 30) {
          vcp::multirelational_graph<30> g;
          file >> g;
          vcp::vcp<3, 30, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        }
      } else if (n == 4) {
        if (r == 1) {
          vcp::graph g;
          file >> g;
          vcp::vcp<4, 1, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 2) {
          vcp::multirelational_graph<2> g;
          file >> g;
          vcp::vcp<4, 2, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        } else if (r == 30) {
          vcp::multirelational_graph<30> g;
          file >> g;
          vcp::vcp<4, 30, false> profiler(g);
          while (std::cin >> v1 >> v2) {
            std::cout << profiler.generate_vector(
                vcp::const_vertex_iterator(g.vertices_begin() + v1),
                vcp::const_vertex_iterator(g.vertices_begin() + v2));
          }
        }
      }
    }

    return 0;
  } catch (std::exception const &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}

/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation
Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_GRAPH_H
#define VCP_GRAPH_H

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace vcp {

/** @brief Integer type used to identify a vertex by its position in the CSR arrays. */
using vertex_id_t = std::size_t;

/** @brief Integer type used to identify an edge by its position in the CSR arrays. */
using edge_id_t = std::size_t;

/**
 * @brief Opaque iterator type over vertices in the CSR layout.
 *
 * Internally a pointer-to-pointer; advance with `++` to move to the next
 * vertex. The difference of two vertex iterators equals the vertex-id gap.
 */
using const_vertex_iterator = void *const *;

/**
 * @brief Opaque iterator type over edges (or neighbor lists) in the CSR layout.
 *
 * Internally a pointer-to-pointer; advance with `++` to move to the next
 * edge slot. The range [neighbors_begin(v), neighbors_end(v)) gives all
 * neighbors of vertex @p v.
 */
using const_edge_iterator = void *const *;

/**
 * @brief Undirected graph in compressed sparse row (CSR) format.
 *
 * Each undirected edge {u, v} is stored twice: once in u's adjacency list
 * and once in v's. The `vertices` array holds `vertex_count() + 1`
 * pointers into `edges`; `edges` holds `2 * edge_count()` pointers back
 * into `vertices`. This pointer-based CSR avoids index arithmetic and
 * enables O(1) neighbor-range access without index recomputation.
 *
 * Graphs are populated by reading from a stream via `operator>>`.
 * The text format is one line per vertex containing space-separated
 * neighbor ids (0-based), with an empty line representing an isolated vertex.
 *
 * Neighbor lists are expected to be sorted; the `vcp<3, 1, false>` and
 * `vcp<4, 1, false>` specializations rely on sorted adjacency for their
 * merge-based enumeration.
 */
class graph {
public:
  /** @brief Construct an empty graph with no vertices or edges. */
  graph();

  /** @brief Copy-construct a graph, deep-copying both CSR arrays. */
  graph(graph const &);

  /** @brief Copy-assign a graph, deep-copying both CSR arrays. */
  graph &operator=(graph const &);

  ~graph();

  /** @brief Return the number of vertices. */
  std::size_t vertex_count() const;

  /** @brief Return the number of undirected edges (each pair counted once). */
  std::size_t edge_count() const;

  /** @brief Return an iterator to the first vertex. */
  const_vertex_iterator vertices_begin() const;

  /** @brief Return a past-the-end iterator for vertices. */
  const_vertex_iterator vertices_end() const;

  /** @brief Return an iterator to the first edge slot in the CSR edge array. */
  const_edge_iterator edges_begin() const;

  /** @brief Return a past-the-end iterator for the CSR edge array. */
  const_edge_iterator edges_end() const;

  /**
   * @brief Return an iterator to the first neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first neighbor edge slot.
   */
  const_edge_iterator neighbors_begin(const_vertex_iterator it) const;

  /**
   * @brief Return a past-the-end iterator for the neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the neighbor range.
   */
  const_edge_iterator neighbors_end(const_vertex_iterator it) const;

  /**
   * @brief Return the integer id of the given vertex.
   *
   * @param it Vertex iterator obtained from this graph.
   * @return Zero-based vertex identifier equal to the iterator's offset from vertices_begin().
   */
  vertex_id_t vertex_id(const_vertex_iterator it) const;

  /**
   * @brief Return the vertex pointed to by an edge iterator.
   *
   * @param it Edge iterator within a neighbor range.
   * @return Iterator to the target vertex of the edge.
   */
  const_vertex_iterator target_of(const_edge_iterator it) const;

  /**
   * @brief Return the integer id of the given edge slot.
   *
   * @param it Edge iterator obtained from this graph.
   * @return Zero-based offset of @p it within the CSR edge array.
   */
  edge_id_t edge_id(const_edge_iterator it) const;

  /**
   * @brief Return true if @p it refers to an existing edge (i.e., is not edges_end()).
   *
   * @param it Edge iterator to test.
   * @return True if @p it != edges_end().
   */
  bool edge_exists(const_edge_iterator it) const;

  /**
   * @brief Find the edge between two vertices and return an iterator to it.
   *
   * Performs a linear search through @p source's neighbor list. Returns
   * edges_end() if no edge exists.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return Iterator to the edge, or edges_end() if absent.
   */
  const_edge_iterator edge(const_vertex_iterator source, const_vertex_iterator target) const;

  /**
   * @brief Return true if an edge exists between @p source and @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the edge is present.
   */
  bool edge_exists(const_vertex_iterator source, const_vertex_iterator target) const;

  friend std::ostream &operator<<(std::ostream &, graph const &);
  friend std::istream &operator>>(std::istream &, graph &);

private:
  std::size_t num_vertices{0};
  std::size_t num_edges{0};
  std::unique_ptr<void *[]> vertices;
  std::unique_ptr<void *[]> edges;
};

inline graph::graph()
    : vertices(std::unique_ptr<void *[]>(new void *[1])),
      edges(std::unique_ptr<void *[]>(new void *[1])) {
  vertices[0] = static_cast<void *>(&edges[0]);
  edges[0] = nullptr;
}

inline graph::graph(graph const &g)
    : num_vertices(g.num_vertices), num_edges(g.num_edges),
      vertices(std::unique_ptr<void *[]>(new void *[g.vertex_count() + 1])),
      edges(std::unique_ptr<void *[]>(new void *[g.num_edges + 1])) {
  for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
    vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.neighbors_begin(it))]);
  }
  vertices[vertex_count()] = static_cast<void *>(&edges[num_edges]);
  for (const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
  }
  edges[num_edges] = nullptr;
}

inline graph::~graph() = default;

inline graph &graph::operator=(graph const &g) {
  if (this != &g) {
    num_vertices = g.num_vertices;
    num_edges = g.num_edges;
    vertices = std::unique_ptr<void *[]>(new void *[g.vertex_count() + 1]);
    edges = std::unique_ptr<void *[]>(new void *[g.num_edges + 1]);
    for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
      vertex_id_t id = g.vertex_id(it);
      vertices[id] = static_cast<void *>(&edges[g.edge_id(g.neighbors_begin(it))]);
    }
    vertices[vertex_count()] = static_cast<void *>(&edges[num_edges]);
    for (const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it) {
      edge_id_t id = g.edge_id(it);
      edges[id] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    }
    edges[num_edges] = nullptr;
  }
  return *this;
}

inline std::size_t graph::vertex_count() const { return num_vertices; }

inline std::size_t graph::edge_count() const { return num_edges / 2; }

inline const_vertex_iterator graph::vertices_begin() const { return &vertices[0]; }

inline const_vertex_iterator graph::vertices_end() const { return &vertices[vertex_count()]; }

inline const_vertex_iterator graph::edges_begin() const { return &edges[0]; }

inline const_vertex_iterator graph::edges_end() const { return &edges[num_edges]; }

inline const_edge_iterator graph::neighbors_begin(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*it);
}

inline const_edge_iterator graph::neighbors_end(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*(it + 1));
}

inline vertex_id_t graph::vertex_id(const_vertex_iterator it) const {
  return it - static_cast<const_vertex_iterator>(vertices.get());
}

inline const_vertex_iterator graph::target_of(const_edge_iterator it) const {
  return static_cast<const_vertex_iterator>(*it);
}

inline edge_id_t graph::edge_id(const_edge_iterator it) const {
  return it - static_cast<const_edge_iterator>(edges.get());
}

inline bool graph::edge_exists(const_edge_iterator it) const { return it != edges_end(); }

inline const_edge_iterator graph::edge(const_vertex_iterator source,
                                       const_vertex_iterator target) const {
  const_edge_iterator it(std::find(neighbors_begin(source), neighbors_end(source), target));
  return it == neighbors_end(source) ? edges_end() : it;
}

inline bool graph::edge_exists(const_vertex_iterator source, const_vertex_iterator target) const {
  return neighbors_end(source) != std::find(neighbors_begin(source), neighbors_end(source), target);
}

/**
 * @brief Write a graph to an output stream.
 *
 * Emits one line per vertex containing space-separated neighbor ids. An
 * isolated vertex produces an empty line. The format is the same as that
 * accepted by `operator>>`.
 *
 * @param os Output stream.
 * @param g  Graph to write.
 * @return Reference to @p os.
 */
inline std::ostream &operator<<(std::ostream &os, graph const &g) {
  for (const_vertex_iterator vIt = g.vertices_begin(); vIt != g.vertices_end(); ++vIt) {
    const_edge_iterator const nBegin = g.neighbors_begin(vIt);
    const_edge_iterator const nEnd = g.neighbors_end(vIt);
    for (const_edge_iterator nIt = nBegin; nIt != nEnd; ++nIt) {
      if (nIt != nBegin) {
        os << ' ';
      }
      os << g.vertex_id(g.target_of(nIt));
    }
    os << '\n';
  }
  return os;
}

/**
 * @brief Read a graph from an input stream.
 *
 * Expects one line per vertex with space-separated neighbor ids (0-based).
 * Each undirected edge must appear in both directions (i.e., if vertex u
 * lists v as a neighbor, v must also list u). An empty line represents an
 * isolated vertex.
 *
 * @param is Input stream positioned at the first vertex line.
 * @param g  Graph object to populate; any previous contents are replaced.
 * @return Reference to @p is.
 * @throws std::invalid_argument If a token cannot be parsed as a vertex id.
 */
inline std::istream &operator>>(std::istream &is, graph &g) {
  std::vector<edge_id_t> v_temp;
  std::vector<vertex_id_t> e_temp;
  std::string s1;
  std::size_t edge_id(0);
  while (std::getline(is, s1, '\n')) {
    v_temp.push_back(edge_id);
    std::istringstream iss(s1);
    std::string s2;
    while (std::getline(iss, s2, ' ')) {
      // atol previously used here silently returns 0 on malformed input and
      // overflows are UB. from_chars gives explicit errors and no locale
      // dependency. Empty tokens (from trailing or doubled spaces) are
      // preserved as zero to match the historical atol-on-empty behavior.
      vertex_id_t parsed(0);
      if (!s2.empty()) {
        auto const [ptr, ec] = std::from_chars(s2.data(), s2.data() + s2.size(), parsed);
        if (ec != std::errc{} || ptr != s2.data() + s2.size()) {
          throw std::invalid_argument("graph parse: invalid vertex id '" + s2 + "'");
        }
      }
      e_temp.push_back(parsed);
      ++edge_id;
    }
  }

  g.num_vertices = v_temp.size();
  g.num_edges = e_temp.size();

  g.vertices = std::unique_ptr<void *[]>(new void *[g.vertex_count() + 1]);
  g.edges = std::unique_ptr<void *[]>(new void *[g.num_edges + 1]);

  for (size_t i = 0; i < g.vertex_count(); ++i) {
    g.vertices[i] = static_cast<void *>(&g.edges[v_temp.at(i)]);
  }
  g.vertices[g.vertex_count()] = static_cast<void *>(&g.edges[g.num_edges]);
  for (size_t i = 0; i < g.num_edges; ++i) {
    g.edges[i] = static_cast<void *>(&g.vertices[e_temp.at(i)]);
  }
  g.edges[g.num_edges] = nullptr;

  return is;
}

} // namespace vcp

#endif

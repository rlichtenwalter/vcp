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

#ifndef VCP_MULTIRELATIONAL_GRAPH
#define VCP_MULTIRELATIONAL_GRAPH

#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vcp/graph.hpp>
#include <vector>

namespace vcp {

using vertex_id_t = std::size_t;
using edge_id_t = std::size_t;
using const_vertex_iterator = void *const *;
using const_edge_iterator = void *const *;

/**
 * @brief Undirected multirelational graph in CSR format with per-edge relation bitmasks.
 *
 * Extends `graph` to support `r` independent edge relations. Each edge stores
 * an `r`-bit bitmask (`connectivity_address_type`) where bit `k` is set if
 * relation `k` is active on that edge. The CSR structure is identical to
 * `graph`; an auxiliary `edge_values` array holds the bitmask for each edge slot.
 *
 * For `r <= CHAR_BIT * sizeof(std::size_t)` (typically 64), the bitmask fits
 * in a `std::size_t` and no external library is needed. For larger `r` a
 * fixed-width Boost.Multiprecision unsigned integer is used instead.
 *
 * The text format extends `graph`'s format: each neighbor is written as
 * `<target_id>,<bitmask>`, with neighbors separated by spaces.
 *
 * @tparam r Maximum number of edge relations (capacity); the actual number
 *           of active relations in a loaded graph may be less.
 */
template <std::size_t r> class multirelational_graph {
public:
  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   *
   * Resolves to `std::size_t` when r fits within a machine word, and to a
   * fixed-width Boost.Multiprecision unsigned integer otherwise.
   */
  using connectivity_address_type = typename std::conditional<
      r <= CHAR_BIT * sizeof(std::size_t), std::size_t,
      boost::multiprecision::number<
          boost::multiprecision::cpp_int_backend<r, r, boost::multiprecision::unsigned_magnitude,
                                                 boost::multiprecision::unchecked, void>>>::type;

  /** @brief Construct an empty graph with no vertices or edges. */
  multirelational_graph();

  /** @brief Copy-construct a graph, deep-copying CSR arrays and edge-value array. */
  multirelational_graph(multirelational_graph const &);

  ~multirelational_graph();

  /** @brief Copy-assign a graph, deep-copying CSR arrays and edge-value array. */
  multirelational_graph &operator=(multirelational_graph const &);

  /** @brief Return the number of vertices. */
  std::size_t vertex_count() const;

  /** @brief Return the number of undirected edges (each pair counted once). */
  std::size_t edge_count() const;

  /**
   * @brief Return the number of active relations in the loaded graph.
   *
   * Counts the minimum number of relation bits needed to represent the
   * highest-valued edge bitmask. This is a property of the data, not of
   * the template parameter `r` (which is the capacity). Returns 0 for an
   * empty graph or one whose edges all have bitmask 0.
   */
  std::size_t relation_count() const;

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
   * @return Zero-based vertex identifier.
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
   * @brief Return the relation bitmask stored on the edge at @p it.
   *
   * Returns 0 if @p it equals edges_end() (the "no edge" sentinel), matching
   * the convention used by VCP algorithms that call this without a null check.
   *
   * @param it Edge iterator within the CSR edge array.
   * @return The r-bit relation bitmask, or 0 if @p it == edges_end().
   */
  connectivity_address_type edge_value(const_edge_iterator it) const;

  /**
   * @brief Return true if an edge exists between @p source and @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the edge is present.
   */
  bool edge_exists(const_vertex_iterator source, const_vertex_iterator target) const;

  template <std::size_t r_>
  friend std::ostream &operator<<(std::ostream &, multirelational_graph<r_> const &);
  template <std::size_t r_>
  friend std::istream &operator>>(std::istream &, multirelational_graph<r_> &);

private:
  std::size_t num_vertices{0};
  std::size_t num_edges{0};
  std::unique_ptr<void *[]> vertices;
  std::unique_ptr<void *[]> edges;
  std::unique_ptr<connectivity_address_type[]> edge_values;
};

template <std::size_t r>
multirelational_graph<r>::multirelational_graph()
    : vertices(std::make_unique<void *[]>(1)), edges(std::make_unique<void *[]>(1)),
      edge_values(
          std::make_unique<typename multirelational_graph<r>::connectivity_address_type[]>(1)) {
  vertices[0] = static_cast<void *>(&edges[0]);
  edges[0] = nullptr;
}

template <std::size_t r>
multirelational_graph<r>::multirelational_graph(multirelational_graph const &g)
    : num_vertices(g.num_vertices), num_edges(g.num_edges),
      vertices(std::make_unique<void *[]>(g.vertex_count() + 1)),
      edges(std::make_unique<void *[]>(g.num_edges + 1)),
      edge_values(std::make_unique<typename multirelational_graph<r>::connectivity_address_type[]>(
          g.num_edges)) {
  for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
    vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.neighbors_begin(it))]);
  }
  vertices[vertex_count()] = static_cast<void *>(&edges[num_edges]);
  for (const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    edge_values[g.edge_id(it)] = g.edge_value(it);
  }
  edges[num_edges] = nullptr;
}

template <std::size_t r> multirelational_graph<r>::~multirelational_graph() = default;

template <std::size_t r>
multirelational_graph<r> &multirelational_graph<r>::operator=(multirelational_graph<r> const &g) {
  if (this != &g) {
    num_vertices = g.num_vertices;
    num_edges = g.num_edges;
    vertices = std::make_unique<void *[]>(g.vertex_count() + 1);
    edges = std::make_unique<void *[]>(g.num_edges + 1);
    edge_values = std::make_unique<typename multirelational_graph<r>::connectivity_address_type[]>(
        g.num_edges);
    for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
      vertex_id_t id = g.vertex_id(it);
      vertices[id] = static_cast<void *>(&edges[g.edge_id(g.neighbors_begin(it))]);
    }
    vertices[vertex_count()] = static_cast<void *>(&edges[num_edges]);
    for (const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it) {
      edge_id_t id = g.edge_id(it);
      edges[id] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
      edge_values[g.edge_id(it)] = g.edge_value(it);
    }
    edges[num_edges] = nullptr;
  }
  return *this;
}

template <std::size_t r> std::size_t multirelational_graph<r>::vertex_count() const {
  return num_vertices;
}

template <std::size_t r> std::size_t multirelational_graph<r>::edge_count() const {
  return num_edges / 2;
}

template <std::size_t r> std::size_t multirelational_graph<r>::relation_count() const {
  // Number of bits required to represent the largest stored edge value —
  // i.e. the index of the highest-numbered relation with at least one edge.
  // This is a property of the data, not the template parameter r (which is
  // the maximum capacity). Returns 0 for an empty graph or one whose edges
  // all have value 0 (no relations set). Bit-counting rather than
  // std::log2 is used because connectivity_address_type may be a
  // boost::multiprecision integer, for which std::log2 is not defined.
  if (num_edges == 0) {
    return 0;
  }
  connectivity_address_type const max_val(
      *std::max_element(&edge_values[0], &edge_values[num_edges]));
  std::size_t bits(0);
  for (connectivity_address_type v(max_val); v > 0; v >>= 1) {
    ++bits;
  }
  return bits;
}

template <std::size_t r> const_vertex_iterator multirelational_graph<r>::vertices_begin() const {
  return &vertices[0];
}

template <std::size_t r> const_vertex_iterator multirelational_graph<r>::vertices_end() const {
  return &vertices[vertex_count()];
}

template <std::size_t r> const_vertex_iterator multirelational_graph<r>::edges_begin() const {
  return &edges[0];
}

template <std::size_t r> const_vertex_iterator multirelational_graph<r>::edges_end() const {
  return &edges[num_edges];
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::neighbors_begin(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*it);
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::neighbors_end(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*(it + 1));
}

template <std::size_t r>
vertex_id_t multirelational_graph<r>::vertex_id(const_vertex_iterator it) const {
  return it - static_cast<const_vertex_iterator>(vertices.get());
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::target_of(const_edge_iterator it) const {
  return static_cast<const_vertex_iterator>(*it);
}

template <std::size_t r> edge_id_t multirelational_graph<r>::edge_id(const_edge_iterator it) const {
  return it - static_cast<const_edge_iterator>(edges.get());
}

template <std::size_t r> bool multirelational_graph<r>::edge_exists(const_edge_iterator it) const {
  return it != edges_end();
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::edge(const_vertex_iterator source,
                                                   const_vertex_iterator target) const {
  const_edge_iterator it(std::find(neighbors_begin(source), neighbors_end(source), target));
  return it == neighbors_end(source) ? edges_end() : it;
}

template <std::size_t r>
typename multirelational_graph<r>::connectivity_address_type
multirelational_graph<r>::edge_value(const_edge_iterator it) const {
  // A non-existent edge has connectivity value 0 (unconnected), matching the
  // contract of graph::edge_value / directed_graph::edge_value. Callers such
  // as vcp::edge_value() pipe edge() directly here without a null check.
  return it == edges_end() ? 0 : edge_values[edge_id(it)];
}

template <std::size_t r>
bool multirelational_graph<r>::edge_exists(const_vertex_iterator source,
                                           const_vertex_iterator target) const {
  return neighbors_end(source) != std::find(neighbors_begin(source), neighbors_end(source), target);
}

/**
 * @brief Write a multirelational graph to an output stream.
 *
 * Emits one line per vertex. Each neighbor is written as `<id>,<bitmask>`
 * with neighbors separated by spaces. An isolated vertex produces an empty line.
 * The format matches the input accepted by `operator>>`.
 *
 * @tparam r Number of edge relation bits.
 * @param os Output stream.
 * @param g  Graph to write.
 * @return Reference to @p os.
 */
template <std::size_t r>
std::ostream &operator<<(std::ostream &os, multirelational_graph<r> const &g) {
  for (const_vertex_iterator vIt = g.vertices_begin(); vIt != g.vertices_end(); ++vIt) {
    const_edge_iterator const nBegin = g.neighbors_begin(vIt);
    const_edge_iterator const nEnd = g.neighbors_end(vIt);
    for (const_edge_iterator nIt = nBegin; nIt != nEnd; ++nIt) {
      if (nIt != nBegin) {
        os << ' ';
      }
      os << g.vertex_id(g.target_of(nIt)) << ',' << g.edge_value(nIt);
    }
    os << '\n';
  }
  return os;
}

/**
 * @brief Read a multirelational graph from an input stream.
 *
 * Expects one line per vertex with space-separated entries of the form
 * `<neighbor_id>,<bitmask>`. Each undirected edge must appear in both
 * directions. An empty line represents an isolated vertex.
 *
 * @tparam r Number of edge relation bits.
 * @param is Input stream positioned at the first vertex line.
 * @param g  Graph object to populate; any previous contents are replaced.
 * @return Reference to @p is.
 */
template <std::size_t r> std::istream &operator>>(std::istream &is, multirelational_graph<r> &g) {
  std::vector<edge_id_t> v_temp;
  std::vector<std::pair<vertex_id_t, typename multirelational_graph<r>::connectivity_address_type>>
      e_temp;
  std::string s1;
  std::size_t edge_id(0);
  while (std::getline(is, s1, '\n')) {
    v_temp.push_back(edge_id);
    std::istringstream iss1(s1);
    std::string s2;
    while (std::getline(iss1, s2, ' ')) {
      std::istringstream iss2(s2);
      std::string neighbor_str;
      std::string value_str;
      getline(iss2, neighbor_str, ',');
      getline(iss2, value_str, ',');
      std::istringstream iss3(neighbor_str);
      std::istringstream iss4(value_str);
      vertex_id_t neighbor;
      typename multirelational_graph<r>::connectivity_address_type value;
      iss3 >> neighbor;
      iss4 >> value;
      e_temp.push_back(std::make_pair(neighbor, value));
      ++edge_id;
    }
  }

  g.num_vertices = v_temp.size();
  g.num_edges = e_temp.size();

  g.vertices = std::make_unique<void *[]>(g.vertex_count() + 1);
  g.edges = std::make_unique<void *[]>(g.num_edges + 1);
  g.edge_values =
      std::make_unique<typename multirelational_graph<r>::connectivity_address_type[]>(g.num_edges);

  for (size_t i(0); i < g.vertex_count(); ++i) {
    g.vertices[i] = static_cast<void *>(&g.edges[v_temp[i]]);
  }
  g.vertices[g.vertex_count()] = static_cast<void *>(&g.edges[g.num_edges]);
  for (size_t i(0); i < g.num_edges; ++i) {
    g.edges[i] = static_cast<void *>(&g.vertices[e_temp[i].first]);
    g.edge_values[i] = e_temp[i].second;
  }
  g.edges[g.num_edges] = nullptr;

  return is;
}

} // namespace vcp

#endif

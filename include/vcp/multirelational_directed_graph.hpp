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

#ifndef VCP_MULTIRELATIONAL_DIRECTED_GRAPH_HPP
#define VCP_MULTIRELATIONAL_DIRECTED_GRAPH_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vcp/multirelational_graph.hpp>
#include <vector>

namespace vcp {

using vertex_id_t = std::size_t;
using edge_id_t = std::size_t;
using const_vertex_iterator = void *const *;
using const_edge_iterator = void *const *;

/**
 * @brief Directed multirelational graph with per-edge relation bitmasks and bidirectional CSR.
 *
 * Combines the bidirectional CSR layout of `directed_graph` with the per-edge
 * relation bitmask of `multirelational_graph`. Out-edges and in-edges are
 * stored in a single flat `edges` array. A parallel `edge_values` array stores
 * the r-bit connectivity bitmask for each edge slot (both out and in blocks).
 *
 * In-edges and their bitmasks are derived by inverting the out-edge list during
 * `operator>>`. The `edge_value()` method returns 0 for the sentinel in_edges_end()
 * iterator, matching the convention used by VCP algorithms.
 *
 * Note: `out_edge_count()` and `in_edge_count()` both return the number of directed
 * out-edges, since every out-edge implies exactly one in-edge.
 *
 * @tparam r Maximum number of edge relations (capacity).
 */
template <std::size_t r> class multirelational_directed_graph {
public:
  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   *
   * Alias of `multirelational_graph<r>::connectivity_address_type`.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /** @brief Construct an empty directed graph with no vertices or edges. */
  multirelational_directed_graph();

  /** @brief Copy-construct a graph, deep-copying CSR arrays and edge-value array. */
  multirelational_directed_graph(multirelational_directed_graph const &);

  ~multirelational_directed_graph();

  /** @brief Copy-assign a graph, deep-copying CSR arrays and edge-value array. */
  multirelational_directed_graph &operator=(multirelational_directed_graph const &);

  /** @brief Return the number of vertices. */
  std::size_t vertex_count() const;

  /** @brief Return the number of directed out-edges. */
  std::size_t out_edge_count() const;

  /**
   * @brief Return the number of directed in-edges.
   *
   * Always equals out_edge_count() since every out-edge implies one in-edge.
   */
  std::size_t in_edge_count() const;

  /**
   * @brief Return the number of active relations in the loaded graph.
   *
   * Counts the minimum number of relation bits needed to represent the
   * highest-valued edge bitmask across all out- and in-edge slots.
   * Returns 0 for an empty graph or one whose edges all have bitmask 0.
   */
  std::size_t relation_count() const;

  /** @brief Return an iterator to the first vertex. */
  const_vertex_iterator vertices_begin() const;

  /** @brief Return a past-the-end iterator for vertices. */
  const_vertex_iterator vertices_end() const;

  /** @brief Return an iterator to the first slot in the out-edge block. */
  const_edge_iterator out_edges_begin() const;

  /** @brief Return a past-the-end iterator for the out-edge block. */
  const_edge_iterator out_edges_end() const;

  /** @brief Return an iterator to the first slot in the in-edge block. */
  const_edge_iterator in_edges_begin() const;

  /** @brief Return a past-the-end iterator for the in-edge block. */
  const_edge_iterator in_edges_end() const;

  /**
   * @brief Return an iterator to the first out-neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first out-neighbor edge slot.
   */
  const_edge_iterator out_neighbors_begin(const_vertex_iterator it) const;

  /**
   * @brief Return a past-the-end iterator for out-neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the out-neighbor range.
   */
  const_edge_iterator out_neighbors_end(const_vertex_iterator it) const;

  /**
   * @brief Return an iterator to the first in-neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first in-neighbor edge slot.
   */
  const_edge_iterator in_neighbors_begin(const_vertex_iterator it) const;

  /**
   * @brief Return a past-the-end iterator for in-neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the in-neighbor range.
   */
  const_edge_iterator in_neighbors_end(const_vertex_iterator it) const;

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
   * @param it Edge iterator within an out- or in-neighbor range.
   * @return Iterator to the target vertex of the edge.
   */
  const_vertex_iterator target_of(const_edge_iterator it) const;

  /**
   * @brief Return the integer id of the given edge slot.
   *
   * @param it Edge iterator obtained from this graph.
   * @return Zero-based offset of @p it within the combined edge array.
   */
  edge_id_t edge_id(const_edge_iterator it) const;

  /**
   * @brief Return true if @p it refers to an existing edge (i.e., is not in_edges_end()).
   *
   * @param it Edge iterator to test.
   * @return True if @p it != in_edges_end().
   */
  bool edge_exists(const_edge_iterator it) const;

  /**
   * @brief Return the relation bitmask stored on the edge at @p it.
   *
   * Returns 0 if @p it equals in_edges_end() (the "no edge" sentinel),
   * matching the convention used by VCP algorithms.
   *
   * @param it Edge iterator within the combined edge array.
   * @return The r-bit relation bitmask, or 0 if @p it == in_edges_end().
   */
  connectivity_address_type edge_value(const_edge_iterator it) const;

  /**
   * @brief Find the directed out-edge from @p source to @p target.
   *
   * Performs a linear search through @p source's out-neighbor list.
   * Returns in_edges_end() if no such edge exists.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return Iterator to the out-edge, or in_edges_end() if absent.
   */
  const_edge_iterator out_edge(const_vertex_iterator source, const_vertex_iterator target) const;

  /**
   * @brief Find the directed in-edge from @p source to @p target.
   *
   * Performs a linear search through @p source's in-neighbor list.
   * Returns in_edges_end() if no such edge exists.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return Iterator to the in-edge, or in_edges_end() if absent.
   */
  const_edge_iterator in_edge(const_vertex_iterator source, const_vertex_iterator target) const;

  /**
   * @brief Return true if a directed out-edge exists from @p source to @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the out-edge is present.
   */
  bool out_edge_exists(const_vertex_iterator source, const_vertex_iterator target) const;

  /**
   * @brief Return true if a directed in-edge exists from @p source to @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the in-edge is present.
   */
  bool in_edge_exists(const_vertex_iterator source, const_vertex_iterator target) const;

  template <std::size_t r_>
  friend std::ostream &operator<<(std::ostream &, multirelational_directed_graph<r_> const &);
  template <std::size_t r_>
  friend std::istream &operator>>(std::istream &, multirelational_directed_graph<r_> &);

private:
  std::size_t num_vertices{0};
  std::size_t num_out_edges{0};
  std::unique_ptr<void *[]> vertices;
  std::unique_ptr<void *[]> edges;
  std::unique_ptr<connectivity_address_type[]> edge_values;
};

template <std::size_t r>
multirelational_directed_graph<r>::multirelational_directed_graph()
    : vertices(std::make_unique<void *[]>(1)), edges(std::make_unique<void *[]>(1)),
      edge_values(
          std::make_unique<typename multirelational_directed_graph<r>::connectivity_address_type[]>(
              1)) {
  vertices[0] = static_cast<void *>(&edges[0]);
  edges[0] = nullptr;
}

template <std::size_t r>
multirelational_directed_graph<r>::multirelational_directed_graph(
    multirelational_directed_graph const &g)
    : num_vertices(g.num_vertices), num_out_edges(g.num_out_edges),
      vertices(std::make_unique<void *[]>(2 * g.vertex_count() + 1)),
      edges(std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1)),
      edge_values(
          std::make_unique<typename multirelational_directed_graph<r>::connectivity_address_type[]>(
              2 * g.num_out_edges)) {
  for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
    vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.out_neighbors_begin(it))]);
    vertices[vertex_count() + g.vertex_id(it)] =
        static_cast<void *>(&edges[g.edge_id(g.in_neighbors_begin(it))]);
  }
  vertices[2 * vertex_count()] = static_cast<void *>(&edges[out_edge_count() + in_edge_count()]);
  for (const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    edge_values[g.edge_id(it)] = g.edge_value(it);
  }
  for (const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    edge_values[g.edge_id(it)] = g.edge_value(it);
  }
  edges[out_edge_count() + in_edge_count()] = nullptr;
}

template <std::size_t r>
multirelational_directed_graph<r>::~multirelational_directed_graph() = default;

template <std::size_t r>
multirelational_directed_graph<r> &
multirelational_directed_graph<r>::operator=(multirelational_directed_graph const &g) {
  if (this != &g) {
    num_vertices = g.num_vertices;
    num_out_edges = g.num_out_edges;
    vertices = std::make_unique<void *[]>(2 * g.vertex_count() + 1);
    edges = std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1);
    edge_values =
        std::make_unique<typename multirelational_directed_graph<r>::connectivity_address_type[]>(
            2 * g.num_out_edges);
    for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
      vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.out_neighbors_begin(it))]);
      vertices[vertex_count() + g.vertex_id(it)] =
          static_cast<void *>(&edges[g.edge_id(g.in_neighbors_begin(it))]);
    }
    vertices[2 * vertex_count()] = static_cast<void *>(&edges[out_edge_count() + in_edge_count()]);
    for (const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it) {
      edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
      edge_values[g.edge_id(it)] = g.edge_value(it);
    }
    for (const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it) {
      edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
      edge_values[g.edge_id(it)] = g.edge_value(it);
    }
    edges[out_edge_count() + in_edge_count()] = nullptr;
  }
  return *this;
}

template <std::size_t r> std::size_t multirelational_directed_graph<r>::vertex_count() const {
  return num_vertices;
}

template <std::size_t r> std::size_t multirelational_directed_graph<r>::out_edge_count() const {
  return num_out_edges;
}

template <std::size_t r> std::size_t multirelational_directed_graph<r>::in_edge_count() const {
  return num_out_edges;
}

template <std::size_t r> std::size_t multirelational_directed_graph<r>::relation_count() const {
  // Number of bits required to represent the largest stored edge value —
  // i.e. the index of the highest-numbered relation with at least one edge.
  // This is a property of the data, not the template parameter r (which is
  // the maximum capacity). Returns 0 for an empty graph or one whose edges
  // all have value 0 (no relations set). Bit-counting rather than
  // std::log2 is used because connectivity_address_type may be a
  // boost::multiprecision integer, for which std::log2 is not defined.
  if (num_out_edges == 0) {
    return 0;
  }
  connectivity_address_type const max_val(
      *std::max_element(&edge_values[0], &edge_values[2 * num_out_edges]));
  std::size_t bits(0);
  for (connectivity_address_type v(max_val); v > 0; v >>= 1) {
    ++bits;
  }
  return bits;
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::vertices_begin() const {
  return &vertices[0];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::vertices_end() const {
  return &vertices[vertex_count()];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::out_edges_begin() const {
  return &edges[0];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::out_edges_end() const {
  return &edges[out_edge_count()];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::in_edges_begin() const {
  return &edges[out_edge_count()];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::in_edges_end() const {
  return &edges[out_edge_count() + in_edge_count()];
}

template <std::size_t r>
const_edge_iterator
multirelational_directed_graph<r>::out_neighbors_begin(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*it);
}

template <std::size_t r>
const_edge_iterator
multirelational_directed_graph<r>::out_neighbors_end(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*(it + 1));
}

template <std::size_t r>
const_edge_iterator
multirelational_directed_graph<r>::in_neighbors_begin(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*(vertex_count() + it));
}

template <std::size_t r>
const_edge_iterator
multirelational_directed_graph<r>::in_neighbors_end(const_vertex_iterator it) const {
  return static_cast<const_edge_iterator>(*(vertex_count() + it + 1));
}

template <std::size_t r>
vertex_id_t multirelational_directed_graph<r>::vertex_id(const_vertex_iterator it) const {
  return it - static_cast<const_vertex_iterator>(vertices.get());
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::target_of(const_edge_iterator it) const {
  return static_cast<const_vertex_iterator>(*it);
}

template <std::size_t r>
edge_id_t multirelational_directed_graph<r>::edge_id(const_edge_iterator it) const {
  return it - static_cast<const_edge_iterator>(edges.get());
}

template <std::size_t r>
bool multirelational_directed_graph<r>::edge_exists(const_edge_iterator it) const {
  return it != in_edges_end();
}

template <std::size_t r>
typename multirelational_directed_graph<r>::connectivity_address_type
multirelational_directed_graph<r>::edge_value(const_edge_iterator it) const {
  // A non-existent edge has connectivity value 0 (unconnected), matching the
  // contract of graph::edge_value / directed_graph::edge_value. Callers such
  // as vcp::edge_value() pipe out_edge()/in_edge() directly here without a
  // null check. Both return in_edges_end() when no edge is found.
  return it == in_edges_end() ? 0 : edge_values[edge_id(it)];
}

template <std::size_t r>
const_edge_iterator
multirelational_directed_graph<r>::out_edge(const_vertex_iterator source,
                                            const_vertex_iterator target) const {
  const_edge_iterator it(std::find(out_neighbors_begin(source), out_neighbors_end(source), target));
  return it == out_neighbors_end(source) ? in_edges_end() : it;
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::in_edge(const_vertex_iterator source,
                                                               const_vertex_iterator target) const {
  const_edge_iterator it(std::find(in_neighbors_begin(source), in_neighbors_end(source), target));
  return it == in_neighbors_end(source) ? in_edges_end() : it;
}

template <std::size_t r>
bool multirelational_directed_graph<r>::out_edge_exists(const_vertex_iterator source,
                                                        const_vertex_iterator target) const {
  return out_neighbors_end(source) !=
         std::find(out_neighbors_begin(source), out_neighbors_end(source), target);
}

template <std::size_t r>
bool multirelational_directed_graph<r>::in_edge_exists(const_vertex_iterator source,
                                                       const_vertex_iterator target) const {
  return in_neighbors_end(source) !=
         std::find(in_neighbors_begin(source), in_neighbors_end(source), target);
}

/**
 * @brief Write a multirelational directed graph to an output stream.
 *
 * Emits one line per vertex. Each out-neighbor is written as `<id>,<bitmask>`
 * with neighbors separated by spaces. An isolated vertex produces an empty line.
 * The format matches the input accepted by `operator>>`.
 *
 * @tparam r Number of edge relation bits.
 * @param os Output stream.
 * @param g  Graph to write.
 * @return Reference to @p os.
 */
template <std::size_t r>
std::ostream &operator<<(std::ostream &os, multirelational_directed_graph<r> const &g) {
  for (const_vertex_iterator vIt = g.vertices_begin(); vIt != g.vertices_end(); ++vIt) {
    const_edge_iterator const nBegin = g.out_neighbors_begin(vIt);
    const_edge_iterator const nEnd = g.out_neighbors_end(vIt);
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
 * @brief Read a multirelational directed graph from an input stream.
 *
 * Expects one line per vertex with space-separated entries of the form
 * `<neighbor_id>,<bitmask>` representing out-edges. In-edges and their
 * bitmasks are derived automatically by inverting the out-edge list. An
 * empty line represents a vertex with no out-edges.
 *
 * @tparam r Number of edge relation bits.
 * @param is Input stream positioned at the first vertex line.
 * @param g  Graph object to populate; any previous contents are replaced.
 * @return Reference to @p is.
 */
template <std::size_t r>
std::istream &operator>>(std::istream &is, multirelational_directed_graph<r> &g) {
  std::vector<edge_id_t> out_v_temp;
  std::vector<
      std::pair<vertex_id_t, typename multirelational_directed_graph<r>::connectivity_address_type>>
      out_e_temp;
  std::string s1;
  std::size_t edge_id(0);

  while (std::getline(is, s1, '\n')) {
    out_v_temp.push_back(edge_id);
    std::istringstream iss(s1);
    std::string s2;
    while (std::getline(iss, s2, ' ')) {
      std::istringstream iss2(s2);
      std::string neighbor_str;
      std::string value_str;
      getline(iss2, neighbor_str, ',');
      getline(iss2, value_str, ',');
      std::istringstream iss3(neighbor_str);
      std::istringstream iss4(value_str);
      vertex_id_t neighbor;
      typename multirelational_directed_graph<r>::connectivity_address_type value;
      iss3 >> neighbor;
      iss4 >> value;
      out_e_temp.push_back(std::make_pair(neighbor, value));
      ++edge_id;
    }
  }

  g.num_vertices = out_v_temp.size();
  g.num_out_edges = out_e_temp.size();

  g.vertices = std::make_unique<void *[]>(2 * g.vertex_count() + 1);
  g.edges = std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1);
  g.edge_values =
      std::make_unique<typename multirelational_directed_graph<r>::connectivity_address_type[]>(
          g.out_edge_count() + g.in_edge_count());

  for (size_t i = 0; i < g.vertex_count(); ++i) {
    g.vertices[i] = static_cast<void *>(&g.edges[out_v_temp[i]]);
  }
  g.vertices[g.vertex_count()] = static_cast<void *>(&g.edges[g.out_edge_count()]);
  for (size_t i(0); i < g.out_edge_count(); ++i) {
    g.edges[i] = static_cast<void *>(&g.vertices[out_e_temp[i].first]);
    g.edge_values[i] = out_e_temp[i].second;
  }

  out_v_temp = std::vector<edge_id_t>();
  out_e_temp = std::vector<std::pair<
      vertex_id_t, typename multirelational_directed_graph<r>::connectivity_address_type>>();

  std::vector<std::vector<std::pair<
      vertex_id_t, typename multirelational_directed_graph<r>::connectivity_address_type>>>
      in_temp(g.vertex_count());
  for (const_vertex_iterator vIt(g.vertices_begin()); vIt != g.vertices_end(); ++vIt) {
    for (const_edge_iterator eIt(g.out_neighbors_begin(vIt)); eIt != g.out_neighbors_end(vIt);
         ++eIt) {
      in_temp[g.vertex_id(g.target_of(eIt))].push_back(
          std::make_pair(g.vertex_id(vIt), g.edge_value(eIt)));
    }
  }

  std::size_t in_count(0);
  for (size_t i(0); i < g.vertex_count(); ++i) {
    g.vertices[g.vertex_count() + i] = static_cast<void *>(&g.edges[g.out_edge_count() + in_count]);
    std::vector<std::pair<
        vertex_id_t, typename multirelational_directed_graph<r>::connectivity_address_type>> const
        &in_neighbors = in_temp[i];
    for (auto it(in_neighbors.begin()); it != in_neighbors.end(); ++it) {
      g.edges[g.out_edge_count() + in_count] = static_cast<void *>(&g.vertices[it->first]);
      g.edge_values[g.out_edge_count() + in_count] = it->second;
      ++in_count;
    }
  }

  g.vertices[2 * g.vertex_count()] =
      static_cast<void *>(&g.edges[g.out_edge_count() + g.in_edge_count()]);
  g.edges[g.out_edge_count() + g.in_edge_count()] = nullptr;

  return is;
}

} // namespace vcp

#endif

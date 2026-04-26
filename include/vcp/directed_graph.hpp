// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_DIRECTED_GRAPH_HPP
#define VCP_DIRECTED_GRAPH_HPP

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
#include <vcp/graph.hpp>
#include <vector>

namespace vcp {

/**
 * @brief Directed graph with both out-edge and in-edge CSR adjacency lists.
 *
 * Stores out-edges and in-edges in a single flat `edges` array. The
 * out-edge block occupies indices [0, out_edge_count()) and the in-edge
 * block occupies [out_edge_count(), out_edge_count() + in_edge_count()).
 * Each vertex has two entries in the `vertices` array: `vertices[id]`
 * points into the out-edge block and `vertices[vertex_count() + id]` points
 * into the in-edge block, enabling O(1) access to both neighbor lists.
 *
 * The out-edge data is read directly from the stream; in-edges are derived
 * by inverting the out-edge list during `operator>>`. Neighbor lists are
 * expected to be sorted by target vertex id; the `vcp<3, 1, true>` and
 * `vcp<4, 1, true>` specializations rely on this ordering for their
 * merge-based enumeration.
 *
 * Note: `out_edge_count()` and `in_edge_count()` both return the same value
 * (the number of directed out-edges), since every out-edge implies exactly
 * one in-edge.
 */
class directed_graph {
public:
  /** @brief Construct an empty directed graph with no vertices or edges. */
  directed_graph();

  /** @brief Copy-construct a directed graph, deep-copying both CSR arrays. */
  directed_graph(directed_graph const &);

  /** @brief Move-construct in O(1) by transferring ownership of the CSR arrays. */
  directed_graph(directed_graph &&) noexcept = default;

  ~directed_graph();

  /** @brief Copy-assign a directed graph, deep-copying both CSR arrays. */
  directed_graph &operator=(directed_graph const &);

  /** @brief Move-assign in O(1) by transferring ownership of the CSR arrays. */
  directed_graph &operator=(directed_graph &&) noexcept = default;

  /** @brief Return the number of vertices. */
  [[nodiscard]] std::size_t vertex_count() const noexcept;

  /** @brief Return the number of directed out-edges. */
  [[nodiscard]] std::size_t out_edge_count() const noexcept;

  /**
   * @brief Return the number of directed in-edges.
   *
   * Because each out-edge implies exactly one in-edge, this always equals
   * out_edge_count().
   */
  [[nodiscard]] std::size_t in_edge_count() const noexcept;

  /** @brief Return an iterator to the first vertex. */
  [[nodiscard]] const_vertex_iterator vertices_begin() const noexcept;

  /** @brief Return a past-the-end iterator for vertices. */
  [[nodiscard]] const_vertex_iterator vertices_end() const noexcept;

  /** @brief Return an iterator to the first slot in the out-edge block. */
  [[nodiscard]] const_edge_iterator out_edges_begin() const noexcept;

  /** @brief Return a past-the-end iterator for the out-edge block. */
  [[nodiscard]] const_edge_iterator out_edges_end() const noexcept;

  /** @brief Return an iterator to the first slot in the in-edge block. */
  [[nodiscard]] const_edge_iterator in_edges_begin() const noexcept;

  /** @brief Return a past-the-end iterator for the in-edge block. */
  [[nodiscard]] const_edge_iterator in_edges_end() const noexcept;

  /**
   * @brief Return an iterator to the first out-neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first out-neighbor edge slot.
   */
  [[nodiscard]] const_edge_iterator out_neighbors_begin(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return a past-the-end iterator for out-neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the out-neighbor range.
   */
  [[nodiscard]] const_edge_iterator out_neighbors_end(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return an iterator to the first in-neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first in-neighbor edge slot.
   */
  [[nodiscard]] const_edge_iterator in_neighbors_begin(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return a past-the-end iterator for in-neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the in-neighbor range.
   */
  [[nodiscard]] const_edge_iterator in_neighbors_end(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return the integer id of the given vertex.
   *
   * @param it Vertex iterator obtained from this graph.
   * @return Zero-based vertex identifier.
   */
  [[nodiscard]] vertex_id_t vertex_id(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return the vertex pointed to by an edge iterator.
   *
   * @param it Edge iterator within an out- or in-neighbor range.
   * @return Iterator to the target vertex of the edge.
   */
  [[nodiscard]] const_vertex_iterator target_of(const_edge_iterator it) const noexcept;

  /**
   * @brief Return the integer id of the given edge slot.
   *
   * @param it Edge iterator obtained from this graph.
   * @return Zero-based offset of @p it within the combined edge array.
   */
  [[nodiscard]] edge_id_t edge_id(const_edge_iterator it) const noexcept;

  /**
   * @brief Return true if @p it refers to an existing edge (i.e., is not in_edges_end()).
   *
   * @param it Edge iterator to test.
   * @return True if @p it != in_edges_end().
   */
  [[nodiscard]] bool edge_exists(const_edge_iterator it) const noexcept;

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
  [[nodiscard]] const_edge_iterator out_edge(const_vertex_iterator source,
                                             const_vertex_iterator target) const;

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
  [[nodiscard]] const_edge_iterator in_edge(const_vertex_iterator source,
                                            const_vertex_iterator target) const;

  /**
   * @brief Return true if a directed out-edge exists from @p source to @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the out-edge is present.
   */
  [[nodiscard]] bool out_edge_exists(const_vertex_iterator source,
                                     const_vertex_iterator target) const;

  /**
   * @brief Return true if a directed in-edge exists from @p source to @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the in-edge is present.
   */
  [[nodiscard]] bool in_edge_exists(const_vertex_iterator source,
                                    const_vertex_iterator target) const;

  friend std::ostream &operator<<(std::ostream &, directed_graph const &);
  friend std::istream &operator>>(std::istream &, directed_graph &);

private:
  std::size_t num_vertices{0};
  std::size_t num_out_edges{0};
  std::unique_ptr<void *[]> vertices;
  std::unique_ptr<void *[]> edges;
};

inline directed_graph::directed_graph()
    : vertices(std::make_unique<void *[]>(1)), edges(std::make_unique<void *[]>(1)) {
  vertices[0] = static_cast<void *>(&edges[0]);
  edges[0] = nullptr;
}

inline directed_graph::directed_graph(directed_graph const &g)
    : num_vertices(g.num_vertices), num_out_edges(g.num_out_edges),
      vertices(std::make_unique<void *[]>(2 * g.vertex_count() + 1)),
      edges(std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1)) {
  for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
    vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.out_neighbors_begin(it))]);
    vertices[vertex_count() + g.vertex_id(it)] =
        static_cast<void *>(&edges[g.edge_id(g.in_neighbors_begin(it))]);
  }
  vertices[2 * vertex_count()] = static_cast<void *>(&edges[out_edge_count() + in_edge_count()]);
  for (const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
  }
  for (const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it) {
    edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
  }
  edges[out_edge_count() + in_edge_count()] = nullptr;
}

inline directed_graph::~directed_graph() = default;

inline directed_graph &directed_graph::operator=(directed_graph const &g) {
  if (this != &g) {
    num_vertices = g.num_vertices;
    num_out_edges = g.num_out_edges;
    vertices = std::make_unique<void *[]>(2 * g.vertex_count() + 1);
    edges = std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1);
    for (const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it) {
      vertices[g.vertex_id(it)] = static_cast<void *>(&edges[g.edge_id(g.out_neighbors_begin(it))]);
      vertices[vertex_count() + g.vertex_id(it)] =
          static_cast<void *>(&edges[g.edge_id(g.in_neighbors_begin(it))]);
    }
    vertices[2 * vertex_count()] = static_cast<void *>(&edges[out_edge_count() + in_edge_count()]);
    for (const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it) {
      edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    }
    for (const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it) {
      edges[g.edge_id(it)] = static_cast<void *>(&vertices[g.vertex_id(g.target_of(it))]);
    }
    edges[out_edge_count() + in_edge_count()] = nullptr;
  }
  return *this;
}

inline std::size_t directed_graph::vertex_count() const noexcept { return num_vertices; }

inline std::size_t directed_graph::out_edge_count() const noexcept { return num_out_edges; }

inline std::size_t directed_graph::in_edge_count() const noexcept { return num_out_edges; }

inline const_vertex_iterator directed_graph::vertices_begin() const noexcept {
  return &vertices[0];
}

inline const_vertex_iterator directed_graph::vertices_end() const noexcept {
  return &vertices[vertex_count()];
}

inline const_edge_iterator directed_graph::out_edges_begin() const noexcept { return &edges[0]; }

inline const_edge_iterator directed_graph::out_edges_end() const noexcept {
  return &edges[out_edge_count()];
}

inline const_edge_iterator directed_graph::in_edges_begin() const noexcept {
  return &edges[out_edge_count()];
}

inline const_edge_iterator directed_graph::in_edges_end() const noexcept {
  return &edges[out_edge_count() + in_edge_count()];
}

inline const_edge_iterator
directed_graph::out_neighbors_begin(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*it);
}

inline const_edge_iterator
directed_graph::out_neighbors_end(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*(it + 1));
}

inline const_edge_iterator
directed_graph::in_neighbors_begin(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*(vertex_count() + it));
}

inline const_edge_iterator
directed_graph::in_neighbors_end(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*(vertex_count() + it + 1));
}

inline vertex_id_t directed_graph::vertex_id(const_vertex_iterator it) const noexcept {
  // Iterator subtraction yields a signed difference_type; the result is
  // non-negative for any valid `it` in [vertices_begin, vertices_end), so
  // the cast to the unsigned vertex_id_t is value-preserving.
  return static_cast<vertex_id_t>(it - static_cast<const_vertex_iterator>(vertices.get()));
}

inline const_vertex_iterator directed_graph::target_of(const_edge_iterator it) const noexcept {
  return static_cast<const_vertex_iterator>(*it);
}

inline edge_id_t directed_graph::edge_id(const_edge_iterator it) const noexcept {
  // See vertex_id() comment above re: signed-to-unsigned cast.
  return static_cast<edge_id_t>(it - static_cast<const_edge_iterator>(edges.get()));
}

inline bool directed_graph::edge_exists(const_edge_iterator it) const noexcept {
  return it != in_edges_end();
}

inline const_edge_iterator directed_graph::out_edge(const_vertex_iterator source,
                                                    const_vertex_iterator target) const {
  const_edge_iterator it(std::find(out_neighbors_begin(source), out_neighbors_end(source), target));
  return it == out_neighbors_end(source) ? in_edges_end() : it;
}

inline const_edge_iterator directed_graph::in_edge(const_vertex_iterator source,
                                                   const_vertex_iterator target) const {
  const_edge_iterator it(std::find(in_neighbors_begin(source), in_neighbors_end(source), target));
  return it == in_neighbors_end(source) ? in_edges_end() : it;
}

inline bool directed_graph::out_edge_exists(const_vertex_iterator source,
                                            const_vertex_iterator target) const {
  return out_neighbors_end(source) !=
         std::find(out_neighbors_begin(source), out_neighbors_end(source), target);
}

inline bool directed_graph::in_edge_exists(const_vertex_iterator source,
                                           const_vertex_iterator target) const {
  return in_neighbors_end(source) !=
         std::find(in_neighbors_begin(source), in_neighbors_end(source), target);
}

/**
 * @brief Write a directed graph to an output stream.
 *
 * Emits one line per vertex containing space-separated out-neighbor ids.
 * An isolated vertex produces an empty line. The format is the same as
 * that accepted by `operator>>`.
 *
 * @param os Output stream.
 * @param g  Directed graph to write.
 * @return Reference to @p os.
 */
inline std::ostream &operator<<(std::ostream &os, directed_graph const &g) {
  for (const_vertex_iterator vIt = g.vertices_begin(); vIt != g.vertices_end(); ++vIt) {
    const_edge_iterator const nBegin = g.out_neighbors_begin(vIt);
    const_edge_iterator const nEnd = g.out_neighbors_end(vIt);
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
 * @brief Read a directed graph from an input stream.
 *
 * Expects one line per vertex with space-separated out-neighbor ids (0-based).
 * In-edges are derived automatically by inverting the out-edge list. An empty
 * line represents a vertex with no out-edges.
 *
 * @param is Input stream positioned at the first vertex line.
 * @param g  Directed graph object to populate; any previous contents are replaced.
 * @return Reference to @p is.
 * @throws std::invalid_argument If a token cannot be parsed as a vertex id.
 */
inline std::istream &operator>>(std::istream &is, directed_graph &g) {
  std::vector<edge_id_t> out_v_temp;
  std::vector<vertex_id_t> out_e_temp;
  std::string s1;
  std::size_t edge_id(0);

  while (std::getline(is, s1, '\n')) {
    out_v_temp.push_back(edge_id);
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
          throw std::invalid_argument("directed_graph parse: invalid vertex id '" + s2 + "'");
        }
      }
      out_e_temp.push_back(parsed);
      ++edge_id;
    }
  }

  g.num_vertices = out_v_temp.size();
  g.num_out_edges = out_e_temp.size();

  g.vertices = std::make_unique<void *[]>(2 * g.vertex_count() + 1);
  g.edges = std::make_unique<void *[]>(g.out_edge_count() + g.in_edge_count() + 1);

  for (std::size_t i = 0; i < g.vertex_count(); ++i) {
    g.vertices[i] = static_cast<void *>(&g.edges[out_v_temp.at(i)]);
  }
  g.vertices[g.vertex_count()] = static_cast<void *>(&g.edges[g.out_edge_count()]);
  for (std::size_t i = 0; i < g.out_edge_count(); ++i) {
    g.edges[i] = static_cast<void *>(&g.vertices[out_e_temp.at(i)]);
  }

  out_v_temp = std::vector<edge_id_t>();
  out_e_temp = std::vector<vertex_id_t>();

  std::vector<std::vector<vertex_id_t>> in_temp(g.vertex_count());
  for (const_vertex_iterator vIt(g.vertices_begin()); vIt != g.vertices_end(); ++vIt) {
    for (const_edge_iterator eIt(g.out_neighbors_begin(vIt)); eIt != g.out_neighbors_end(vIt);
         ++eIt) {
      in_temp[g.vertex_id(g.target_of(eIt))].push_back(g.vertex_id(vIt));
    }
  }

  std::size_t in_count(0);
  for (std::size_t i = 0; i < g.vertex_count(); ++i) {
    g.vertices[g.vertex_count() + i] = static_cast<void *>(&g.edges[g.out_edge_count() + in_count]);
    std::vector<vertex_id_t> const &in_neighbors = in_temp[i];
    for (unsigned long in_neighbor : in_neighbors) {
      g.edges[g.out_edge_count() + in_count++] = static_cast<void *>(&g.vertices[in_neighbor]);
    }
  }

  g.vertices[2 * g.vertex_count()] =
      static_cast<void *>(&g.edges[g.out_edge_count() + g.in_edge_count()]);
  g.edges[g.out_edge_count() + g.in_edge_count()] = nullptr;

  return is;
}

} // namespace vcp

#endif

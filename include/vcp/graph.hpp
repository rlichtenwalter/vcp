// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_GRAPH_HPP
#define VCP_GRAPH_HPP

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

  /** @brief Move-construct a graph in O(1) by transferring ownership of the CSR arrays. */
  graph(graph &&) noexcept = default;

  /** @brief Copy-assign a graph, deep-copying both CSR arrays. */
  graph &operator=(graph const &);

  /** @brief Move-assign a graph in O(1) by transferring ownership of the CSR arrays. */
  graph &operator=(graph &&) noexcept = default;

  ~graph();

  /** @brief Return the number of vertices. */
  [[nodiscard]] std::size_t vertex_count() const noexcept;

  /** @brief Return the number of undirected edges (each pair counted once). */
  [[nodiscard]] std::size_t edge_count() const noexcept;

  /** @brief Return an iterator to the first vertex. */
  [[nodiscard]] const_vertex_iterator vertices_begin() const noexcept;

  /** @brief Return a past-the-end iterator for vertices. */
  [[nodiscard]] const_vertex_iterator vertices_end() const noexcept;

  /** @brief Return an iterator to the first edge slot in the CSR edge array. */
  [[nodiscard]] const_edge_iterator edges_begin() const noexcept;

  /** @brief Return a past-the-end iterator for the CSR edge array. */
  [[nodiscard]] const_edge_iterator edges_end() const noexcept;

  /**
   * @brief Return an iterator to the first neighbor of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Iterator to the first neighbor edge slot.
   */
  [[nodiscard]] const_edge_iterator neighbors_begin(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return a past-the-end iterator for the neighbors of the given vertex.
   *
   * @param it Iterator to a vertex in this graph.
   * @return Past-the-end iterator for the neighbor range.
   */
  [[nodiscard]] const_edge_iterator neighbors_end(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return the integer id of the given vertex.
   *
   * @param it Vertex iterator obtained from this graph.
   * @return Zero-based vertex identifier equal to the iterator's offset from vertices_begin().
   */
  [[nodiscard]] vertex_id_t vertex_id(const_vertex_iterator it) const noexcept;

  /**
   * @brief Return the vertex pointed to by an edge iterator.
   *
   * @param it Edge iterator within a neighbor range.
   * @return Iterator to the target vertex of the edge.
   */
  [[nodiscard]] const_vertex_iterator target_of(const_edge_iterator it) const noexcept;

  /**
   * @brief Return the integer id of the given edge slot.
   *
   * @param it Edge iterator obtained from this graph.
   * @return Zero-based offset of @p it within the CSR edge array.
   */
  [[nodiscard]] edge_id_t edge_id(const_edge_iterator it) const noexcept;

  /**
   * @brief Return true if @p it refers to an existing edge (i.e., is not edges_end()).
   *
   * @param it Edge iterator to test.
   * @return True if @p it != edges_end().
   */
  [[nodiscard]] bool edge_exists(const_edge_iterator it) const noexcept;

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
  [[nodiscard]] const_edge_iterator edge(const_vertex_iterator source,
                                         const_vertex_iterator target) const;

  /**
   * @brief Return true if an edge exists between @p source and @p target.
   *
   * @param source Iterator to the source vertex.
   * @param target Iterator to the target vertex.
   * @return True if the edge is present.
   */
  [[nodiscard]] bool edge_exists(const_vertex_iterator source, const_vertex_iterator target) const;

  friend std::ostream &operator<<(std::ostream &, graph const &);
  friend std::istream &operator>>(std::istream &, graph &);

private:
  std::size_t num_vertices{0};
  std::size_t num_edges{0};
  std::unique_ptr<void *[]> vertices;
  std::unique_ptr<void *[]> edges;
};

inline graph::graph()
    : vertices(std::make_unique<void *[]>(1)), edges(std::make_unique<void *[]>(1)) {
  vertices[0] = static_cast<void *>(&edges[0]);
  edges[0] = nullptr;
}

inline graph::graph(graph const &g)
    : num_vertices(g.num_vertices), num_edges(g.num_edges),
      vertices(std::make_unique<void *[]>(g.vertex_count() + 1)),
      edges(std::make_unique<void *[]>(g.num_edges + 1)) {
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
    vertices = std::make_unique<void *[]>(g.vertex_count() + 1);
    edges = std::make_unique<void *[]>(g.num_edges + 1);
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

inline std::size_t graph::vertex_count() const noexcept { return num_vertices; }

inline std::size_t graph::edge_count() const noexcept { return num_edges / 2; }

inline const_vertex_iterator graph::vertices_begin() const noexcept { return &vertices[0]; }

inline const_vertex_iterator graph::vertices_end() const noexcept {
  return &vertices[vertex_count()];
}

inline const_edge_iterator graph::edges_begin() const noexcept { return &edges[0]; }

inline const_edge_iterator graph::edges_end() const noexcept { return &edges[num_edges]; }

inline const_edge_iterator graph::neighbors_begin(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*it);
}

inline const_edge_iterator graph::neighbors_end(const_vertex_iterator it) const noexcept {
  return static_cast<const_edge_iterator>(*(it + 1));
}

inline vertex_id_t graph::vertex_id(const_vertex_iterator it) const noexcept {
  return it - static_cast<const_vertex_iterator>(vertices.get());
}

inline const_vertex_iterator graph::target_of(const_edge_iterator it) const noexcept {
  return static_cast<const_vertex_iterator>(*it);
}

inline edge_id_t graph::edge_id(const_edge_iterator it) const noexcept {
  return it - static_cast<const_edge_iterator>(edges.get());
}

inline bool graph::edge_exists(const_edge_iterator it) const noexcept { return it != edges_end(); }

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

  g.vertices = std::make_unique<void *[]>(g.vertex_count() + 1);
  g.edges = std::make_unique<void *[]>(g.num_edges + 1);

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

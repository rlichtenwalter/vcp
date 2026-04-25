// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_HPP
#define VCP_VCP_HPP

#include <cstddef>
#include <map>
#include <type_traits>
#include <utility>
#include <vcp/directed_graph.hpp>
#include <vcp/graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_3_1_0.hpp>
#include <vcp/vcp_3_1_1.hpp>
#include <vcp/vcp_3_r_0.hpp>
#include <vcp/vcp_3_r_1.hpp>
#include <vcp/vcp_4_1_0.hpp>
#include <vcp/vcp_4_1_1.hpp>
#include <vcp/vcp_4_r_0.hpp>
#include <vcp/vcp_4_r_1.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {

/**
 * @brief Generic Vertex Collocation Profile (VCP) calculator.
 *
 * Computes the VCP feature vector for a pivot vertex pair (v1, v2) in a
 * graph. The VCP counts how many n-vertex subgraphs containing v1 and v2
 * exhibit each possible isomorphism class of internal connectivity, as
 * defined in Lichtenwalter & Chawla (2012, 2014).
 *
 * This primary template handles arbitrary (n, r, d) combinations using
 * `vcp_dynamic_mapper` for canonicalization. For (n=3, r=1), (n=4, r=1),
 * and (n=3, r>1), (n=4, r>1), specialized implementations provide
 * faster execution and are selected automatically by the compiler.
 *
 * Thread safety: this template is NOT safe for shared-instance concurrent
 * calls because `generate_vector` uses mutable intermediate state. Construct
 * one `vcp` instance per thread. The underlying graph may be shared for read.
 *
 * @tparam n Number of vertices in each enumerated subgraph (pivot pair + n-2 others).
 * @tparam r Number of edge relations.
 * @tparam d True for directed graphs; false for undirected.
 */
template <std::size_t n, std::size_t r, bool d> class vcp {
public:
  /**
   * @brief Graph type selected automatically based on the (r, d) parameters.
   *
   * Resolves to `directed_graph` (r=1, d=true), `graph` (r=1, d=false),
   * `multirelational_directed_graph<r>` (r>1, d=true), or
   * `multirelational_graph<r>` (r>1, d=false).
   */
  using graph_type = typename std::conditional<
      d,
      typename std::conditional<(r > 1), multirelational_directed_graph<r>, directed_graph>::type,
      typename std::conditional<(r > 1), multirelational_graph<r>, graph>::type>::type;

  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   *
   * Alias of `multirelational_graph<r>::connectivity_address_type`.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete n-vertex subgraph.
   *
   * Alias of `vcp_dynamic_mapper<n, r, d>::subgraph_address_type`.
   */
  using subgraph_address_type = typename vcp_dynamic_mapper<n, r, d>::subgraph_address_type;

  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * The graph must outlive this object. The constructor does not copy the
   * graph; it stores a const reference.
   *
   * @param g Graph to analyze.
   */
  vcp(graph_type const &g);

  /**
   * @brief Compute the VCP feature vector for the pivot pair (v1, v2).
   *
   * Enumerates all n-vertex subgraphs that include v1 and v2, classifies each
   * by its canonical element address, and returns a sparse map from element
   * address to count. Entries with count 0 are omitted. The sum of all counts
   * equals C(|V|-2, n-2), where |V| is the graph's vertex count.
   *
   * @param v1 Iterator to the first pivot vertex.
   * @param v2 Iterator to the second pivot vertex.
   * @return Sparse map from canonical subgraph address to occurrence count.
   */
  std::map<subgraph_address_type, unsigned long> const generate_vector(const_vertex_iterator v1,
                                                                       const_vertex_iterator v2);

private:
  graph_type const &g;
  vcp_dynamic_mapper<n, r, d> mapper;
  void helper(std::array<const_vertex_iterator, n> &vertices, std::size_t current_index,
              square_matrix<connectivity_address_type, n> &connectivity,
              std::map<subgraph_address_type, unsigned long> &counts);
};

const_edge_iterator edge(graph const &g, const_vertex_iterator v1, const_vertex_iterator v2);
const_edge_iterator edge(directed_graph const &g, const_vertex_iterator v1,
                         const_vertex_iterator v2);
template <std::size_t r>
const_edge_iterator edge(multirelational_graph<r> const &g, const_vertex_iterator v1,
                         const_vertex_iterator v2);
template <std::size_t r>
const_edge_iterator edge(multirelational_directed_graph<r> const &g, const_vertex_iterator v1,
                         const_vertex_iterator v2);
bool edge_value(graph const &g, const_edge_iterator edge);
bool edge_value(directed_graph const &g, const_edge_iterator edge);
template <std::size_t r>
typename multirelational_graph<r>::connectivity_address_type
edge_value(multirelational_graph<r> const &g, const_edge_iterator edge);
template <std::size_t r>
typename multirelational_graph<r>::connectivity_address_type
edge_value(multirelational_directed_graph<r> const &g, const_edge_iterator edge);

inline const_edge_iterator edge(graph const &g, const_vertex_iterator v1,
                                const_vertex_iterator v2) {
  return g.edge(v1, v2);
}

inline const_edge_iterator edge(directed_graph const &g, const_vertex_iterator v1,
                                const_vertex_iterator v2) {
  return g.out_edge(v1, v2);
}

template <std::size_t r>
const_edge_iterator edge(multirelational_graph<r> const &g, const_vertex_iterator v1,
                         const_vertex_iterator v2) {
  return g.edge(v1, v2);
}
template <std::size_t r>
const_edge_iterator edge(multirelational_directed_graph<r> const &g, const_vertex_iterator v1,
                         const_vertex_iterator v2) {
  return g.out_edge(v1, v2);
}

inline bool edge_value(graph const &g, const_edge_iterator edge) { return g.edge_exists(edge); }

inline bool edge_value(directed_graph const &g, const_edge_iterator edge) {
  return g.edge_exists(edge);
}

template <typename std::size_t r>
typename multirelational_graph<r>::connectivity_address_type
edge_value(multirelational_graph<r> const &g, const_edge_iterator edge) {
  return g.edge_value(edge);
}

template <typename std::size_t r>
typename multirelational_graph<r>::connectivity_address_type
edge_value(multirelational_directed_graph<r> const &g, const_edge_iterator edge) {
  return g.edge_value(edge);
}

template <std::size_t n, std::size_t r, bool d> vcp<n, r, d>::vcp(graph_type const &g) : g(g) {}

template <std::size_t n, std::size_t r, bool d>
void vcp<n, r, d>::helper(std::array<const_vertex_iterator, n> &vertices, std::size_t current_index,
                          square_matrix<connectivity_address_type, n> &connectivity,
                          std::map<subgraph_address_type, unsigned long> &counts) {
  if (vertices[current_index] == vertices[0] || vertices[current_index] == vertices[1]) {
    return;
  }
  for (std::size_t index(0); index < current_index; ++index) {
    connectivity(index, current_index) =
        edge_value(g, edge(g, vertices[index], vertices[current_index]));
    if (d) {
      connectivity(current_index, index) =
          edge_value(g, edge(g, vertices[current_index], vertices[index]));
    }
  }
  if (n == current_index + 1) {
    ++counts.insert(std::make_pair(mapper.canonical_subgraph_address(connectivity), 0))
          .first->second;
  } else {
    for (const_vertex_iterator v(vertices[current_index] + 1); v != g.vertices_end(); ++v) {
      vertices[current_index + 1] = v;
      helper(vertices, current_index + 1, connectivity, counts);
    }
  }
  return;
}

template <std::size_t n, std::size_t r, bool d>
std::map<typename vcp<n, r, d>::subgraph_address_type, unsigned long> const
vcp<n, r, d>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::map<subgraph_address_type, unsigned long> counts;
  std::array<const_vertex_iterator, n> vertices;
  square_matrix<connectivity_address_type, n> connectivity;
  connectivity(0, 1) = edge_value(g, edge(g, v1, v2));
  if (d) {
    connectivity(1, 0) = edge_value(g, edge(g, v2, v1));
  }
  vertices[0] = v1;
  vertices[1] = v2;
  for (const_vertex_iterator v3(g.vertices_begin()); v3 != g.vertices_end(); ++v3) {
    vertices[2] = v3;
    helper(vertices, 2, connectivity, counts);
  }
  return counts;
}

} // namespace vcp

#endif

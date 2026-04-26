// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_3_R_0_HPP
#define VCP_VCP_3_R_0_HPP

#include <cstddef>
#include <map>
#include <utility>
#include <vcp/graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {

template <std::size_t n, std::size_t r, bool d>
  requires(n >= 2 && r >= 1)
class vcp;

// Safe for shared-instance concurrent calls: this specialization holds
// no mutable per-call state (all state in `generate_vector` is stack-local,
// and the sole class member is a const reference to the graph). The
// vcp<4, ...> specializations do NOT have this property because they carry
// a scratch buffer. See README "Thread safety" for the supported pattern
// and the per-specialization contract.

/**
 * @brief Partial specialization of vcp for 3-vertex subgraphs on a multirelational undirected
 * graph.
 *
 * Handles any number of edge relations r >= 1 using a sorted-merge approach
 * over undirected neighbor lists. The result is a sparse map from canonical
 * subgraph address to count; sparsity increases with r because most of the
 * 2^(3r) possible edge-relation combinations appear rarely or never.
 *
 * This specialization is safe for shared-instance concurrent calls because
 * `generate_vector` holds no mutable class state.
 *
 * @tparam r Number of edge relations (bits per undirected edge).
 */
template <std::size_t r> class vcp<3, r, false> {
public:
  /**
   * @brief Graph type: `multirelational_graph<r>` for r > 1, or `graph` for r == 1.
   */
  using graph_type = typename std::conditional<(r > 1), multirelational_graph<r>, graph>::type;

  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete 3-vertex subgraph address.
   */
  using subgraph_address_type = typename vcp_dynamic_mapper<3, r, false>::subgraph_address_type;

  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * @param g Undirected graph to analyze; must outlive this object.
   */
  vcp(graph_type const &g);

  /**
   * @brief Compute the VCP feature vector for the pivot pair (v1, v2).
   *
   * Returns a sparse map from canonical subgraph address to count. The sum
   * of all counts equals |V| - 2.
   *
   * @param v1 Iterator to the first pivot vertex.
   * @param v2 Iterator to the second pivot vertex.
   * @return Sparse map from canonical subgraph address to occurrence count.
   */
  [[nodiscard]] std::map<subgraph_address_type, unsigned long>
  generate_vector(const_vertex_iterator v1, const_vertex_iterator v2);

private:
  enum connectivity_value : std::size_t { // in bit shifting terms
    V1V2 = 0 * r,
    V1V3 = 1 * r,
    V2V3 = 2 * r
  };
  graph_type const &g;
};

template <std::size_t r> vcp<3, r, false>::vcp(graph_type const &g) : g(g) {}

template <std::size_t r>
std::map<typename vcp<3, r, false>::subgraph_address_type, unsigned long>
vcp<3, r, false>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::map<subgraph_address_type, unsigned long> counts;

  subgraph_address_type v1v2(subgraph_address_type(g.edge_value(g.edge(v1, v2)))
                             << static_cast<std::size_t>(V1V2));

  const_edge_iterator v1_it(g.neighbors_begin(v1));
  const_edge_iterator v1_end(g.neighbors_end(v1));
  const_edge_iterator v2_it(g.neighbors_begin(v2));
  const_edge_iterator v2_end(g.neighbors_end(v2));

  unsigned long union_cardinality(0);
  while (v1_it != v1_end && v2_it != v2_end) {
    if (g.target_of(v1_it) == v2) {
      ++v1_it;
    } else if (g.target_of(v2_it) == v1) {
      ++v2_it;
    } else {
      ++union_cardinality;
      if (g.target_of(v1_it) < g.target_of(v2_it)) {
        ++counts
              .insert(std::make_pair(v1v2 + (subgraph_address_type(g.edge_value(v1_it))
                                             << static_cast<std::size_t>(V1V3)),
                                     0))
              .first->second;
        ++v1_it;
      } else if (g.target_of(v1_it) > g.target_of(v2_it)) {
        ++counts
              .insert(std::make_pair(v1v2 + (subgraph_address_type(g.edge_value(v2_it))
                                             << static_cast<std::size_t>(V2V3)),
                                     0))
              .first->second;
        ++v2_it;
      } else {
        ++counts
              .insert(std::make_pair(v1v2 +
                                         (subgraph_address_type(g.edge_value(v1_it))
                                          << static_cast<std::size_t>(V1V3)) +
                                         (subgraph_address_type(g.edge_value(v2_it))
                                          << static_cast<std::size_t>(V2V3)),
                                     0))
              .first->second;
        ++v1_it;
        ++v2_it;
      }
    }
  }
  while (v1_it != v1_end) {
    if (g.target_of(v1_it) != v2) {
      ++union_cardinality;
      ++counts
            .insert(std::make_pair(v1v2 + (subgraph_address_type(g.edge_value(v1_it))
                                           << static_cast<std::size_t>(V1V3)),
                                   0))
            .first->second;
    }
    ++v1_it;
  }
  while (v2_it != v2_end) {
    if (g.target_of(v2_it) != v1) {
      ++union_cardinality;
      ++counts
            .insert(std::make_pair(v1v2 + (subgraph_address_type(g.edge_value(v2_it))
                                           << static_cast<std::size_t>(V2V3)),
                                   0))
            .first->second;
    }
    ++v2_it;
  }
  counts.insert(std::make_pair(v1v2, g.vertex_count() - 2 - union_cardinality));

  return counts;
}

} // namespace vcp

#endif

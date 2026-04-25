// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_DETAIL_V3_BUFFER_BOUND_HPP
#define VCP_DETAIL_V3_BUFFER_BOUND_HPP

#include <cstddef>
#include <vcp/graph.hpp>

namespace vcp::detail {

/**
 * @brief Sum of the top two neighborhood sizes in an undirected graph.
 *
 * The `vcp<4, r, false>` profilers fill their `v3Vertices` scratch
 * buffer by merging the neighbor lists of the two pivot vertices. The
 * worst-case occupancy for any pivot pair (v1, v2) is bounded above by
 * deg(v1) + deg(v2); the universal upper bound across all pivot pairs
 * is therefore the sum of the two largest degrees in the graph. Top-2
 * is achievable: pick the two highest-degree vertices as the pivot
 * pair and the bound becomes tight.
 *
 * Computed in a single O(V) pass.
 *
 * @tparam UndirectedGraph Any class exposing `vertices_begin`,
 *   `vertices_end`, `neighbors_begin(it)`, and `neighbors_end(it)`.
 */
template <class UndirectedGraph>
std::size_t top_two_neighborhood_size_undirected(UndirectedGraph const &g) {
  std::size_t top1 = 0;
  std::size_t top2 = 0;
  for (const_vertex_iterator it(g.vertices_begin()); it != g.vertices_end(); ++it) {
    std::size_t const deg = static_cast<std::size_t>(g.neighbors_end(it) - g.neighbors_begin(it));
    if (deg >= top1) {
      top2 = top1;
      top1 = deg;
    } else if (deg > top2) {
      top2 = deg;
    }
  }
  return top1 + top2;
}

/**
 * @brief Sum of the top two unique-neighborhood sizes in a directed graph.
 *
 * The `vcp<4, r, true>` profilers fill `v3Vertices` with one entry per
 * distinct neighbor of each pivot vertex. A target appearing in both
 * the out- and in-neighbor lists is enumerated once via
 * `next_union_element`, not twice. The worst-case occupancy for any
 * pivot pair (v1, v2) is therefore
 * |N_out(v1) ∪ N_in(v1)| + |N_out(v2) ∪ N_in(v2)|, and the universal
 * upper bound across all pivot pairs is the sum of the two largest
 * such union sizes.
 *
 * Computed via a two-finger merge of the sorted out- and in-neighbor
 * lists per vertex. Total cost: O(V + E). The looser
 * out_deg(v) + in_deg(v) bound is up to 2× larger on graphs with many
 * mutual edges (e.g. an undirected dataset re-encoded as directed);
 * the O(E) extra work is negligible relative to the subsequent VCP
 * enumeration cost, which is at least O(V² × per-pair merge work).
 *
 * @tparam DirectedGraph Any class exposing `vertices_begin`,
 *   `vertices_end`, `out_neighbors_begin(it)`, `out_neighbors_end(it)`,
 *   `in_neighbors_begin(it)`, `in_neighbors_end(it)`, and `target_of`.
 */
template <class DirectedGraph>
std::size_t top_two_neighborhood_size_directed(DirectedGraph const &g) {
  std::size_t top1 = 0;
  std::size_t top2 = 0;
  for (const_vertex_iterator it(g.vertices_begin()); it != g.vertices_end(); ++it) {
    const_edge_iterator outIt = g.out_neighbors_begin(it);
    const_edge_iterator outEnd = g.out_neighbors_end(it);
    const_edge_iterator inIt = g.in_neighbors_begin(it);
    const_edge_iterator inEnd = g.in_neighbors_end(it);
    std::size_t unique = 0;
    while (outIt != outEnd && inIt != inEnd) {
      if (g.target_of(outIt) < g.target_of(inIt)) {
        ++outIt;
      } else if (g.target_of(outIt) > g.target_of(inIt)) {
        ++inIt;
      } else {
        ++outIt;
        ++inIt;
      }
      ++unique;
    }
    unique += static_cast<std::size_t>(outEnd - outIt) + static_cast<std::size_t>(inEnd - inIt);
    if (unique >= top1) {
      top2 = top1;
      top1 = unique;
    } else if (unique > top2) {
      top2 = unique;
    }
  }
  return top1 + top2;
}

} // namespace vcp::detail

#endif

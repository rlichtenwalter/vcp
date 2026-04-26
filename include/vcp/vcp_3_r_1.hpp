// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_3_R_1_HPP
#define VCP_VCP_3_R_1_HPP

#include <cstddef>
#include <map>
#include <utility>
#include <vcp/directed_graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
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
 * @brief Partial specialization of vcp for 3-vertex subgraphs on a multirelational directed graph.
 *
 * Handles any number of edge relations r >= 1 on directed graphs. Each directed
 * edge pair encodes both the out-direction and in-direction bitmasks separately,
 * so the subgraph address space grows as 2^(6r). The algorithm merges out- and
 * in-neighbor lists to build a sorted union in O(degree) time per pivot pair.
 *
 * This specialization is safe for shared-instance concurrent calls because
 * `generate_vector` holds no mutable class state.
 *
 * @tparam r Number of edge relations (bits per directed edge).
 */
template <std::size_t r> class vcp<3, r, true> {
public:
  /**
   * @brief Graph type: `multirelational_directed_graph<r>` for r > 1, or `directed_graph` for r
   * == 1.
   */
  using graph_type =
      typename std::conditional<(r > 1), multirelational_directed_graph<r>, directed_graph>::type;

  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete 3-vertex directed subgraph address.
   */
  using subgraph_address_type = typename vcp_dynamic_mapper<3, r, true>::subgraph_address_type;

  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * @param g Directed graph to analyze; must outlive this object.
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
  // Bit-shift positions (not mutually-exclusive enum values): the V*V*
  // constants index which vertex pair, the OUT/IN constants index which
  // direction within the pair. They are summed at use sites to form a
  // composite shift amount. Defined as static constexpr (not enum) because
  // C++20 deprecates arithmetic between values of distinct enum types.
  static constexpr std::size_t V1V2 = 0 * r;
  static constexpr std::size_t V1V3 = 2 * r;
  static constexpr std::size_t V2V3 = 4 * r;
  static constexpr std::size_t OUT = 0 * r;
  static constexpr std::size_t IN = 1 * r;
  graph_type const &g;
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
  next_union_element(const_edge_iterator &it1, const_edge_iterator end1, const_edge_iterator &it2,
                     const_edge_iterator end2) const;
};

template <std::size_t r> vcp<3, r, true>::vcp(graph_type const &g) : g(g) {}

template <std::size_t r>
std::pair<const_edge_iterator, std::pair<typename vcp<3, r, true>::connectivity_address_type,
                                         typename vcp<3, r, true>::connectivity_address_type>>
vcp<3, r, true>::next_union_element(
    const_edge_iterator &it1, const_edge_iterator end1, const_edge_iterator &it2,
    const_edge_iterator end2) const { // out-neighbor iterators should always come first
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      temp;
  if (it1 != end1 && it2 != end2) {
    if (g.target_of(it1) < g.target_of(it2)) {
      temp = std::make_pair(it1, std::make_pair(g.edge_value(it1), 0));
      ++it1;
    } else if (g.target_of(it1) > g.target_of(it2)) {
      temp = std::make_pair(it2, std::make_pair(0, g.edge_value(it2)));
      ++it2;
    } else {
      temp = std::make_pair(it1, std::make_pair(g.edge_value(it1), g.edge_value(it2)));
      ++it1;
      ++it2;
    }
  } else if (it1 != end1) {
    temp = std::make_pair(it1, std::make_pair(g.edge_value(it1), 0));
    ++it1;
  } else if (it2 != end2) {
    temp = std::make_pair(it2, std::make_pair(0, g.edge_value(it2)));
    ++it2;
  } else {
    return std::make_pair(
        end2, std::make_pair(0, 0)); // the second element of this pair should never be used
  }
  return temp;
}

template <std::size_t r>
std::map<typename vcp<3, r, true>::subgraph_address_type, unsigned long>
vcp<3, r, true>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::map<subgraph_address_type, unsigned long> counts;

  subgraph_address_type v1v2(
      (subgraph_address_type(g.edge_value(g.out_edge(v1, v2))) << (V1V2 + OUT)) +
      (subgraph_address_type(g.edge_value(g.in_edge(v1, v2))) << (V1V2 + IN)));

  const_edge_iterator v1_out_neighbors_it(g.out_neighbors_begin(v1));
  const_edge_iterator v1_out_neighbors_end(g.out_neighbors_end(v1));
  const_edge_iterator v1_in_neighbors_it(g.in_neighbors_begin(v1));
  const_edge_iterator v1_in_neighbors_end(g.in_neighbors_end(v1));
  const_edge_iterator v2_out_neighbors_it(g.out_neighbors_begin(v2));
  const_edge_iterator v2_out_neighbors_end(g.out_neighbors_end(v2));
  const_edge_iterator v2_in_neighbors_it(g.in_neighbors_begin(v2));
  const_edge_iterator v2_in_neighbors_end(g.in_neighbors_end(v2));
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      min1(next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                              v1_in_neighbors_end));
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      min2(next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                              v2_in_neighbors_end));
  unsigned long union_cardinality(0);
  while (min1.first != v1_in_neighbors_end && min2.first != v2_in_neighbors_end) {
    if (g.target_of(min1.first) < g.target_of(min2.first)) {
      if (g.target_of(min1.first) != v2) {
        ++union_cardinality;
        ++counts
              .insert(
                  std::make_pair(v1v2 + (subgraph_address_type(min1.second.first) << (V1V3 + OUT)) +
                                     (subgraph_address_type(min1.second.second) << (V1V3 + IN)),
                                 0))
              .first->second;
      }
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
    } else if (g.target_of(min1.first) > g.target_of(min2.first)) {
      if (g.target_of(min2.first) != v1) {
        ++union_cardinality;
        ++counts
              .insert(
                  std::make_pair(v1v2 + (subgraph_address_type(min2.second.first) << (V2V3 + OUT)) +
                                     (subgraph_address_type(min2.second.second) << (V2V3 + IN)),
                                 0))
              .first->second;
      }
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    } else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not
             // need to check to exclude it
      ++union_cardinality;
      ++counts
            .insert(
                std::make_pair(v1v2 + ((subgraph_address_type(min1.second.first)) << (V1V3 + OUT)) +
                                   (subgraph_address_type(min1.second.second) << (V1V3 + IN)) +
                                   (subgraph_address_type(min2.second.first) << (V2V3 + OUT)) +
                                   (subgraph_address_type(min2.second.second) << (V2V3 + IN)),
                               0))
            .first->second;
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    }
  }
  while (min1.first != v1_in_neighbors_end) {
    if (g.target_of(min1.first) != v2) {
      ++union_cardinality;
      ++counts
            .insert(std::make_pair(v1v2 +
                                       (subgraph_address_type(min1.second.first) << (V1V3 + OUT)) +
                                       (subgraph_address_type(min1.second.second) << (V1V3 + IN)),
                                   0))
            .first->second;
    }
    min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                              v1_in_neighbors_end);
  }
  while (min2.first != v2_in_neighbors_end) {
    if (g.target_of(min2.first) != v1) {
      ++union_cardinality;
      ++counts
            .insert(std::make_pair(v1v2 +
                                       (subgraph_address_type(min2.second.first) << (V2V3 + OUT)) +
                                       (subgraph_address_type(min2.second.second) << (V2V3 + IN)),
                                   0))
            .first->second;
    }
    min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                              v2_in_neighbors_end);
  }

  counts.insert(std::make_pair(v1v2, g.vertex_count() - 2 - union_cardinality));

  return counts;
}

} // namespace vcp

#endif

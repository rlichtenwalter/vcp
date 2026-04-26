// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_3_1_1_HPP
#define VCP_VCP_3_1_1_HPP

#include <array>
#include <cstddef>
#include <vcp/directed_graph.hpp>

namespace vcp {

template <std::size_t n, std::size_t r, bool d>
  requires(n >= 2 && r >= 1)
class vcp;

/**
 * @brief Full specialization of vcp for 3-vertex subgraphs on a single-relation directed graph.
 *
 * Each directed edge can be absent, outgoing, incoming, or mutual (3 states × 3 pairs
 * = 64 element classes). The algorithm merges out- and in-neighbor lists into a sorted
 * union in O(degree) time per pivot pair.
 *
 * This specialization is safe for shared-instance concurrent calls because
 * `generate_vector` holds no mutable class state.
 */
template <> class vcp<3, 1, true> {
private:
  const static unsigned long num_elements = 64;

public:
  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * @param g Directed graph to analyze; must outlive this object.
   */
  vcp(directed_graph const &g);

  /**
   * @brief Return the number of VCP element classes for this specialization.
   *
   * Always 64 for (n=3, r=1, d=true).
   *
   * @return 64.
   */
  [[nodiscard]] constexpr static std::size_t element_count() noexcept;

  /**
   * @brief Compute the VCP feature vector for the pivot pair (v1, v2).
   *
   * Returns a fixed-size array of 64 counts, one per element class. The sum
   * of all counts equals |V| - 2.
   *
   * @param v1 Iterator to the first pivot vertex.
   * @param v2 Iterator to the second pivot vertex.
   * @return Array of 64 occurrence counts indexed by element address.
   */
  [[nodiscard]] std::array<unsigned long, num_elements> generate_vector(const_vertex_iterator v1,
                                                                        const_vertex_iterator v2);

private:
  enum directedness_value : std::size_t { OUT = 1, IN = 2, BOTH = 3 };
  // C++20 deprecates arithmetic between distinct enum types; the
  // V1V2/V1V3/V2V3 pair-stride constants are summed with directedness_value
  // values at use sites, so they are static constexpr (not enum) to keep
  // one operand of every cross-axis sum a plain integer.
  static constexpr std::size_t V1V2 = 1;
  static constexpr std::size_t V1V3 = 4;
  static constexpr std::size_t V2V3 = 16;
  directed_graph const &g;
  std::pair<const_edge_iterator, directedness_value>
  next_union_element(const_edge_iterator &it1, const_edge_iterator end1, const_edge_iterator &it2,
                     const_edge_iterator end2) const;
};

constexpr std::size_t vcp<3, 1, true>::element_count() noexcept { return num_elements; }

inline vcp<3, 1, true>::vcp(directed_graph const &graph) : g(graph) {}

std::pair<const_edge_iterator, vcp<3, 1, true>::directedness_value> inline vcp<
    3, 1, true>::next_union_element(const_edge_iterator &it1, const_edge_iterator end1,
                                    const_edge_iterator &it2,
                                    const_edge_iterator end2)
    const { // out-neighbor iterators should always come first
  if (it1 != end1 && it2 != end2) {
    if (g.target_of(it1) < g.target_of(it2)) {
      return std::pair{it1++, OUT};
    } else if (g.target_of(it1) > g.target_of(it2)) {
      return std::pair{it2++, IN};
    } else {
      ++it2;
      return std::pair{it1++, BOTH};
    }
  } else if (it1 != end1) {
    return std::pair{it1++, OUT};
  } else if (it2 != end2) {
    return std::pair{it2++, IN};
  } else {
    return std::pair{end2, OUT}; // the second element of this pair should never be used
  }
}

inline std::array<unsigned long, vcp<3, 1, true>::element_count()>
vcp<3, 1, true>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::array<unsigned long, element_count()> counts = {{0}};

  std::size_t v1v2(g.out_edge_exists(v1, v2) * OUT + g.in_edge_exists(v1, v2) * IN);

  const_edge_iterator v1_out_neighbors_it(g.out_neighbors_begin(v1));
  const_edge_iterator v1_out_neighbors_end(g.out_neighbors_end(v1));
  const_edge_iterator v1_in_neighbors_it(g.in_neighbors_begin(v1));
  const_edge_iterator v1_in_neighbors_end(g.in_neighbors_end(v1));
  const_edge_iterator v2_out_neighbors_it(g.out_neighbors_begin(v2));
  const_edge_iterator v2_out_neighbors_end(g.out_neighbors_end(v2));
  const_edge_iterator v2_in_neighbors_it(g.in_neighbors_begin(v2));
  const_edge_iterator v2_in_neighbors_end(g.in_neighbors_end(v2));
  std::pair<const_edge_iterator, directedness_value> min1(next_union_element(
      v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end));
  std::pair<const_edge_iterator, directedness_value> min2(next_union_element(
      v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end));
  unsigned long union_cardinality(0);
  while (min1.first != v1_in_neighbors_end && min2.first != v2_in_neighbors_end) {
    if (g.target_of(min1.first) < g.target_of(min2.first)) {
      if (g.target_of(min1.first) != v2) {
        ++union_cardinality;
        ++counts[v1v2 + min1.second * V1V3];
      }
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
    } else if (g.target_of(min1.first) > g.target_of(min2.first)) {
      if (g.target_of(min2.first) != v1) {
        ++union_cardinality;
        ++counts[v1v2 + min2.second * V2V3];
      }
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    } else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not
             // need to check to exclude it
      ++union_cardinality;
      ++counts[v1v2 + min1.second * V1V3 + min2.second * V2V3];
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    }
  }
  while (min1.first != v1_in_neighbors_end) {
    if (g.target_of(min1.first) != v2) {
      ++union_cardinality;
      ++counts[v1v2 + min1.second * V1V3];
    }
    min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                              v1_in_neighbors_end);
  }
  while (min2.first != v2_in_neighbors_end) {
    if (g.target_of(min2.first) != v1) {
      ++union_cardinality;
      ++counts[v1v2 + min2.second * V2V3];
    }
    min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                              v2_in_neighbors_end);
  }

  counts[v1v2] = g.vertex_count() - 2 - union_cardinality;

  return counts;
}

} // namespace vcp

#endif

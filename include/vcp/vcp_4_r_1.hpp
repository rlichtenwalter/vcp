// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_4_R_1_HPP
#define VCP_VCP_4_R_1_HPP

#include <cstddef>
#include <map>
#include <utility>
#include <vcp/detail/dense_or_sparse_map.hpp>
#include <vcp/detail/v3_buffer_bound.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

#ifdef VCP_INSTRUMENT_K
#include <k_probe.hpp>
#endif

namespace vcp {

template <std::size_t n, std::size_t r, bool d>
  requires(n >= 2 && r >= 1)
class vcp;

/**
 * @brief Partial specialization of vcp for 4-vertex subgraphs on a multirelational directed graph.
 *
 * Extends the `vcp<4, 1, true>` two-phase algorithm to multirelational directed
 * graphs. The edge-type histogram is keyed by sorted (lo, hi) pairs of r-bit
 * bitmasks representing the two directed directions between a vertex pair. A
 * two-tier `detail::dense_or_sparse_map` handles both small and large r efficiently.
 *
 * This specialization is NOT safe for shared-instance concurrent calls because
 * `generate_vector` uses the `v3Vertices` scratch buffer and `temp_edge_types`
 * map held as class members. Construct one instance per thread.
 *
 * @tparam r Number of edge relations.
 */
template <std::size_t r> class vcp<4, r, true> {
public:
  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   */
  using connectivity_address_type =
      typename multirelational_directed_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete 4-vertex directed subgraph address.
   */
  using subgraph_address_type = typename vcp_dynamic_mapper<4, r, true>::subgraph_address_type;

  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * Builds the global directed edge-type histogram and allocates the v3 scratch buffer.
   *
   * @param g Multirelational directed graph to analyze; must outlive this object.
   */
  vcp(multirelational_directed_graph<r> const &g);

  /**
   * @brief Compute the VCP feature vector for the pivot pair (v1, v2).
   *
   * Returns a sparse map from canonical subgraph address to count. The sum
   * of all counts equals C(|V|-2, 2).
   *
   * @param v1 Iterator to the first pivot vertex.
   * @param v2 Iterator to the second pivot vertex.
   * @return Sparse map from canonical subgraph address to occurrence count.
   */
  [[nodiscard]] std::map<subgraph_address_type, unsigned long>
  generate_vector(const_vertex_iterator v1, const_vertex_iterator v2);

private:
  using connectivity_matrix = square_matrix<connectivity_address_type, 4>;
  // Edge-types map: key is a sorted (lo, hi) pair of r-bit connectivity
  // addresses (the directed edge's out/in classes), so the effective
  // bit count is 2r and the total key space is 2^(2r). For small r
  // (r ≤ 10 at sizeof(unsigned long)=8) the dense tier fires via
  // pack_pair_key_by_bits<r> packing; for larger r the sparse tier
  // uses std::unordered_map keyed on the raw std::pair via pair_hash.
  // See include/vcp/detail/dense_or_sparse_map.hpp.
  using edge_key_type = std::pair<connectivity_address_type, connectivity_address_type>;
  using edge_types_map = detail::dense_or_sparse_map<edge_key_type, unsigned long,
                                                     detail::key_count_from_bits<2 * r>(),
                                                     detail::pack_pair_key_by_bits<r>>;
  multirelational_directed_graph<r> const &g;
  vcp_dynamic_mapper<4, r, true> mapper;
  edge_types_map edge_types;
  // temp_edge_types promoted to class member (see vcp_4_r_0.hpp for
  // rationale). Cleared at the start of each generate_vector.
  edge_types_map temp_edge_types;
  // Scratch buffer for v3 candidates accumulated during `generate_vector`.
  // Allocated once per instance, reused across every call — which makes
  // this class NOT thread-safe for shared-instance concurrent calls.
  // Construct one `vcp<4, r, true>` per thread if parallelism is needed;
  // the underlying graph is safe to share for read. See README "Thread
  // safety" for the supported pattern.
  std::unique_ptr<std::pair<const_vertex_iterator, connectivity_matrix>[]> v3Vertices;
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
  next_union_element(const_edge_iterator &, const_edge_iterator, const_edge_iterator &,
                     const_edge_iterator) const;
};

template <std::size_t r>
vcp<4, r, true>::vcp(multirelational_directed_graph<r> const &graph)
    : g(graph), mapper(),
      v3Vertices(std::make_unique<std::pair<const_vertex_iterator, connectivity_matrix>[]>(
          detail::top_two_neighborhood_size_directed(graph))) {
  // (0, 0) is the unconnected class — seed its count with the number
  // of unordered pairs and decrement per real edge, same pattern as
  // vcp_4_r_0's undirected ctor.
  unsigned long &gaps(edge_types.insert_or_zero(
      std::pair{connectivity_address_type(0), connectivity_address_type(0)}));
  gaps = graph.vertex_count() * (graph.vertex_count() - 1) / 2;
  for (const_vertex_iterator it(graph.vertices_begin()); it != graph.vertices_end(); ++it) {
    const_edge_iterator outIt = graph.out_neighbors_begin(it);
    const_edge_iterator outEnd = graph.out_neighbors_end(it);
    const_edge_iterator inIt = graph.in_neighbors_begin(it);
    const_edge_iterator inEnd = graph.in_neighbors_end(it);
    while (outIt != outEnd && graph.target_of(outIt) <= it) {
      ++outIt;
    }
    while (inIt != inEnd && graph.target_of(inIt) <= it) {
      ++inIt;
    }
    while (outIt != outEnd && inIt != inEnd) {
      if (graph.target_of(outIt) < graph.target_of(inIt)) {
        ++edge_types.insert_or_zero(
            std::pair{connectivity_address_type(0), graph.edge_value(outIt)});
        --gaps;
        ++outIt;
      } else if (graph.target_of(outIt) > graph.target_of(inIt)) {
        ++edge_types.insert_or_zero(
            std::pair{connectivity_address_type(0), graph.edge_value(inIt)});
        --gaps;
        ++inIt;
      } else {
        ++edge_types.insert_or_zero(
            graph.edge_value(outIt) < graph.edge_value(inIt)
                ? std::pair{graph.edge_value(outIt), graph.edge_value(inIt)}
                : std::pair{graph.edge_value(inIt), graph.edge_value(outIt)});
        --gaps;
        ++outIt;
        ++inIt;
      }
    }
    while (outIt != outEnd) {
      ++edge_types.insert_or_zero(std::pair{connectivity_address_type(0), graph.edge_value(outIt)});
      --gaps;
      ++outIt;
    }
    while (inIt != inEnd) {
      ++edge_types.insert_or_zero(std::pair{connectivity_address_type(0), graph.edge_value(inIt)});
      --gaps;
      ++inIt;
    }
  }
}

template <std::size_t r>
std::pair<const_edge_iterator, std::pair<typename vcp<4, r, true>::connectivity_address_type,
                                         typename vcp<4, r, true>::connectivity_address_type>>
vcp<4, r, true>::next_union_element(const_edge_iterator &it1, const_edge_iterator end1,
                                    const_edge_iterator &it2, const_edge_iterator end2) const {
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      temp;
  if (it1 != end1 && it2 != end2) {
    if (g.target_of(it1) < g.target_of(it2)) {
      temp = std::pair{it1, std::pair{g.edge_value(it1), 0}};
      ++it1;
    } else if (g.target_of(it1) > g.target_of(it2)) {
      temp = std::pair{it2, std::pair{0, g.edge_value(it2)}};
      ++it2;
    } else {
      temp = std::pair{it1, std::pair{g.edge_value(it1), g.edge_value(it2)}};
      ++it1;
      ++it2;
    }
  } else if (it1 != end1) {
    temp = std::pair{it1, std::pair{g.edge_value(it1), 0}};
    ++it1;
  } else if (it2 != end2) {
    temp = std::pair{it2, std::pair{0, g.edge_value(it2)}};
    ++it2;
  } else {
    // Both input iterator pairs exhausted. The returned iterator is `end2`
    // by convention — callers compare `min.first` against their own
    // `in_neighbors_end` (which they passed as end2) to detect exhaustion.
    // The v1/v2 value pair is unused; it is (0, 0) only as a defensive
    // default. If this exhaustion convention ever changes, every loop
    // condition of the form `min.first != v*_in_neighbors_end` in
    // `generate_vector` must be updated in lockstep.
    return std::pair{end2, std::pair{0, 0}};
  }
  return temp;
}

template <std::size_t r>
std::map<typename vcp<4, r, true>::subgraph_address_type, unsigned long>
vcp<4, r, true>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::map<subgraph_address_type, unsigned long> counts;
  // temp_edge_types is a class member; see header rationale. O(k) clear.
  temp_edge_types.clear();

  connectivity_matrix connectivity;
  connectivity(0, 1) = g.edge_value(g.out_edge(v1, v2));
  connectivity(1, 0) = g.edge_value(g.in_edge(v1, v2));

  unsigned long &gaps(temp_edge_types.insert_or_zero(
      std::pair{connectivity_address_type(0), connectivity_address_type(0)}));

  // compose ordered list of v3 candidates
  const_edge_iterator v1_out_neighbors_it(g.out_neighbors_begin(v1));
  const_edge_iterator v1_out_neighbors_end(g.out_neighbors_end(v1));
  const_edge_iterator v1_in_neighbors_it(g.in_neighbors_begin(v1));
  const_edge_iterator v1_in_neighbors_end(g.in_neighbors_end(v1));
  const_edge_iterator v2_out_neighbors_it(g.out_neighbors_begin(v2));
  const_edge_iterator v2_out_neighbors_end(g.out_neighbors_end(v2));
  const_edge_iterator v2_in_neighbors_it(g.in_neighbors_begin(v2));
  const_edge_iterator v2_in_neighbors_end(g.in_neighbors_end(v2));
  std::pair<const_vertex_iterator, connectivity_matrix> *v3Vertices_begin(&v3Vertices[0]);
  std::pair<const_vertex_iterator, connectivity_matrix> *v3Vertices_end(&v3Vertices[0]);
  // `next_union_element` returns end2 on exhaustion (see its final else
  // branch). We alias v{1,2}_in_neighbors_end to exhaustion-sentinel names
  // so the loop conditions below read as a contract check rather than an
  // accidental iterator comparison.
  const_edge_iterator const min1_exhausted_sentinel(v1_in_neighbors_end);
  const_edge_iterator const min2_exhausted_sentinel(v2_in_neighbors_end);
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      min1(next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                              v1_in_neighbors_end));
  std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
      min2(next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                              v2_in_neighbors_end));
  while (min1.first != min1_exhausted_sentinel && min2.first != min2_exhausted_sentinel) {
    if (g.target_of(min1.first) < g.target_of(min2.first)) {
      if (g.target_of(min1.first) != v2) {
        ++temp_edge_types.insert_or_zero(min1.second.first < min1.second.second
                                             ? std::pair{min1.second.first, min1.second.second}
                                             : std::pair{min1.second.second, min1.second.first});
        ++gaps;
        v3Vertices_end->first = g.target_of(min1.first);
        v3Vertices_end->second = connectivity;
        v3Vertices_end->second(0, 2) = min1.second.first;
        v3Vertices_end->second(2, 0) = min1.second.second;
        ++v3Vertices_end;
      }
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
    } else if (g.target_of(min1.first) > g.target_of(min2.first)) {
      if (g.target_of(min2.first) != v1) {
        ++temp_edge_types.insert_or_zero(min2.second.first < min2.second.second
                                             ? std::pair{min2.second.first, min2.second.second}
                                             : std::pair{min2.second.second, min2.second.first});
        ++gaps;
        v3Vertices_end->first = g.target_of(min2.first);
        v3Vertices_end->second = connectivity;
        v3Vertices_end->second(1, 2) = min2.second.first;
        v3Vertices_end->second(2, 1) = min2.second.second;
        ++v3Vertices_end;
      }
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    } else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not
             // need to check to exclude it
      ++temp_edge_types.insert_or_zero(min1.second.first < min1.second.second
                                           ? std::pair{min1.second.first, min1.second.second}
                                           : std::pair{min1.second.second, min1.second.first});
      ++temp_edge_types.insert_or_zero(min2.second.first < min2.second.second
                                           ? std::pair{min2.second.first, min2.second.second}
                                           : std::pair{min2.second.second, min2.second.first});
      v3Vertices_end->first = g.target_of(min1.first);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(0, 2) = min1.second.first;
      v3Vertices_end->second(2, 0) = min1.second.second;
      v3Vertices_end->second(1, 2) = min2.second.first;
      v3Vertices_end->second(2, 1) = min2.second.second;
      ++v3Vertices_end;
      min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                                v1_in_neighbors_end);
      min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                                v2_in_neighbors_end);
    }
  }
  while (min1.first != min1_exhausted_sentinel) {
    if (g.target_of(min1.first) != v2) {
      ++temp_edge_types.insert_or_zero(min1.second.first < min1.second.second
                                           ? std::pair{min1.second.first, min1.second.second}
                                           : std::pair{min1.second.second, min1.second.first});
      ++gaps;
      v3Vertices_end->first = g.target_of(min1.first);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(0, 2) = min1.second.first;
      v3Vertices_end->second(2, 0) = min1.second.second;
      ++v3Vertices_end;
    }
    min1 = next_union_element(v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it,
                              v1_in_neighbors_end);
  }
  while (min2.first != min2_exhausted_sentinel) {
    if (g.target_of(min2.first) != v1) {
      ++temp_edge_types.insert_or_zero(min2.second.first < min2.second.second
                                           ? std::pair{min2.second.first, min2.second.second}
                                           : std::pair{min2.second.second, min2.second.first});
      ++gaps;
      v3Vertices_end->first = g.target_of(min2.first);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(1, 2) = min2.second.first;
      v3Vertices_end->second(2, 1) = min2.second.second;
      ++v3Vertices_end;
    }
    min2 = next_union_element(v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it,
                              v2_in_neighbors_end);
  }

  auto v3_count = static_cast<std::size_t>(v3Vertices_end - v3Vertices_begin);
  std::size_t v4_count(0);
  for (std::pair<const_vertex_iterator, connectivity_matrix> *it1(v3Vertices_begin);
       it1 != v3Vertices_end; ++it1) { // for each v3 vertex computed above
    const_edge_iterator v3_out_neighbors_it(g.out_neighbors_begin(it1->first));
    const_edge_iterator v3_out_neighbors_end(g.out_neighbors_end(it1->first));
    const_edge_iterator v3_in_neighbors_it(g.in_neighbors_begin(it1->first));
    const_edge_iterator v3_in_neighbors_end(g.in_neighbors_end(it1->first));
    // See comment on min{1,2}_exhausted_sentinel above.
    const_edge_iterator const min_exhausted_sentinel(v3_in_neighbors_end);
    unsigned long v4_local_count(
        0); // keep track of how many v4 vertices are only the result of the neighbors of this v3
    std::pair<const_edge_iterator, std::pair<connectivity_address_type, connectivity_address_type>>
        min(next_union_element(v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it,
                               v3_in_neighbors_end));
    for (std::pair<const_vertex_iterator, connectivity_matrix> *it2(v3Vertices_begin);
         it2 != v3Vertices_end; ++it2) {
      while (min.first != min_exhausted_sentinel && g.target_of(min.first) < it2->first) {
        if (g.target_of(min.first) != v1 && g.target_of(min.first) != v2) {
          ++temp_edge_types.insert_or_zero(min.second.first < min.second.second
                                               ? std::pair{min.second.first, min.second.second}
                                               : std::pair{min.second.second, min.second.first});
          ++v4_local_count;
          it1->second(0, 3) = 0;
          it1->second(3, 0) = 0;
          it1->second(1, 3) = 0;
          it1->second(3, 1) = 0;
          it1->second(2, 3) = min.second.first;
          it1->second(3, 2) = min.second.second;
          ++counts.insert(std::pair{mapper.canonical_subgraph_address(it1->second), 0})
                .first->second;
        }
        min = next_union_element(v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it,
                                 v3_in_neighbors_end);
      }
      if (min.first == min_exhausted_sentinel || g.target_of(min.first) > it2->first) {
        if (it1->first < it2->first) {
          ++gaps;
          it1->second(0, 3) = it2->second(0, 2);
          it1->second(3, 0) = it2->second(2, 0);
          it1->second(1, 3) = it2->second(1, 2);
          it1->second(3, 1) = it2->second(2, 1);
          it1->second(2, 3) = 0;
          it1->second(3, 2) = 0;
          ++counts.insert(std::pair{mapper.canonical_subgraph_address(it1->second), 0})
                .first->second;
        }
      } else {
        if (it1->first < it2->first) {
          ++temp_edge_types.insert_or_zero(min.second.first < min.second.second
                                               ? std::pair{min.second.first, min.second.second}
                                               : std::pair{min.second.second, min.second.first});
          it1->second(0, 3) = it2->second(0, 2);
          it1->second(3, 0) = it2->second(2, 0);
          it1->second(1, 3) = it2->second(1, 2);
          it1->second(3, 1) = it2->second(2, 1);
          it1->second(2, 3) = min.second.first;
          it1->second(3, 2) = min.second.second;
          ++counts.insert(std::pair{mapper.canonical_subgraph_address(it1->second), 0})
                .first->second;
        }
        min = next_union_element(v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it,
                                 v3_in_neighbors_end);
      }
    }
    while (min.first != min_exhausted_sentinel) {
      if (g.target_of(min.first) != v1 && g.target_of(min.first) != v2) {
        ++temp_edge_types.insert_or_zero(min.second.first < min.second.second
                                             ? std::pair{min.second.first, min.second.second}
                                             : std::pair{min.second.second, min.second.first});
        ++v4_local_count;
        it1->second(0, 3) = 0;
        it1->second(3, 0) = 0;
        it1->second(1, 3) = 0;
        it1->second(3, 1) = 0;
        it1->second(2, 3) = min.second.first;
        it1->second(3, 2) = min.second.second;
        ++counts.insert(std::pair{mapper.canonical_subgraph_address(it1->second), 0}).first->second;
      }
      min = next_union_element(v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it,
                               v3_in_neighbors_end);
    }
    v4_count += v4_local_count;
    gaps += 2 * v4_local_count;
    it1->second(0, 3) = 0;
    it1->second(3, 0) = 0;
    it1->second(1, 3) = 0;
    it1->second(3, 1) = 0;
    it1->second(2, 3) = 0;
    it1->second(3, 2) = 0;
    counts.insert(std::pair{mapper.canonical_subgraph_address(it1->second), 0}).first->second +=
        g.vertex_count() - 2 -
        (static_cast<unsigned long>(v3Vertices_end - v3Vertices_begin) + v4_local_count);
  }

  // v1v2's own directionality class, sorted (lo, hi) to match edge_types keys.
  connectivity_address_type const v1v2_lo(std::min(connectivity(0, 1), connectivity(1, 0)));
  connectivity_address_type const v1v2_hi(std::max(connectivity(0, 1), connectivity(1, 0)));
  edge_types.for_each([&](edge_key_type const &key, unsigned long edge_count) {
    connectivity(2, 3) = key.first;
    connectivity(3, 2) = key.second;
    unsigned long const *temp_val = temp_edge_types.find(key);
    unsigned long count(edge_count);
    if (temp_val != nullptr) {
      count -= *temp_val;
      if (key.first == 0 && key.second == 0) {
        // For the unconnected class, !bool(c01+c10) is 1 iff v1v2 is itself
        // unconnected — the self-correction for that class. The trailing
        // terms subtract pairs that can never be enumerated (both endpoints
        // in V \ ({v1,v2} ∪ Γ(v1) ∪ Γ(v2))).
        count -= !static_cast<bool>(connectivity(0, 1) + connectivity(1, 0)) +
                 (2 + v3_count) * (g.vertex_count() - 2 - v3_count) - 3 * v4_count;
      }
    }
    // Self-correction for non-unconnected classes: edge_types counts v1v2
    // itself if it has any directedness, so when the current class matches
    // v1v2's own class, subtract 1. The unconnected class handles its own
    // self-correction above via !bool(c01+c10).
    if ((key.first != 0 || key.second != 0) && key.first == v1v2_lo && key.second == v1v2_hi) {
      count -= 1;
    }
    counts.insert(std::pair{mapper.canonical_subgraph_address(connectivity), 0}).first->second +=
        count;
  });

#ifdef VCP_INSTRUMENT_K
  k_probe::record("temp_edge_types", temp_edge_types.size());
  k_probe::record("edge_types", edge_types.size());
  k_probe::record("counts", counts.size());
#endif

  return counts;
}

} // namespace vcp

#endif

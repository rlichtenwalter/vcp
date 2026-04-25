// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_4_R_0_HPP
#define VCP_VCP_4_R_0_HPP

#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vcp/detail/dense_or_sparse_map.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

#ifdef VCP_INSTRUMENT_K
#include <k_probe.hpp>
#endif

namespace vcp {

template <std::size_t n, std::size_t r, bool d> class vcp;

/**
 * @brief Partial specialization of vcp for 4-vertex subgraphs on a multirelational undirected
 *        graph.
 *
 * Extends the `vcp<4, 1, false>` two-phase algorithm to multirelational graphs.
 * The constructor builds a global edge-type histogram keyed by r-bit relation
 * bitmask; `generate_vector` uses it to fill the least-connected element classes
 * by subtraction. An internal `vcp_dynamic_mapper` canonicalizes subgraph addresses.
 *
 * The `edge_types` and `temp_edge_types` maps use a two-tier
 * `detail::dense_or_sparse_map` for performance: a heap-backed array for small r
 * (r <= 20) and an unordered map for larger r.
 *
 * This specialization is NOT safe for shared-instance concurrent calls because
 * `generate_vector` uses the `v3Vertices` scratch buffer and `temp_edge_types`
 * map held as class members. Construct one instance per thread.
 *
 * @tparam r Number of edge relations.
 */
template <std::size_t r> class vcp<4, r, false> {
public:
  /**
   * @brief Unsigned integer type that holds an r-bit edge connectivity bitmask.
   */
  using connectivity_address_type = typename multirelational_graph<r>::connectivity_address_type;

  /**
   * @brief Unsigned integer type that encodes a complete 4-vertex undirected subgraph address.
   */
  using subgraph_address_type = typename vcp_dynamic_mapper<4, r, false>::subgraph_address_type;

  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * Builds the global edge-type histogram and allocates the v3 scratch buffer.
   *
   * @param g Multirelational undirected graph to analyze; must outlive this object.
   */
  vcp(multirelational_graph<r> const &g);

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
  std::map<subgraph_address_type, unsigned long> const generate_vector(const_vertex_iterator v1,
                                                                       const_vertex_iterator v2);

private:
  using connectivity_matrix = square_matrix<connectivity_address_type, 4>;
  // Edge-types map: key is the r-bit connectivity address, so the
  // total key space is 2^r. The detail::dense_or_sparse_map facade
  // picks std::array for small r (through r=20 at sizeof(unsigned
  // long)=8) and std::unordered_map above that. Replaces the
  // std::map<connectivity_address_type, unsigned long> that used to
  // dominate this hot path's allocator profile. See
  // include/vcp/detail/dense_or_sparse_map.hpp for the tier rationale
  // and benchmark/small_map_study/results/summary_phase_C.md for the
  // numbers behind it.
  using edge_types_map = detail::dense_or_sparse_map<connectivity_address_type, unsigned long,
                                                     detail::key_count_from_bits<r>()>;
  multirelational_graph<r> const &g;
  vcp_dynamic_mapper<4, r, false> mapper;
  edge_types_map edge_types;
  // temp_edge_types: promoted from function-local to class member so
  // the dense_or_sparse_map allocation survives across calls. clear()
  // at the start of each generate_vector is O(k). At the r=20 dense
  // regime, allocating a 1M-slot array per call via mmap would cost
  // ~50 µs × 20k calls = 1 s of pure malloc overhead, hence the move
  // to a reusable class-member instance. This does NOT change the
  // class's thread-safety story: it was already not thread-safe
  // per-instance because of v3Vertices.
  edge_types_map temp_edge_types;
  // Scratch buffer for v3 candidates accumulated during `generate_vector`.
  // Allocated once per instance, reused across every call — which makes
  // this class NOT thread-safe for shared-instance concurrent calls.
  // Construct one `vcp<4, r, false>` per thread if parallelism is needed;
  // the underlying graph is safe to share for read. See README "Thread
  // safety" for the supported pattern.
  std::unique_ptr<std::pair<const_vertex_iterator, connectivity_matrix>[]> v3Vertices;
};

template <std::size_t r>
vcp<4, r, false>::vcp(multirelational_graph<r> const &g)
    : g(g), mapper(),
      v3Vertices(std::make_unique<std::pair<const_vertex_iterator, connectivity_matrix>[]>(
          MAX_NEIGHBORS)) {
  // Edge-type key 0 is the "no relation" class — initialize its count
  // to the number of unordered pairs, then decrement as real edges
  // are added below. Each real edge shifts one pair out of the
  // unconnected class into its actual relation class.
  unsigned long &gaps(edge_types.insert_or_zero(0));
  gaps = g.vertex_count() * (g.vertex_count() - 1) / 2;
  for (const_vertex_iterator it(g.vertices_begin()); it != g.vertices_end(); ++it) {
    for (const_edge_iterator eIt(g.neighbors_begin(it)); eIt != g.neighbors_end(it); ++eIt) {
      if (it < g.target_of(eIt)) {
        ++edge_types.insert_or_zero(g.edge_value(eIt));
        --gaps;
      }
    }
  }
}

template <std::size_t r>
std::map<typename vcp<4, r, false>::subgraph_address_type, unsigned long> const
vcp<4, r, false>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::map<subgraph_address_type, unsigned long> counts;
  // temp_edge_types is a class member (see header comment). Reset its
  // logical contents before use; the O(1)/O(k) clear keeps the hot
  // path fast while reusing the dense-tier storage.
  temp_edge_types.clear();

  connectivity_matrix connectivity;
  connectivity(0, 1) = g.edge_value(g.edge(v1, v2));

  unsigned long &gaps(temp_edge_types.insert_or_zero(0));

  const_edge_iterator v1_neighbors_it(g.neighbors_begin(v1));
  const_edge_iterator v1_neighbors_end(g.neighbors_end(v1));
  const_edge_iterator v2_neighbors_it(g.neighbors_begin(v2));
  const_edge_iterator v2_neighbors_end(g.neighbors_end(v2));
  assert(MAX_NEIGHBORS >
         (v1_neighbors_end - v1_neighbors_it) +
             (v2_neighbors_end -
              v2_neighbors_it)); // this should always be contiguous storage; we can only
                                 // over-allocate by a factor of 2, which is of much lower cost than
                                 // maintaining a doubly-linked list; there exists a strict upper
                                 // bound on the final size of v3Vertices
  std::pair<const_vertex_iterator, connectivity_matrix> *v3Vertices_begin(&v3Vertices[0]);
  std::pair<const_vertex_iterator, connectivity_matrix> *v3Vertices_end(&v3Vertices[0]);
  while (v1_neighbors_it != v1_neighbors_end && v2_neighbors_it != v2_neighbors_end) {
    if (g.target_of(v1_neighbors_it) < g.target_of(v2_neighbors_it)) {
      if (g.target_of(v1_neighbors_it) != v2) {
        ++temp_edge_types.insert_or_zero(g.edge_value(v1_neighbors_it));
        ++gaps;
        v3Vertices_end->first = g.target_of(v1_neighbors_it);
        v3Vertices_end->second = connectivity;
        v3Vertices_end->second(0, 2) = g.edge_value(v1_neighbors_it);
        ++v3Vertices_end;
      }
      ++v1_neighbors_it;
    } else if (g.target_of(v1_neighbors_it) > g.target_of(v2_neighbors_it)) {
      if (g.target_of(v2_neighbors_it) != v1) {
        ++temp_edge_types.insert_or_zero(g.edge_value(v2_neighbors_it));
        ++gaps;
        v3Vertices_end->first = g.target_of(v2_neighbors_it);
        v3Vertices_end->second = connectivity;
        v3Vertices_end->second(1, 2) = g.edge_value(v2_neighbors_it);
        ++v3Vertices_end;
      }
      ++v2_neighbors_it;
    } else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not
             // need to check to exclude it
      ++temp_edge_types.insert_or_zero(g.edge_value(v1_neighbors_it));
      ++temp_edge_types.insert_or_zero(g.edge_value(v2_neighbors_it));
      v3Vertices_end->first = g.target_of(v1_neighbors_it);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(0, 2) = g.edge_value(v1_neighbors_it);
      v3Vertices_end->second(1, 2) = g.edge_value(v2_neighbors_it);
      ++v3Vertices_end;
      ++v1_neighbors_it;
      ++v2_neighbors_it;
    }
  }
  while (v1_neighbors_it != v1_neighbors_end) {
    if (g.target_of(v1_neighbors_it) != v2) {
      ++temp_edge_types.insert_or_zero(g.edge_value(v1_neighbors_it));
      ++gaps;
      v3Vertices_end->first = g.target_of(v1_neighbors_it);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(0, 2) = g.edge_value(v1_neighbors_it);
      ++v3Vertices_end;
    }
    ++v1_neighbors_it;
  }
  while (v2_neighbors_it != v2_neighbors_end) {
    if (g.target_of(v2_neighbors_it) != v1) {
      ++temp_edge_types.insert_or_zero(g.edge_value(v2_neighbors_it));
      ++gaps;
      v3Vertices_end->first = g.target_of(v2_neighbors_it);
      v3Vertices_end->second = connectivity;
      v3Vertices_end->second(1, 2) = g.edge_value(v2_neighbors_it);
      ++v3Vertices_end;
    }
    ++v2_neighbors_it;
  }

  std::size_t v3_count(v3Vertices_end - v3Vertices_begin);
  std::size_t v4_count(0);
  for (std::pair<const_vertex_iterator, connectivity_matrix> *it1(v3Vertices_begin);
       it1 != v3Vertices_end; ++it1) { // for each v3 vertex computed above
    const_edge_iterator v3_neighbors_it(g.neighbors_begin(it1->first));
    const_edge_iterator v3_neighbors_end(g.neighbors_end(it1->first));
    unsigned long v4_local_count(
        0); // keep track of how many v4 vertices are only the result of the neighbors of this v3
    for (std::pair<const_vertex_iterator, connectivity_matrix> *it2(v3Vertices_begin);
         it2 != v3Vertices_end; ++it2) { // consider other v3 vertices as candidate v4 vertices
      while (v3_neighbors_it != v3_neighbors_end &&
             g.target_of(v3_neighbors_it) <
                 it2->first) { // the v3 neighbor is exclusively a v4 vertex
        if (g.target_of(v3_neighbors_it) != v1 &&
            g.target_of(v3_neighbors_it) != v2) { // if this exclusively v4 vertex is not v1 or v2
          ++temp_edge_types.insert_or_zero(g.edge_value(v3_neighbors_it));
          ++v4_local_count;
          it1->second(0, 3) = 0;
          it1->second(1, 3) = 0;
          it1->second(2, 3) = g.edge_value(v3_neighbors_it);
          ++counts.insert(std::make_pair(mapper.canonical_subgraph_address(it1->second), 0))
                .first->second;
        }
        ++v3_neighbors_it;
      }
      if (v3_neighbors_it == v3_neighbors_end ||
          g.target_of(v3_neighbors_it) > it2->first) { // there is no edge between the v3 vertex and
                                                       // the other v3 vertex serving as a v4 vertex
        if (it1->first < it2->first) { // to be a candidate vertex, the other v3 vertex must be
                                       // greater to avoid double counting
          ++gaps;
          it1->second(0, 3) = it2->second(0, 2);
          it1->second(1, 3) = it2->second(1, 2);
          it1->second(2, 3) = 0;
          ++counts.insert(std::make_pair(mapper.canonical_subgraph_address(it1->second), 0))
                .first->second;
        }
      } else { // there is an edge between the v3 vertex and the other v3 vertex serving as a v4
               // vertex
        if (it1->first < it2->first) { // to be a candidate vertex, the other v3 vertex must be
                                       // greater to avoid double counting
          ++temp_edge_types.insert_or_zero(g.edge_value(v3_neighbors_it));
          it1->second(0, 3) = it2->second(0, 2);
          it1->second(1, 3) = it2->second(1, 2);
          it1->second(2, 3) = g.edge_value(v3_neighbors_it);
          ++counts.insert(std::make_pair(mapper.canonical_subgraph_address(it1->second), 0))
                .first->second;
        }
        ++v3_neighbors_it;
      }
    }
    while (v3_neighbors_it !=
           v3_neighbors_end) { // we have to be sure to go through the rest of the neighbors of v3,
                               // and all of these are exclusively v4
      if (g.target_of(v3_neighbors_it) != v1 && g.target_of(v3_neighbors_it) != v2) {
        ++temp_edge_types.insert_or_zero(g.edge_value(v3_neighbors_it));
        ++v4_local_count;
        it1->second(0, 3) = 0;
        it1->second(1, 3) = 0;
        it1->second(2, 3) = g.edge_value(v3_neighbors_it);
        ++counts.insert(std::make_pair(mapper.canonical_subgraph_address(it1->second), 0))
              .first->second;
      }
      ++v3_neighbors_it;
    }
    v4_count += v4_local_count;
    gaps += 2 * v4_local_count;
    it1->second(0, 3) = 0;
    it1->second(1, 3) = 0;
    it1->second(2, 3) = 0;
    counts.insert(std::make_pair(mapper.canonical_subgraph_address(it1->second), 0))
        .first->second += g.vertex_count() - 2 - v3_count - v4_local_count;
  }

  // account for the least connected substructures
  edge_types.for_each([&](connectivity_address_type key, unsigned long edge_count) {
    connectivity(2, 3) = key;
    unsigned long count = edge_count;
    unsigned long const *temp_val = temp_edge_types.find(key);
    if (temp_val != nullptr) {
      count -= *temp_val;
      if (key == 0) {
        count -= !static_cast<bool>(connectivity(0, 1)) +
                 (2 + v3_count) * (g.vertex_count() - 2 - v3_count) - 3 * v4_count;
      }
    }
    counts.insert(std::make_pair(mapper.canonical_subgraph_address(connectivity), 0))
        .first->second += count;
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

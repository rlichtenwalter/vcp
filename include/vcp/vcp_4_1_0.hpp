// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#ifndef VCP_VCP_4_1_0_HPP
#define VCP_VCP_4_1_0_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vcp/detail/v3_buffer_bound.hpp>
#include <vcp/graph.hpp>

namespace vcp {

template <std::size_t n, std::size_t r, bool d>
  requires(n >= 2 && r >= 1)
class vcp;

/**
 * @brief Full specialization of vcp for 4-vertex subgraphs on a single-relation undirected graph.
 *
 * Reduces the 4-vertex enumeration to a two-phase algorithm: first collect
 * the v3 candidate set by merging the sorted neighbor lists of v1 and v2,
 * then iterate v3 candidates to discover v4 vertices. Exactly 40 element
 * classes exist for this configuration.
 *
 * This specialization is NOT safe for shared-instance concurrent calls
 * because `generate_vector` uses the `v3Vertices` scratch buffer held as a
 * class member. Construct one `vcp<4, 1, false>` per thread.
 */
template <> class vcp<4, 1, false> {
private:
  constexpr static const std::size_t num_elements = 40;

public:
  /**
   * @brief Construct the VCP calculator bound to the given graph.
   *
   * Precomputes the global count of unconnected vertex pairs, which is
   * used to fill the least-connected element slot in `generate_vector`.
   *
   * @param g Undirected graph to analyze; must outlive this object.
   */
  vcp(graph const &g);

  /**
   * @brief Return the number of VCP element classes for this specialization.
   *
   * Always 40 for (n=4, r=1, d=false).
   *
   * @return 40.
   */
  [[nodiscard]] constexpr static std::size_t element_count() noexcept;

  /**
   * @brief Compute the VCP feature vector for the pivot pair (v1, v2).
   *
   * Returns a fixed-size array of 40 counts. The sum of all counts equals
   * C(|V|-2, 2) (the number of unordered pairs from the non-pivot vertices).
   *
   * @param v1 Iterator to the first pivot vertex.
   * @param v2 Iterator to the second pivot vertex.
   * @return Array of 40 occurrence counts indexed by element address.
   */
  [[nodiscard]] std::array<unsigned long, num_elements> generate_vector(const_vertex_iterator v1,
                                                                        const_vertex_iterator v2);

private:
  enum connectivity_value : std::size_t {
    V1V2 = 1,
    V1V3 = 2,
    V1V4 = 4,
    V2V3 = 8,
    V2V4 = 16,
    V3V4 = 32
  };
  graph const &g;
  constexpr static const std::size_t num_structures = 64;
  [[nodiscard]] static std::size_t element_address(std::size_t subgraph_address) noexcept;
  unsigned long unconnected_pairs;
  // Scratch buffer for v3 candidates accumulated during `generate_vector`.
  // Allocated once per instance, reused across every call — which makes
  // this class NOT thread-safe for shared-instance concurrent calls.
  // Construct one `vcp<4, 1, false>` per thread if parallelism is needed;
  // the underlying graph is safe to share for read. See README "Thread
  // safety" for the supported pattern.
  std::unique_ptr<std::pair<const_vertex_iterator, unsigned char>[]> v3Vertices;
};

inline std::size_t vcp<4, 1, false>::element_address(std::size_t subgraph_address) noexcept {
  constexpr static const std::array<std::size_t, num_structures> map = {
      {0,  1,  2,  3,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 6,  7,  10, 11, 8,  9,
       12, 13, 14, 15, 16, 17, 16, 17, 18, 19, 20, 21, 22, 23, 22, 23, 24, 25, 26, 27, 28, 29,
       30, 31, 32, 33, 26, 27, 30, 31, 28, 29, 32, 33, 34, 35, 36, 37, 36, 37, 38, 39}};
  return map[subgraph_address];
}

constexpr std::size_t vcp<4, 1, false>::element_count() noexcept { return num_elements; }

inline vcp<4, 1, false>::vcp(graph const &graph)
    : g(graph), unconnected_pairs((graph.vertex_count() * (graph.vertex_count() - 1) / 2) -
                                  graph.edge_count()),
      v3Vertices(std::make_unique<std::pair<const_vertex_iterator, unsigned char>[]>(
          detail::top_two_neighborhood_size_undirected(graph))) {}

inline std::array<unsigned long, vcp<4, 1, false>::element_count()>
vcp<4, 1, false>::generate_vector(const_vertex_iterator v1, const_vertex_iterator v2) {
  std::array<unsigned long, element_count()> counts = {{0}};

  // v1v2 is a small bitmask carrying the V1V2 connectivity bit.
  // unsigned char is wide enough; matches the storage type of
  // v3Vertices_end->second (the destination for v1v2+V1V3 etc. below).
  unsigned char const v1v2(static_cast<unsigned char>(V1V2 * g.edge_exists(v1, v2)));

  unsigned long connections(0);
  unsigned long gaps(0);

  const_edge_iterator v1_neighbors_it(g.neighbors_begin(v1));
  const_edge_iterator v1_neighbors_end(g.neighbors_end(v1));
  const_edge_iterator v2_neighbors_it(g.neighbors_begin(v2));
  const_edge_iterator v2_neighbors_end(g.neighbors_end(v2));
  std::pair<const_vertex_iterator, unsigned char> *v3Vertices_begin(&v3Vertices[0]);
  std::pair<const_vertex_iterator, unsigned char> *v3Vertices_end(&v3Vertices[0]);
  while (v1_neighbors_it != v1_neighbors_end && v2_neighbors_it != v2_neighbors_end) {
    if (g.target_of(v1_neighbors_it) < g.target_of(v2_neighbors_it)) {
      if (g.target_of(v1_neighbors_it) != v2) {
        ++connections;
        ++gaps;
        v3Vertices_end->first = g.target_of(v1_neighbors_it);
        v3Vertices_end->second = static_cast<unsigned char>(v1v2 + V1V3);
        ++v3Vertices_end;
      }
      ++v1_neighbors_it;
    } else if (g.target_of(v1_neighbors_it) > g.target_of(v2_neighbors_it)) {
      if (g.target_of(v2_neighbors_it) != v1) {
        ++connections;
        ++gaps;
        v3Vertices_end->first = g.target_of(v2_neighbors_it);
        v3Vertices_end->second = static_cast<unsigned char>(v1v2 + V2V3);
        ++v3Vertices_end;
      }
      ++v2_neighbors_it;
    } else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not
             // need to check to exclude it
      connections += 2;
      v3Vertices_end->first = g.target_of(v1_neighbors_it);
      v3Vertices_end->second = static_cast<unsigned char>(v1v2 + V1V3 + V2V3);
      ++v3Vertices_end;
      ++v1_neighbors_it;
      ++v2_neighbors_it;
    }
  }
  while (v1_neighbors_it != v1_neighbors_end) {
    if (g.target_of(v1_neighbors_it) != v2) {
      ++connections;
      ++gaps;
      v3Vertices_end->first = g.target_of(v1_neighbors_it);
      v3Vertices_end->second = static_cast<unsigned char>(v1v2 + V1V3);
      ++v3Vertices_end;
    }
    ++v1_neighbors_it;
  }
  while (v2_neighbors_it != v2_neighbors_end) {
    if (g.target_of(v2_neighbors_it) != v1) {
      ++connections;
      ++gaps;
      v3Vertices_end->first = g.target_of(v2_neighbors_it);
      v3Vertices_end->second = static_cast<unsigned char>(v1v2 + V2V3);
      ++v3Vertices_end;
    }
    ++v2_neighbors_it;
  }

  unsigned long v3_count(static_cast<unsigned long>(v3Vertices_end - v3Vertices_begin));
  unsigned long v4_count(0);
  for (std::pair<const_vertex_iterator, unsigned char> *it1(v3Vertices_begin);
       it1 != v3Vertices_end; ++it1) { // for each v3 vertex computed above
    const_edge_iterator v3_neighbors_it(g.neighbors_begin(it1->first));
    const_edge_iterator v3_neighbors_end(g.neighbors_end(it1->first));
    unsigned long v4_local_count(
        0); // keep track of how many v4 vertices are only the result of the neighbors of this v3
    for (std::pair<const_vertex_iterator, unsigned char> *it2(v3Vertices_begin);
         it2 != v3Vertices_end; ++it2) { // consider other v3 vertices as candidate v4 vertices
      while (v3_neighbors_it != v3_neighbors_end &&
             g.target_of(v3_neighbors_it) <
                 it2->first) { // the v3 neighbor is exclusively a v4 vertex
        if (g.target_of(v3_neighbors_it) != v1 &&
            g.target_of(v3_neighbors_it) != v2) { // if this exclusively v4 vertex is not v1 or v2
          ++v4_local_count;
        }
        ++v3_neighbors_it;
      }
      if (v3_neighbors_it == v3_neighbors_end ||
          g.target_of(v3_neighbors_it) > it2->first) { // there is no edge between the v3 vertex and
                                                       // the other v3 vertex serving as a v4 vertex
        if (it1->first < it2->first) { // to be a candidate vertex, the other v3 vertex must be
                                       // greater to avoid double counting
          ++gaps;
          // it2->second carries V1V3/V2V3 bits; promoting to v4-role requires
          // re-encoding into V1V4/V2V4 slots. The equivalence
          // `(it2->second - v1v2) << 1` holds because V1V4/V1V3 == V2V4/V2V3 == 2
          // — the per-slot stride for (r=1, d=0) — so a single shift covers all
          // three cases (V1V3-only → V1V4, V2V3-only → V2V4, both → V1V4+V2V4).
          // Stride here is 2 (shift by 1) because the d=0 enum uses 1 bit per
          // pair; vcp_4_1_1.hpp shifts by 2 because d=1 uses 2 bits per pair.
          // If the connectivity_value enum layout ever changes, replace with
          // the explicit form `V1V4 * ((it2->second - v1v2) / V1V3) + ...`.
          // it2->second >= v1v2 by domain (V1V3/V2V3 carry the same V1V2 bit
          // as v1v2), so the cast to size_t is value-preserving.
          ++counts[element_address(static_cast<std::size_t>(it1->second) +
                                   (static_cast<std::size_t>(it2->second - v1v2) << 1))];
        }
      } else { // there is an edge between the v3 vertex and the other v3 vertex serving as a v4
               // vertex
        if (it1->first < it2->first) { // to be a candidate vertex, the other v3 vertex must be
                                       // greater to avoid double counting
          ++connections;
          ++counts[element_address(static_cast<std::size_t>(it1->second) +
                                   (static_cast<std::size_t>(it2->second - v1v2) << 1) + V3V4)];
        }
        ++v3_neighbors_it;
      }
    }
    while (v3_neighbors_it !=
           v3_neighbors_end) { // we have to be sure to go through the rest of the neighbors of v3,
                               // and all of these are exclusively v4
      if (g.target_of(v3_neighbors_it) != v1 && g.target_of(v3_neighbors_it) != v2) {
        ++v4_local_count;
      }
      ++v3_neighbors_it;
    }
    v4_count += v4_local_count;
    connections += v4_local_count;
    gaps += 2 * v4_local_count;
    counts[element_address(it1->second + V3V4)] += v4_local_count;
    counts[element_address(it1->second)] += g.vertex_count() - 2 - v3_count - v4_local_count;
  }

  // account for the least connected substructures
  counts[element_address(v1v2 + V3V4)] = g.edge_count() - (connections + static_cast<bool>(v1v2));

  // the final computation is somewhat complicated
  // we start with the number of total gaps in the network, subtract the number of directly observed
  // gaps, then subtract the number of gaps contributed by triangles that we do not directly observe
  // that result from vertices that we do not directly observe the number of unobserved triangles is
  // recorded above as g.vertex_count() - 2 - v3_count - v4_local_count outside the loop this is
  // expressed as v3_count * (g.vertex_count() - 2 - v3_count) - v4_count this is also the number of
  // unobserved gaps that result from the v3 vertices in combination with vertices we do not
  // directly observe but we also have to account for the unobserved gaps between v1 and v2 and
  // these unobserved vertices this is pretty clearly 2*(g.vertex_count() - 2 - v3_count -
  // v4_count), a gap for each pairing of unobserved vertex with v1 or v2 the computation below is
  // thus equivalent to: unconnected_pairs - (gaps + !static_cast<bool>(v1v2)) - 2*(g.vertex_count()
  // - 2 - v3_count - v4_count) - v3_count*(g.vertex_count() - 2 - v3_count) + v4_count we can
  // simplify this expression as below
  // Debug-only guards for the K4-bidirectional-class underflow bug (see
  // CHANGELOG entry for the K4 fix on the prior branch). Every operand is
  // unsigned; a miscounted v3_count or v4_count would silently wrap one
  // of the sub-expressions and poison this slot with a ~2**64 value.
  // Guard each dangerous sub-expression — including `vertex_count - 2 - v3_count`
  // which can itself underflow if v3_count > V-2 — then guard the combined
  // RHS. Release build compiles these away under NDEBUG.
#ifndef NDEBUG
  {
    assert(g.vertex_count() >= 2 + v3_count &&
           "v3_count exceeds vertex_count - 2 (upstream miscounting)");
    auto const V = static_cast<long long>(g.vertex_count());
    long long const positive =
        static_cast<long long>(unconnected_pairs) + static_cast<long long>(3 * v4_count);
    long long const negative =
        static_cast<long long>(gaps + !static_cast<bool>(v1v2)) +
        (static_cast<long long>(2 + v3_count) * (V - 2 - static_cast<long long>(v3_count)));
    assert(positive >= negative && "unconnected_pairs subtraction chain underflowed");
    assert(positive - negative <= static_cast<long long>(unconnected_pairs) &&
           "unconnected_pairs subtraction chain overshot its upper bound");
  }
#endif
  counts[element_address(v1v2)] = unconnected_pairs - (gaps + !static_cast<bool>(v1v2)) -
                                  (2 + v3_count) * (g.vertex_count() - 2 - v3_count) + 3 * v4_count;

  return counts;
}

} // namespace vcp

#endif

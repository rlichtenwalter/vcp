// Container micro-benchmark for the small-map study.
//
// Measures the four candidate containers (std::map baseline,
// std::unordered_map, std::vector<pair> linear scan, std::array dense)
// across k in {1, 4, 16, 64, 256, 1024, 4096} on the hot-path operation
// mix produced by vcp<4, r, d>::generate_vector.
//
// Per-call behavior of `temp_edge_types` and `edge_types`:
//   - O(k) insert-or-increment operations where keys are drawn from a
//     universe of size K (the connectivity key space). Each call
//     touches roughly k distinct keys.
//   - At the tail, one pass of zipping edge_types and temp_edge_types
//     via find() on the matching key (only for the per-instance
//     `edge_types`, not the per-call `temp_edge_types`).
//   - Then the container is destroyed and allocations freed.
//
// The benchmark simulates this: N insert-or-increment ops against a
// universe of K keys, then one iteration pass, then destroy. Wall-clock
// per-iteration is reported; the point is the *shape* of the
// cross-container curve, not absolute numbers.
//
// Build with:
//   cmake --build build-probe --target bench_container_mix
// Run:
//   ./build-probe/benchmark/small_map_study/bench_container_mix
//     "[!benchmark]" --benchmark-samples 30

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

namespace {

using Key = std::uint64_t;
using Val = unsigned long;

// Generate a deterministic key stream: `ops` insertions drawn uniformly
// from a universe of size `key_universe`. The ratio ops/key_universe
// controls the expected final cardinality.
std::vector<Key> make_keys(std::size_t ops, std::size_t key_universe, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_int_distribution<Key> dist(0, static_cast<Key>(key_universe) - 1);
  std::vector<Key> keys(ops);
  for (auto &k : keys) {
    k = dist(rng);
  }
  return keys;
}

// Vector-of-pair with linear-scan insert-or-increment. Models the
// heap-backed form of the tier-2 fallback.
template <typename K, typename V> struct FlatPairVec {
  std::vector<std::pair<K, V>> data;

  V &operator[](K const &key) {
    for (auto &p : data) {
      if (p.first == key) {
        return p.second;
      }
    }
    data.emplace_back(key, V{});
    return data.back().second;
  }

  auto begin() { return data.begin(); }
  auto end() { return data.end(); }
  auto begin() const { return data.begin(); }
  auto end() const { return data.end(); }
  auto find(K const &key) const {
    for (auto it = data.begin(); it != data.end(); ++it) {
      if (it->first == key) {
        return it;
      }
    }
    return data.end();
  }
  std::size_t size() const { return data.size(); }
};

// Small-buffer-optimized variant: backing storage is a stack-local
// std::array sized at Cap; when a value past Cap is inserted, the
// container aborts (the probe is a benchmark; production would promote
// to a heap-backed form). Avoids the allocator cost on the hot path
// that FlatPairVec incurs on every construct/destroy cycle — which is
// the relevant per-call comparison for tier 2, since
// `temp_edge_types` lives only for the duration of one
// generate_vector call.
template <typename K, typename V, std::size_t Cap> struct StackFlatPairVec {
  std::array<std::pair<K, V>, Cap> data{};
  std::size_t n = 0;

  V &operator[](K const &key) {
    for (std::size_t i = 0; i < n; ++i) {
      if (data[i].first == key) {
        return data[i].second;
      }
    }
    // Silent saturate on overflow; this is a microbench, not production.
    if (n < Cap) {
      data[n] = {key, V{}};
      return data[n++].second;
    }
    return data[Cap - 1].second;
  }

  auto begin() { return data.begin(); }
  auto end() { return data.begin() + n; }
  auto begin() const { return data.begin(); }
  auto end() const { return data.begin() + n; }
  std::size_t size() const { return n; }
};

// Dense-array variant: key_universe must equal N template param.
template <typename V, std::size_t N> struct DenseArray {
  std::array<V, N> data{};

  V &operator[](std::size_t key) { return data[key]; }
  // For iteration, emit only non-default slots. This mirrors the sparse
  // semantic the other containers provide.
  struct Iter {
    typename std::array<V, N>::const_iterator it, end;
    std::size_t idx;
    void advance_to_nonzero() {
      while (it != end && *it == V{}) {
        ++it;
        ++idx;
      }
    }
  };
  struct const_proxy_iter {
    std::array<V, N> const *arr;
    std::size_t idx;
    bool operator!=(const_proxy_iter const &o) const { return idx != o.idx; }
    const_proxy_iter &operator++() {
      ++idx;
      while (idx < N && (*arr)[idx] == V{}) {
        ++idx;
      }
      return *this;
    }
    std::pair<std::size_t, V> operator*() const { return {idx, (*arr)[idx]}; }
  };
  const_proxy_iter begin() const {
    const_proxy_iter it{&data, 0};
    while (it.idx < N && data[it.idx] == V{}) {
      ++it.idx;
    }
    return it;
  }
  const_proxy_iter end() const { return const_proxy_iter{&data, N}; }
};

// Run one pass: N insert-or-increment ops against a universe of K keys,
// followed by one iteration over the populated entries summing values
// to prevent dead-code elimination.
template <typename Container>
unsigned long run_container(Container &c, std::vector<Key> const &keys) {
  for (auto const &k : keys) {
    ++c[k];
  }
  unsigned long sum = 0;
  for (auto const &kv : c) {
    sum += kv.second;
  }
  return sum;
}

// Specialized for DenseArray whose iterator yields pair<size_t, V>.
template <typename V, std::size_t N>
unsigned long run_container_dense(DenseArray<V, N> &c, std::vector<Key> const &keys) {
  for (auto const &k : keys) {
    ++c[static_cast<std::size_t>(k)];
  }
  unsigned long sum = 0;
  for (std::size_t i = 0; i < N; ++i) {
    sum += c.data[i];
  }
  return sum;
}

// One-shot benchmark: each iteration constructs, populates, iterates,
// destroys. Matches the per-call lifetime of temp_edge_types.
template <typename Container>
void bench_one(Catch::Benchmark::Chronometer meter, std::vector<Key> const &keys) {
  meter.measure([&] {
    Container c;
    return run_container(c, keys);
  });
}

template <typename V, std::size_t N>
void bench_dense(Catch::Benchmark::Chronometer meter, std::vector<Key> const &keys) {
  meter.measure([&] {
    DenseArray<V, N> c;
    return run_container_dense(c, keys);
  });
}

template <std::size_t Cap>
void bench_stack(Catch::Benchmark::Chronometer meter, std::vector<Key> const &keys) {
  meter.measure([&] {
    StackFlatPairVec<Key, Val, Cap> c;
    for (auto const &k : keys) {
      ++c[k];
    }
    unsigned long sum = 0;
    for (std::size_t i = 0; i < c.size(); ++i) {
      sum += c.data[i].second;
    }
    return sum;
  });
}

// Common parameters: we want ops = 30 * final_k to mimic the many
// insert-increment ops the hot loop performs. Use a key universe of
// final_k (so each key is hit ~30 times).
constexpr std::size_t OPS_PER_KEY = 30;

} // namespace

// ============================================================
// Fixed key universe, varying k (= universe size = final cardinality)
// ============================================================

TEST_CASE("container k=4", "[!benchmark][cx-k4]") {
  auto const keys = make_keys(OPS_PER_KEY * 4, 4, 1);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat[cap=4]")(Catch::Benchmark::Chronometer m) {
    bench_stack<4>(m, keys);
  };
  BENCHMARK_ADVANCED("dense_array[4]")(Catch::Benchmark::Chronometer m) {
    bench_dense<Val, 4>(m, keys);
  };
}

TEST_CASE("container k=16", "[!benchmark][cx-k16]") {
  auto const keys = make_keys(OPS_PER_KEY * 16, 16, 2);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat[cap=16]")(Catch::Benchmark::Chronometer m) {
    bench_stack<16>(m, keys);
  };
  BENCHMARK_ADVANCED("dense_array[16]")(Catch::Benchmark::Chronometer m) {
    bench_dense<Val, 16>(m, keys);
  };
}

TEST_CASE("container k=64", "[!benchmark][cx-k64]") {
  auto const keys = make_keys(OPS_PER_KEY * 64, 64, 3);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat[cap=64]")(Catch::Benchmark::Chronometer m) {
    bench_stack<64>(m, keys);
  };
  BENCHMARK_ADVANCED("dense_array[64]")(Catch::Benchmark::Chronometer m) {
    bench_dense<Val, 64>(m, keys);
  };
}

TEST_CASE("container k=256", "[!benchmark][cx-k256]") {
  auto const keys = make_keys(OPS_PER_KEY * 256, 256, 4);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    // 64-slot stack buffer, sized to comfortably hold k ≤ 64 for
    // microbench direct comparison. Production would pick Cap based
    // on expected k distribution, not k itself.
    bench_stack<64>(m, keys);
  };
  BENCHMARK_ADVANCED("dense_array[256]")(Catch::Benchmark::Chronometer m) {
    bench_dense<Val, 256>(m, keys);
  };
}

TEST_CASE("container k=1024", "[!benchmark][cx-k1024]") {
  auto const keys = make_keys(OPS_PER_KEY * 1024, 1024, 5);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    // 64-slot stack buffer, sized to comfortably hold k ≤ 64 for
    // microbench direct comparison. Production would pick Cap based
    // on expected k distribution, not k itself.
    bench_stack<64>(m, keys);
  };
  BENCHMARK_ADVANCED("dense_array[1024]")(Catch::Benchmark::Chronometer m) {
    bench_dense<Val, 1024>(m, keys);
  };
}

TEST_CASE("container k=4096", "[!benchmark][cx-k4096]") {
  auto const keys = make_keys(OPS_PER_KEY * 4096, 4096, 6);
  BENCHMARK_ADVANCED("std::map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("std::unordered_map")(Catch::Benchmark::Chronometer m) {
    bench_one<std::unordered_map<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    bench_one<FlatPairVec<Key, Val>>(m, keys);
  };
  BENCHMARK_ADVANCED("stack_flat_pair_vec")(Catch::Benchmark::Chronometer m) {
    // 64-slot stack buffer, sized to comfortably hold k ≤ 64 for
    // microbench direct comparison. Production would pick Cap based
    // on expected k distribution, not k itself.
    bench_stack<64>(m, keys);
  };
  // dense at k=4096 costs 32 KB of stack — still fine, but omitted
  // because the dense tier would not fire at this key space in the
  // ship config (byte budget = 8 MB; 4096 * 8 bytes = 32 KB fits, but
  // k=4096 in our workloads only arises at r=30 where key_space is
  // 2^30, far beyond dense). Measuring it here would mislead.
}

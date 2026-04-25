// Container micro-benchmark for the `counts` return-value of
// vcp<4, r, d>::generate_vector.
//
// `counts` is semantically different from the hot-loop scratch maps
// (temp_edge_types, edge_types). It is populated once per call via
// insert-or-increment and ++ operations and then iterated once in
// key-sorted order at the caller boundary (today: std::map's ordered
// iteration; tomorrow: a materialized sorted vector<pair>).
//
// The candidates below model five replacement designs plus the std::map
// baseline:
//
//   StdMap           - current implementation: std::map<Key, Val>
//   HashSortedKeys   - unordered_map + sorted vector<Key> cache;
//                      iteration looks up each key in the hash
//   HashSortedIters  - unordered_map + sorted vector<iterator> cache;
//                      iteration dereferences iterators (valid only if
//                      rehash never occurred since cache built)
//   HashSortedPairs  - unordered_map + materialized sorted
//                      vector<pair<Key,Val>> cache
//   SortedVector     - sorted vector<pair<Key,Val>>; binary-search
//                      insert shifts tail on miss
//   BatchedDedup     - append-only vector<pair<Key,Val>>; final sort
//                      + consecutive-key reduction at finalize
//
// Each benchmark iteration: construct, run all increment operations,
// then finalize to a sorted vector<pair<Key,Val>>, then destroy. The
// sum over values is accumulated to prevent DCE. Measurements are total
// wall-clock per iteration.
//
// Scenarios span the (N, K, distribution) space that real generate_vector
// workloads produce:
//
//   tiny-hot            K=16     N=100k   uniform      (small-graph hot)
//   small-zipf          K=256    N=1M     Zipf(s=1.2)  (typical run)
//   medium-uniform      K=4096   N=500k   uniform      (transition)
//   large-sparse        K=100k   N=200k   uniform      (large-r)
//   huge-oneshot        K=1M     N=1M     uniform      (pure dedup)
//   adversarial-desc    K=5k     N=5k     descending   (SortedVec worst)
//
// Build with:
//   cmake -S . -B build-counts -DVCP_BUILD_BENCHMARKS=ON
//     -DVCP_BUILD_COUNTS_BENCH=ON -DCMAKE_BUILD_TYPE=Release
//   cmake --build build-counts --target bench_counts_containers
// Run:
//   ./build-counts/benchmark/counts_map_study/bench_counts_containers
//     "[!benchmark]" --benchmark-samples 15

#include <algorithm>
#include <cmath>
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
using Pair = std::pair<Key, Val>;

// ============================================================
// Key-stream generators.
// ============================================================

// Uniform random keys in [0, K).
std::vector<Key> make_uniform(std::size_t n, std::size_t k, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Key> dist(0, static_cast<Key>(k) - 1);
  std::vector<Key> out(n);
  for (auto &x : out) {
    x = dist(rng);
  }
  return out;
}

// Zipf-distributed keys via inverse CDF sampling. Approximate but
// deterministic: we build the CDF of rank probabilities proportional to
// 1/rank^s, then sample.
std::vector<Key> make_zipf(std::size_t n, std::size_t k, double s, unsigned seed) {
  std::vector<double> cdf(k);
  double accum = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    accum += 1.0 / std::pow(static_cast<double>(i + 1), s);
    cdf[i] = accum;
  }
  double const total = accum;
  for (auto &c : cdf) {
    c /= total;
  }
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<Key> out(n);
  for (auto &x : out) {
    double const r = u(rng);
    auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
    x = static_cast<Key>(it - cdf.begin());
  }
  return out;
}

// Descending keys: N distinct, inserted in strictly decreasing order
// (K-1, K-2, ..., 0). Worst case for SortedVector (every insert shifts
// the entire tail).
std::vector<Key> make_descending(std::size_t k) {
  std::vector<Key> out(k);
  for (std::size_t i = 0; i < k; ++i) {
    out[i] = static_cast<Key>(k - 1 - i);
  }
  return out;
}

// ============================================================
// Candidates.
// ============================================================

struct StdMap {
  std::map<Key, Val> m;
  void add(Key const &k, Val w) { m[k] += w; }
  std::vector<Pair> finalize_sorted() const { return {m.begin(), m.end()}; }
};

// Hash + sorted vector<Key>; iteration hashes each key for its value.
struct HashSortedKeys {
  std::unordered_map<Key, Val> h;
  void add(Key const &k, Val w) { h[k] += w; }
  std::vector<Pair> finalize_sorted() {
    std::vector<Key> keys;
    keys.reserve(h.size());
    for (auto const &kv : h) {
      keys.push_back(kv.first);
    }
    std::sort(keys.begin(), keys.end());
    std::vector<Pair> out;
    out.reserve(keys.size());
    for (auto const &k : keys) {
      out.emplace_back(k, h.at(k));
    }
    return out;
  }
};

// Hash + sorted vector<iterator>; iteration dereferences.
struct HashSortedIters {
  std::unordered_map<Key, Val> h;
  void add(Key const &k, Val w) { h[k] += w; }
  std::vector<Pair> finalize_sorted() {
    using It = std::unordered_map<Key, Val>::iterator;
    std::vector<It> iters;
    iters.reserve(h.size());
    for (auto it = h.begin(); it != h.end(); ++it) {
      iters.push_back(it);
    }
    std::sort(iters.begin(), iters.end(),
              [](It const &a, It const &b) { return a->first < b->first; });
    std::vector<Pair> out;
    out.reserve(iters.size());
    for (auto const &it : iters) {
      out.emplace_back(it->first, it->second);
    }
    return out;
  }
};

// Hash + materialized sorted vector<pair>; built at finalize.
struct HashSortedPairs {
  std::unordered_map<Key, Val> h;
  void add(Key const &k, Val w) { h[k] += w; }
  std::vector<Pair> finalize_sorted() {
    std::vector<Pair> out(h.begin(), h.end());
    std::sort(out.begin(), out.end(),
              [](Pair const &a, Pair const &b) { return a.first < b.first; });
    return out;
  }
};

// Sorted vector kept sorted at every insert via lower_bound + shift.
struct SortedVector {
  std::vector<Pair> v;
  void add(Key const &k, Val w) {
    auto it = std::lower_bound(v.begin(), v.end(), k,
                               [](Pair const &p, Key const &x) { return p.first < x; });
    if (it != v.end() && it->first == k) {
      it->second += w;
      return;
    }
    v.insert(it, Pair{k, w});
  }
  std::vector<Pair> finalize_sorted() { return v; }
};

// Append-only; final sort + scan-reduce at finalize.
struct BatchedDedup {
  std::vector<Pair> records;
  void add(Key const &k, Val w) { records.emplace_back(k, w); }
  std::vector<Pair> finalize_sorted() {
    if (records.empty()) {
      return {};
    }
    std::sort(records.begin(), records.end(),
              [](Pair const &a, Pair const &b) { return a.first < b.first; });
    std::vector<Pair> out;
    out.reserve(records.size());
    Key cur = records[0].first;
    Val sum = records[0].second;
    for (std::size_t i = 1; i < records.size(); ++i) {
      if (records[i].first == cur) {
        sum += records[i].second;
      } else {
        out.emplace_back(cur, sum);
        cur = records[i].first;
        sum = records[i].second;
      }
    }
    out.emplace_back(cur, sum);
    return out;
  }
};

// ============================================================
// Benchmark kernel. One iteration = construct, run all increments,
// finalize, consume (sum), destroy. Matches the lifetime of a single
// generate_vector call's `counts` map.
// ============================================================

template <typename Container> Val bench_kernel(std::vector<Key> const &keys) {
  Container c;
  for (auto const &k : keys) {
    c.add(k, 1);
  }
  auto sorted = c.finalize_sorted();
  Val sum = 0;
  for (auto const &p : sorted) {
    sum += p.second;
  }
  return sum;
}

template <typename Container>
void bench_one(Catch::Benchmark::Chronometer meter, std::vector<Key> const &keys) {
  meter.measure([&] { return bench_kernel<Container>(keys); });
}

// Macro for stamping out the six candidates per scenario.
#define RUN_ALL_CANDIDATES(keys_var)                                                               \
  BENCHMARK_ADVANCED("StdMap")(Catch::Benchmark::Chronometer m) {                                  \
    bench_one<StdMap>(m, (keys_var));                                                              \
  };                                                                                               \
  BENCHMARK_ADVANCED("HashSortedKeys")(Catch::Benchmark::Chronometer m) {                          \
    bench_one<HashSortedKeys>(m, (keys_var));                                                      \
  };                                                                                               \
  BENCHMARK_ADVANCED("HashSortedIters")(Catch::Benchmark::Chronometer m) {                         \
    bench_one<HashSortedIters>(m, (keys_var));                                                     \
  };                                                                                               \
  BENCHMARK_ADVANCED("HashSortedPairs")(Catch::Benchmark::Chronometer m) {                         \
    bench_one<HashSortedPairs>(m, (keys_var));                                                     \
  };                                                                                               \
  BENCHMARK_ADVANCED("SortedVector")(Catch::Benchmark::Chronometer m) {                            \
    bench_one<SortedVector>(m, (keys_var));                                                        \
  };                                                                                               \
  BENCHMARK_ADVANCED("BatchedDedup")(Catch::Benchmark::Chronometer m) {                            \
    bench_one<BatchedDedup>(m, (keys_var));                                                        \
  }

} // namespace

// ============================================================
// Scenarios.
// ============================================================

TEST_CASE("counts tiny-hot (K=16, N=100k, uniform)", "[!benchmark][cn-tiny-hot]") {
  auto const keys = make_uniform(100'000, 16, 1);
  RUN_ALL_CANDIDATES(keys);
}

TEST_CASE("counts small-zipf (K=256, N=1M, Zipf s=1.2)", "[!benchmark][cn-small-zipf]") {
  auto const keys = make_zipf(1'000'000, 256, 1.2, 2);
  RUN_ALL_CANDIDATES(keys);
}

TEST_CASE("counts medium-uniform (K=4096, N=500k, uniform)", "[!benchmark][cn-medium-uniform]") {
  auto const keys = make_uniform(500'000, 4096, 3);
  RUN_ALL_CANDIDATES(keys);
}

TEST_CASE("counts large-sparse (K=100k, N=200k, uniform)", "[!benchmark][cn-large-sparse]") {
  auto const keys = make_uniform(200'000, 100'000, 4);
  RUN_ALL_CANDIDATES(keys);
}

TEST_CASE("counts huge-oneshot (K=1M, N=1M, uniform)", "[!benchmark][cn-huge-oneshot]") {
  auto const keys = make_uniform(1'000'000, 1'000'000, 5);
  RUN_ALL_CANDIDATES(keys);
}

TEST_CASE("counts adversarial-desc (K=5k, N=5k, descending)", "[!benchmark][cn-adversarial-desc]") {
  auto const keys = make_descending(5'000);
  RUN_ALL_CANDIDATES(keys);
}

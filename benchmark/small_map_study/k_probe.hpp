// Instrumentation hooks for the std::map → detail::dense_or_sparse_map study.
// Private to benchmark/small_map_study — NOT part of the public library.
//
// When `VCP_INSTRUMENT_K` is defined at compile time (set by the CMake
// option `VCP_BUILD_K_PROBE`), vcp_4_r_0.hpp and vcp_4_r_1.hpp call into
// `k_probe::record()` at the tail of `generate_vector` to log the final
// cardinality of each hot map. Linear histograms get cramped for the
// long-tailed distributions we expect (most calls have k in single
// digits; a few on hub vertices may reach 10^3+), so we use log2 buckets.
//
// No-op when `VCP_INSTRUMENT_K` is undefined: the headers simply do not
// `#include` this file, and the probe calls disappear.
//
// Output: a CSV written to the path in `VCP_K_PROBE_OUTPUT` (env var) at
// atexit, or `/tmp/k_probe.csv` by default. Columns:
//   map_name, workload, bucket_log2, bucket_count, bucket_lo, bucket_hi
//
// This file is deliberately non-thread-safe. The probe is intended for
// single-threaded runs of the bench tool; the point is to measure k, not
// to validate the thread-safety story (see Finding 9 for that, deferred
// to a separate PR).

#ifndef VCP_SMALL_MAP_STUDY_K_PROBE
#define VCP_SMALL_MAP_STUDY_K_PROBE

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>

namespace k_probe {

// 64 buckets covers k in [0, 2^63]. In practice the interesting range is
// [0, 2^14] for our workloads (MAX_NEIGHBORS is 16384).
inline constexpr std::size_t BUCKETS = 64;

struct Histogram {
  std::array<std::size_t, BUCKETS> counts{};
  std::size_t total_samples = 0;
  std::size_t max_observed = 0;
};

inline std::unordered_map<std::string, Histogram> &histograms() {
  static std::unordered_map<std::string, Histogram> h;
  return h;
}

inline std::string &current_workload() {
  static std::string w;
  return w;
}

inline void set_workload(std::string const &w) { current_workload() = w; }

// log2 bucket: for k = 0 bucket = 0, for k in [2^(b-1), 2^b - 1] bucket = b.
// Not std::bit_width to keep the hot path dependency-free (no <bit>).
inline std::size_t bucket_of(std::size_t k) {
  if (k == 0) {
    return 0;
  }
  std::size_t b = 1;
  while ((std::size_t{1} << b) <= k) {
    ++b;
  }
  return b;
}

inline void flush_csv(char const *path);

inline void ensure_flush_registered() {
  // Lazy registration: first call triggers the once-only atexit hook.
  // Two footguns dealt with here:
  //   (1) Inline-variable dynamic init may be skipped if the variable is
  //       not odr-used from anywhere — the first version of this file had
  //       a top-level `inline int const _flush_registered = ...;` that
  //       never fired because the compiler saw nothing touching it. The
  //       static-local-inside-a-function pattern below is odr-used by
  //       every record() call and cannot be elided.
  //   (2) Static destruction order. If `atexit` is registered BEFORE
  //       `histograms()` is first touched, then at program exit the
  //       histograms destructor runs first (reverse of construction
  //       order) and the atexit handler reads a destroyed container.
  //       Symptom: an empty CSV despite record() firing thousands of
  //       times. The force-touch `(void)histograms();` below constructs
  //       the map BEFORE atexit, so the atexit handler always runs
  //       first and observes a live map.
  static int const _once = []() {
    (void)histograms();
    std::atexit([]() {
      char const *p = std::getenv("VCP_K_PROBE_OUTPUT");
      flush_csv(p ? p : "/tmp/k_probe.csv");
    });
    return 0;
  }();
  (void)_once;
}

inline void record(char const *map_name, std::size_t k) {
  ensure_flush_registered();
  std::string const key = std::string(map_name) + "|" + current_workload();
  Histogram &h = histograms()[key];
  ++h.total_samples;
  ++h.counts[bucket_of(k)];
  if (k > h.max_observed) {
    h.max_observed = k;
  }
}

inline void flush_csv(char const *path) {
  std::FILE *f = std::fopen(path, "w");
  if (!f) {
    std::fprintf(stderr, "k_probe: cannot open %s for writing\n", path);
    return;
  }
  std::fprintf(f, "map_name,workload,bucket_log2,bucket_count,bucket_lo,bucket_hi,"
                  "total_samples,max_observed\n");
  for (auto const &[key, h] : histograms()) {
    auto const pipe = key.find('|');
    std::string const map_name = key.substr(0, pipe);
    std::string const workload = key.substr(pipe + 1);
    for (std::size_t b = 0; b < BUCKETS; ++b) {
      if (h.counts[b] == 0) {
        continue;
      }
      std::size_t const lo = (b == 0) ? 0 : (std::size_t{1} << (b - 1));
      std::size_t const hi = (b == 0) ? 0 : ((std::size_t{1} << b) - 1);
      std::fprintf(f, "%s,%s,%zu,%zu,%zu,%zu,%zu,%zu\n", map_name.c_str(), workload.c_str(), b,
                   h.counts[b], lo, hi, h.total_samples, h.max_observed);
    }
  }
  std::fclose(f);
}

} // namespace k_probe

#endif

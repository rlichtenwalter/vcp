// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2013-2026 Ryan N. Lichtenwalter

#include <cstddef>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <vcp/directed_graph.hpp>
#include <vcp/graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/vcp.hpp>

// Library-level micro-benchmarks.
//
// Complements the tool-level cross-ref benchmarks under benchmark/run.sh.
// This target only builds against the modernized CMake/Catch2 setup, so
// it cannot measure the 2012 legacy code. Its role is fine-grained work
// on the modern headers - tracking the effect of a single change, or a
// series of changes, within a single build of the library.
//
// All benchmarks generate deterministic Erdos-Renyi graphs in a fixed
// stringstream and feed them through the library's operator>> parser.
// This exercises both the I/O path and the algorithm, matching what the
// CLI tools do end-to-end.
//
// Usage (after building with -DVCP_BUILD_BENCHMARKS=ON):
//
//   ./build/benchmark/bench_vcp "[!benchmark]"
//   ./build/benchmark/bench_vcp "[!benchmark][io]"
//   ./build/benchmark/bench_vcp "[!benchmark][vcp-4-1-0]"
//   ./build/benchmark/bench_vcp "[!benchmark]" --benchmark-samples 50

namespace {

std::string erdos_renyi_undirected_text(std::size_t n, double p, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<std::vector<std::size_t>> adj(n);
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t v = u + 1; v < n; ++v) {
      if (dist(rng) < p) {
        adj[u].push_back(v);
        adj[v].push_back(u);
      }
    }
  }
  std::ostringstream os;
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t i = 0; i < adj[u].size(); ++i) {
      if (i > 0) {
        os << ' ';
      }
      os << adj[u][i];
    }
    os << '\n';
  }
  return os.str();
}

std::string erdos_renyi_bidirectional_text(std::size_t n, double p, unsigned seed) {
  // Every sampled pair becomes two arcs so the graph is symmetric. Matches
  // the fixture generation in generate_fixtures.py (dsym_*) for consistent
  // comparisons across the two benchmark layers.
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<std::vector<std::size_t>> out_adj(n);
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t v = u + 1; v < n; ++v) {
      if (dist(rng) < p) {
        out_adj[u].push_back(v);
        out_adj[v].push_back(u);
      }
    }
  }
  std::ostringstream os;
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t i = 0; i < out_adj[u].size(); ++i) {
      if (i > 0) {
        os << ' ';
      }
      os << out_adj[u][i];
    }
    os << '\n';
  }
  return os.str();
}

std::string erdos_renyi_multirelational_text(std::size_t n, double p, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> fdist(0.0, 1.0);
  std::uniform_int_distribution<int> bits(1, 3);
  std::vector<std::vector<std::pair<std::size_t, int>>> adj(n);
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t v = u + 1; v < n; ++v) {
      if (fdist(rng) < p) {
        int const b = bits(rng);
        adj[u].emplace_back(v, b);
        adj[v].emplace_back(u, b);
      }
    }
  }
  std::ostringstream os;
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t i = 0; i < adj[u].size(); ++i) {
      if (i > 0) {
        os << ' ';
      }
      os << adj[u][i].first << ',' << adj[u][i].second;
    }
    os << '\n';
  }
  return os.str();
}

template <typename Graph> Graph parse_graph(std::string const &text) {
  std::istringstream is(text);
  Graph g;
  is >> g;
  return g;
}

template <std::size_t n_, std::size_t r_, bool d_, typename Graph>
unsigned long vcp_all_pairs_sum(Graph const &g) {
  // Sum the first component of each VCP vector to defeat dead-code
  // elimination without printing to stdout (which would dominate the
  // measurement).
  vcp::vcp<n_, r_, d_> profiler(g);
  unsigned long accumulator = 0;
  auto vBegin = g.vertices_begin();
  std::size_t const vertex_count = g.vertex_count();
  for (std::size_t i = 0; i < vertex_count; ++i) {
    for (std::size_t j = i + 1; j < vertex_count; ++j) {
      auto vec = profiler.generate_vector(vBegin + i, vBegin + j);
      accumulator += vec[0];
    }
  }
  return accumulator;
}

} // namespace

// ============================================================
// Graph I/O benchmarks
// ============================================================

TEST_CASE("bench: graph parse", "[!benchmark][io]") {
  auto const text_n500 = erdos_renyi_undirected_text(500, 0.05, 42);
  auto const text_n1000 = erdos_renyi_undirected_text(1000, 0.02, 42);
  auto const text_n2000 = erdos_renyi_undirected_text(2000, 0.01, 42);

  BENCHMARK_ADVANCED("undirected n=500 p=0.05")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return parse_graph<vcp::graph>(text_n500); });
  };

  BENCHMARK_ADVANCED("undirected n=1000 p=0.02")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return parse_graph<vcp::graph>(text_n1000); });
  };

  BENCHMARK_ADVANCED("undirected n=2000 p=0.01")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return parse_graph<vcp::graph>(text_n2000); });
  };
}

TEST_CASE("bench: graph serialize", "[!benchmark][io]") {
  auto const text_n1000 = erdos_renyi_undirected_text(1000, 0.02, 42);
  auto const g = parse_graph<vcp::graph>(text_n1000);

  BENCHMARK_ADVANCED("undirected n=1000 operator<<")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      std::ostringstream os;
      os << g;
      return os.str().size();
    });
  };
}

// ============================================================
// Undirected VCP enumeration (the primary workload)
// ============================================================

TEST_CASE("bench: vcp<3,1,0> all pairs", "[!benchmark][vcp-3-1-0]") {
  auto const g200 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(200, 0.10, 42));
  auto const g500 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(500, 0.05, 42));
  auto const g1000 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(1000, 0.02, 42));

  BENCHMARK_ADVANCED("n=200 p=0.10")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 1, false>)(g200); });
  };

  BENCHMARK_ADVANCED("n=500 p=0.05")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 1, false>)(g500); });
  };

  BENCHMARK_ADVANCED("n=1000 p=0.02")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 1, false>)(g1000); });
  };
}

TEST_CASE("bench: vcp<4,1,0> all pairs", "[!benchmark][vcp-4-1-0]") {
  auto const g200 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(200, 0.10, 42));
  auto const g500 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(500, 0.05, 42));
  auto const g1000 = parse_graph<vcp::graph>(erdos_renyi_undirected_text(1000, 0.02, 42));

  BENCHMARK_ADVANCED("n=200 p=0.10")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 1, false>)(g200); });
  };

  BENCHMARK_ADVANCED("n=500 p=0.05")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 1, false>)(g500); });
  };

  BENCHMARK_ADVANCED("n=1000 p=0.02")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 1, false>)(g1000); });
  };
}

// ============================================================
// Directed VCP enumeration
// ============================================================

TEST_CASE("bench: vcp<3,1,1> all pairs", "[!benchmark][vcp-3-1-1]") {
  auto const g200 = parse_graph<vcp::directed_graph>(erdos_renyi_bidirectional_text(200, 0.10, 42));
  auto const g500 = parse_graph<vcp::directed_graph>(erdos_renyi_bidirectional_text(500, 0.05, 42));

  BENCHMARK_ADVANCED("bidirectional n=200 p=0.10")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 1, true>)(g200); });
  };

  BENCHMARK_ADVANCED("bidirectional n=500 p=0.05")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 1, true>)(g500); });
  };
}

TEST_CASE("bench: vcp<4,1,1> all pairs", "[!benchmark][vcp-4-1-1]") {
  auto const g200 = parse_graph<vcp::directed_graph>(erdos_renyi_bidirectional_text(200, 0.10, 42));
  auto const g500 = parse_graph<vcp::directed_graph>(erdos_renyi_bidirectional_text(500, 0.05, 42));

  BENCHMARK_ADVANCED("bidirectional n=200 p=0.10")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 1, true>)(g200); });
  };

  BENCHMARK_ADVANCED("bidirectional n=500 p=0.05")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 1, true>)(g500); });
  };
}

// ============================================================
// Multirelational VCP (r=2)
// ============================================================

TEST_CASE("bench: vcp<3,2,0> all pairs", "[!benchmark][vcp-3-2-0]") {
  auto const g100 =
      parse_graph<vcp::multirelational_graph<2>>(erdos_renyi_multirelational_text(100, 0.30, 42));
  auto const g200 =
      parse_graph<vcp::multirelational_graph<2>>(erdos_renyi_multirelational_text(200, 0.15, 42));

  BENCHMARK_ADVANCED("n=100 p=0.30 r=2")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 2, false>)(g100); });
  };

  BENCHMARK_ADVANCED("n=200 p=0.15 r=2")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<3, 2, false>)(g200); });
  };
}

TEST_CASE("bench: vcp<4,2,0> all pairs", "[!benchmark][vcp-4-2-0]") {
  auto const g100 =
      parse_graph<vcp::multirelational_graph<2>>(erdos_renyi_multirelational_text(100, 0.30, 42));
  auto const g200 =
      parse_graph<vcp::multirelational_graph<2>>(erdos_renyi_multirelational_text(200, 0.15, 42));

  BENCHMARK_ADVANCED("n=100 p=0.30 r=2")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 2, false>)(g100); });
  };

  BENCHMARK_ADVANCED("n=200 p=0.15 r=2")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return (vcp_all_pairs_sum<4, 2, false>)(g200); });
  };
}

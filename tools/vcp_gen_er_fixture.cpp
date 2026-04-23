// Generate a sparse Erdős–Rényi random graph fixture for VCP benchmarks.
//
// Produces three outputs in an output directory:
//
//   {dir}/undirected/er_n{N}_d{D}_s{S}.txt       adjacency, one row per vertex
//   {dir}/directed/dsym_n{N}_d{D}_s{S}.txt       hardlink to the undirected file
//   {dir}/pairs/sample_n{N}_d{D}_s{S}_k{K}.pairs uniformly-sampled pair list
//
// The bidirectional-directed fixture is hardlinked from the undirected one
// because the on-disk representation (line per vertex, sorted neighbor IDs)
// is byte-identical under both interpretations: for an undirected edge {u,v}
// each endpoint lists the other, and for a bidirectionally-symmetric directed
// graph each arc u->v plus v->u yields exactly the same per-vertex listing.
// Single-file generation saves ~400 MB of disk and a whole generation pass.
//
// Algorithm (scales with edge count, not vertex count):
//   Pass 1: geometric-gap sampling of edges in ascending-u order; accumulate
//           per-vertex degrees into a flat uint32 array.
//   Pass 2: reseed the RNG, replay the same edge sequence, and write each
//           endpoint into the pre-sized slot in a flat CSR neighbors array.
//   Pass 3: stream per-vertex rows to disk. The natural sort invariant
//           (predecessors in ascending order, followed by ascending-v
//           successors) means no explicit sort step is needed.
//
// This replaces the earlier Python streaming generator. Peak memory is
// ~(2n + 2m) * 4 bytes ~= 320 MB for n=10M, avg_deg=5, and the whole run
// takes a few seconds on modern hardware.

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

namespace {

// Format helpers kept local so the generator stays a single translation unit.
std::string format_stem(std::uint64_t n, std::uint64_t avg_degree, std::uint64_t seed) {
  return "n" + std::to_string(n) + "_d" + std::to_string(avg_degree) + "_s" + std::to_string(seed);
}

// Two-pass CSR-style ER sampling. Writes each vertex's adjacency as a
// space-separated line of neighbor IDs to the given path. Determinism is
// guaranteed by reseeding the RNG identically between pass 1 (degree
// count) and pass 2 (slot fill).
void stream_er_sparse(std::uint64_t n, double p, std::uint64_t seed,
                      std::filesystem::path const &path) {
  if (!(p > 0.0 && p < 1.0)) {
    throw std::runtime_error("p must be in (0, 1); got " + std::to_string(p));
  }
  double const log1mp = std::log1p(-p);

  std::cerr << "  pass 1 (degrees): n=" << n << " p=" << p << "\n";
  std::vector<std::uint32_t> degrees(n, 0);
  {
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (std::uint64_t u = 0; u + 1 < n; ++u) {
      std::uint64_t v = u;
      while (true) {
        double r = u01(rng);
        // u01 produces values in [0, 1); guard against log(0) without
        // biasing the distribution in any observable way.
        if (r <= 0.0) {
          r = std::nextafter(0.0, 1.0);
        }
        auto gap = static_cast<std::uint64_t>(std::log(r) / log1mp) + 1;
        v += gap;
        if (v >= n) {
          break;
        }
        ++degrees[u];
        ++degrees[v];
      }
    }
  }

  std::uint64_t total = 0;
  std::vector<std::uint64_t> offsets(n + 1);
  for (std::uint64_t u = 0; u < n; ++u) {
    offsets[u] = total;
    total += degrees[u];
  }
  offsets[n] = total;
  std::cerr << "  edges counted: m=" << (total / 2) << " (listings=" << total << ")\n";

  std::cerr << "  pass 2 (fill CSR)\n";
  std::vector<std::uint32_t> neighbors(total);
  std::vector<std::uint64_t> cursor(offsets.begin(),
                                    offsets.begin() + static_cast<std::ptrdiff_t>(n));
  {
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (std::uint64_t u = 0; u + 1 < n; ++u) {
      std::uint64_t v = u;
      while (true) {
        double r = u01(rng);
        if (r <= 0.0) {
          r = std::nextafter(0.0, 1.0);
        }
        auto gap = static_cast<std::uint64_t>(std::log(r) / log1mp) + 1;
        v += gap;
        if (v >= n) {
          break;
        }
        neighbors[cursor[u]++] = static_cast<std::uint32_t>(v);
        neighbors[cursor[v]++] = static_cast<std::uint32_t>(u);
      }
    }
  }

  std::cerr << "  pass 3 (stream rows to disk)\n";
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot open " + path.string() + " for writing");
  }
  // Buffer output: ofstream's default buffer is small; a 1 MiB buffer on
  // a 400 MB write keeps us I/O-bound instead of syscall-bound.
  std::vector<char> iobuf(1 << 20);
  out.rdbuf()->pubsetbuf(iobuf.data(), static_cast<std::streamsize>(iobuf.size()));

  for (std::uint64_t u = 0; u < n; ++u) {
    std::uint64_t const start = offsets[u];
    std::uint64_t const end = offsets[u + 1];
    for (std::uint64_t k = start; k < end; ++k) {
      if (k > start) {
        out.put(' ');
      }
      out << neighbors[k];
    }
    out.put('\n');
  }
}

// Sample k distinct unordered pairs (u, v) with u < v, uniformly at random.
// A std::set guarantees distinctness; 100k pairs is trivial for this
// structure, and the birthday-style growth is bounded by k << C(n,2).
void write_sampled_pairs(std::uint64_t n, std::uint64_t k, std::uint64_t seed,
                         std::filesystem::path const &path) {
  if (k > (n * (n - 1)) / 2) {
    throw std::runtime_error("k exceeds C(n, 2)");
  }
  std::cerr << "  sampling " << k << " pairs (seed=" << seed << ")\n";
  std::mt19937_64 rng{seed};
  std::uniform_int_distribution<std::uint64_t> pick(0, n - 1);
  std::set<std::pair<std::uint64_t, std::uint64_t>> pairs;
  while (pairs.size() < k) {
    auto a = pick(rng);
    auto b = pick(rng);
    if (a == b) {
      continue;
    }
    if (a > b) {
      std::swap(a, b);
    }
    pairs.emplace(a, b);
  }

  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot open " + path.string() + " for writing");
  }
  for (auto const &[u, v] : pairs) {
    out << u << ' ' << v << '\n';
  }
}

// Hardlink source -> target. Fails loudly if the link cannot be created;
// does nothing if the target already exists.
void hardlink_if_missing(std::filesystem::path const &source, std::filesystem::path const &target) {
  if (std::filesystem::exists(target)) {
    std::cerr << "  (exists) " << target.filename().string() << "\n";
    return;
  }
  std::filesystem::create_directories(target.parent_path());
  std::error_code ec;
  std::filesystem::create_hard_link(source, target, ec);
  if (ec) {
    throw std::runtime_error("hardlink " + source.string() + " -> " + target.string() +
                             " failed: " + ec.message());
  }
  std::cerr << "  hardlinked " << target.filename().string() << " -> " << source.filename().string()
            << "\n";
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    CLI::App app{"Generate a sparse Erdős–Rényi random graph fixture "
                 "(undirected + bidirectional-directed hardlink + sampled pairs)."};
    app.set_version_flag("--version", std::string("vcp_gen_er_fixture by Ryan N. "
                                                  "Lichtenwalter v") +
                                          VCP_VERSION);

    std::uint64_t n = 10'000'000;
    std::uint64_t avg_degree = 5;
    std::uint64_t seed = 42;
    std::uint64_t pair_sample = 100'000;
    std::string output_dir = "benchmark/fixtures";

    app.add_option("-n,--vertices", n, "Number of vertices")->capture_default_str();
    app.add_option("-d,--avg-degree", avg_degree, "Target average degree (integer; p = d / (n-1))")
        ->capture_default_str();
    app.add_option("-s,--seed", seed, "RNG seed")->capture_default_str();
    app.add_option("-k,--pair-sample", pair_sample, "Number of sampled pairs")
        ->capture_default_str();
    app.add_option("-o,--output-dir", output_dir,
                   "Root output directory; creates undirected/, directed/, pairs/ "
                   "subdirectories")
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    if (n < 2) {
      std::cerr << "error: n must be >= 2\n";
      return 1;
    }
    if (avg_degree == 0 || avg_degree >= n - 1) {
      std::cerr << "error: avg-degree must be in (0, n-1)\n";
      return 1;
    }

    double const p = static_cast<double>(avg_degree) / static_cast<double>(n - 1);
    std::string const stem = format_stem(n, avg_degree, seed);

    std::filesystem::path base(output_dir);
    std::filesystem::path undirected_path = base / "undirected" / ("er_" + stem + ".txt");
    std::filesystem::path directed_path = base / "directed" / ("dsym_" + stem + ".txt");
    std::filesystem::path pairs_path =
        base / "pairs" / ("sample_" + stem + "_k" + std::to_string(pair_sample) + ".pairs");

    // Undirected (the primary, source-of-truth graph).
    if (std::filesystem::exists(undirected_path)) {
      std::cerr << "  (exists) " << undirected_path.filename().string() << "\n";
    } else {
      std::cerr << "generating " << undirected_path.filename().string() << "\n";
      stream_er_sparse(n, p, seed, undirected_path);
    }

    // Directed (bidirectional -> byte-identical to undirected, hardlink).
    hardlink_if_missing(undirected_path, directed_path);

    // Sampled pairs. Use seed+1 so pair RNG stream is disjoint from graph RNG.
    if (std::filesystem::exists(pairs_path)) {
      std::cerr << "  (exists) " << pairs_path.filename().string() << "\n";
    } else {
      std::cerr << "generating " << pairs_path.filename().string() << "\n";
      write_sampled_pairs(n, pair_sample, seed + 1, pairs_path);
    }

    return 0;
  } catch (std::exception const &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}

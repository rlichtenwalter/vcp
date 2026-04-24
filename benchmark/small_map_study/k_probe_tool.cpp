// k-probe harness — runs vcp<4, r, d>::generate_vector over a pair file and
// lets the VCP_INSTRUMENT_K hook record the observed cardinality of the
// three hot maps (temp_edge_types, edge_types, counts) per call.
//
// This tool is NOT a public CLI; it exists only for the small_map study
// and is built when `-DVCP_BUILD_K_PROBE=ON` is passed to cmake. The
// instrumentation header is `benchmark/small_map_study/k_probe.hpp`
// (added to the include path by this target's CMakeLists contribution).
//
// Usage:
//   VCP_K_PROBE_OUTPUT=/tmp/k_probe_<label>.csv
//     ./build/benchmark/small_map_study/k_probe_tool
//       --workload <label> -r <2|30> -d <0|1>
//       --graph <path> --pairs <path>
//
// Output CSV columns:
//   map_name, workload, bucket_log2, bucket_count, bucket_lo, bucket_hi,
//   total_samples, max_observed
//
// The tool does not print the VCP vectors — the point is to record k, not
// to compute the profile. Sinking to stdout would dominate the runtime on
// large fixtures and does not contribute to the study.

#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include <CLI/CLI.hpp>

#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/vcp.hpp>

#include "k_probe.hpp"

namespace {

template <std::size_t r, bool d>
void probe_multirelational(std::string const &graph_file, std::string const &pairs_file) {
  using graph_type =
      std::conditional_t<d, vcp::multirelational_directed_graph<r>, vcp::multirelational_graph<r>>;
  std::ifstream gfs(graph_file);
  if (!gfs) {
    std::cerr << "cannot open graph: " << graph_file << "\n";
    std::exit(1);
  }
  graph_type g;
  gfs >> g;

  std::ifstream pfs(pairs_file);
  if (!pfs) {
    std::cerr << "cannot open pairs: " << pairs_file << "\n";
    std::exit(1);
  }

  vcp::vcp<4, r, d> profiler(g);
  vcp::vertex_id_t v1;
  vcp::vertex_id_t v2;
  std::size_t pair_count = 0;
  while (pfs >> v1 >> v2) {
    (void)profiler.generate_vector(vcp::const_vertex_iterator(g.vertices_begin() + v1),
                                   vcp::const_vertex_iterator(g.vertices_begin() + v2));
    ++pair_count;
  }
  std::cerr << "processed " << pair_count << " pairs\n";
}

} // namespace

int main(int argc, char *argv[]) {
  CLI::App app{"k-probe harness for vcp<4, r, d>::generate_vector hot maps"};
  std::string workload;
  std::size_t r = 0;
  std::size_t d = 0;
  std::string graph_file;
  std::string pairs_file;

  app.add_option("--workload", workload, "Workload label (recorded in CSV)")->required();
  app.add_option("-r", r, "Number of relations (2 or 30)")
      ->required()
      ->check(CLI::IsMember({2, 30}));
  app.add_option("-d", d, "Directedness (0 or 1)")->required()->check(CLI::IsMember({0, 1}));
  app.add_option("--graph", graph_file, "Graph file")->required()->check(CLI::ExistingFile);
  app.add_option("--pairs", pairs_file, "Pairs file")->required()->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  k_probe::set_workload(workload);

  if (r == 2 && d == 0) {
    probe_multirelational<2, false>(graph_file, pairs_file);
  } else if (r == 2 && d == 1) {
    probe_multirelational<2, true>(graph_file, pairs_file);
  } else if (r == 30 && d == 0) {
    probe_multirelational<30, false>(graph_file, pairs_file);
  } else if (r == 30 && d == 1) {
    probe_multirelational<30, true>(graph_file, pairs_file);
  }

  return 0;
}

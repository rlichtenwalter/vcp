# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- CMake build system (minimum 3.21) replacing the hand-written Makefile
- CMake install support with `find_package(vcp)` for downstream consumers
- pkg-config support for non-CMake consumers
- VERSION file, sourced into the CMake project version
- CTest integration with a Catch2 v3 test scaffold (real tests added in a subsequent release)
- clang-format configuration (LLVM base style, 100-column)
- clang-tidy configuration for static analysis
- `.pre-commit-config.yaml` with file hygiene, clang-format, and branch protection hooks
- Gitea Actions CI workflow: build-and-test (Release/Debug matrix), quality (pre-commit), lint (clang-tidy)
- CLI `--version` output for each tool, sourced from VERSION via CMake

### Changed
- **BREAKING**: Move headers from `inc/vcp/` to `include/vcp/` (use `#include <vcp/vcp.hpp>`)
- **BREAKING**: Move CLI tool sources from `src/` to `tools/` following header-only library conventions
- **BREAKING**: Minimum C++ standard raised from C++11 to C++14
- Replace vendored Boost 1.53 with FetchContent of a modern Boost.Multiprecision
- Replace vendored TCLAP with CLI11 for command-line argument parsing
- Rename `COPYING` to `LICENSE` (GPL-3.0 content unchanged)
- Rename `CHANGES` to `CHANGELOG.md` with Keep-a-Changelog formatting

### Removed
- `build_local_gcc.sh` — obsolete gcc-4.8 bootstrap script from 2013
- `AUTHORS` — redundant with license header and source-file copyright notices
- Vendored `lib/boost_1_53_0/` (~23 MB, 2013 vintage)
- Vendored `lib/tclap-1.2.1/`

### Fixed
- `multirelational_graph::edge_value` and `multirelational_directed_graph::edge_value` now return 0 (unconnected) when passed the `edges_end()` sentinel, matching the contract of the simple `graph` and `directed_graph` siblings. Callers in `vcp::edge_value()` pipe `edge()` results directly into `edge_value()` without a null check; the previous unconditional array lookup caused an out-of-bounds read on any multirelational graph with missing edges (e.g. paths).
- `directed_to_undirected` no longer infinite-loops on directed graphs whose out-neighbor and in-neighbor sets differ. The two tail merge loops were missing their iterator advances. Also skip the output block for isolated vertices (empty neighbor list), whose `neighbors.end() - 1` formed a pointer before `begin()`.
- `multirelational_graph::operator<<` and `multirelational_directed_graph::operator<<` no longer compute `neighbors_end(vIt) - 1` on an isolated vertex at the start of the shared edges array, which formed a pointer before the array base (undefined behavior per C++ [expr.add]/4, flagged by UBSan's pointer-overflow check). Both operators are rewritten using the separator-before-non-first idiom, producing byte-identical output on all non-pathological cases.

## [1.0.0] - 2012-05-01

### Added
- Initial commit of source files to the repository
- Header-only library for Vertex Collocation Profile computation
- Compressed sparse row (CSR) graph classes: `graph`, `directed_graph`, `multirelational_graph`, `multirelational_directed_graph`
- VCP algorithm template with partial specializations for (n=3,r=1) and (n=4,r=1), directed and undirected
- Static and dynamic subgraph-to-element mappers
- CLI tools: `vcp_generate`, `vcp_map`, `directed_to_undirected`, `ell_2_pairs`

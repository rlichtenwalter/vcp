# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- Pre-commit `mixed-line-ending` hook now normalises every commit to LF via `args: [--fix=lf]`

## [2.0.0] - 2026-05-05

### Added
- `check-json` pre-commit hook validating `CMakePresets.json` and any future JSON files at commit time
- Gitea Actions workflow mirroring Gitea releases to GitHub, closing the push mirror's release-metadata gap
  - Idempotent, with a `workflow_dispatch` path (`tag` input) for manual backfill
- Ignore `.env` and `.env.*` in `.gitignore` while still permitting a committed `.env.example` template
- Doxygen docstrings on all 17 public headers, covering every public class, typedef, constructor, and method
  - Thread-safety contracts and per-specialization element-count invariants called out explicitly
- Unit-test coverage backfill for foundation classes and specializations that previously had zero direct coverage; ctest grows from 60 to 118
- `detail::dense_or_sparse_map` — two-tier compile-time-dispatched container replacing `std::map` for the `vcp<4, r, d>` edge-type maps
  - Dense tier under an 8 MB byte budget (r ≤ 20 undirected, r ≤ 10 directed); `std::unordered_map` fallback above
- `test/test_dense_or_sparse_map.cpp` covering tier selection, insert-or-zero semantics, and both pair-keyed paths
- `benchmark/small_map_study/` — opt-in (`-DVCP_BUILD_SMALL_MAP_STUDY=ON`) microbenchmark harness and recorded results behind the container-strategy decision
- `benchmark/counts_map_study/` — opt-in (`-DVCP_BUILD_COUNTS_BENCH=ON`) study of six candidate `counts`-map replacements; findings in `benchmark/counts_map_study/results/summary.md`
- `vcp_gen_er_fixture` CLI tool (`tools/vcp_gen_er_fixture.cpp`) generating sparse Erdős–Rényi fixtures (~6 s at n=10M vs ~1-2 min for the earlier Python prototype)
- Large-tier benchmark workload (`./benchmark/run.sh --large`) on a 10M-node sparse ER graph, exercising the four unirelational VCP procedures
- `regression/run.sh` Phase L runs legacy-vs-current byte-diff checks on the benchmark's 10M-node fixture, with legacy's known `vcp 4 1 1` divergence encoded as expected
- Benchmark suite under `benchmark/`: a cross-revision tool-level runner (`benchmark/run.sh`) and Catch2 micro-benchmarks gated on `VCP_BUILD_BENCHMARKS`
  - The runner builds arbitrary git revisions in isolated worktrees and compares against the 2012 baseline by default
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

- `CITATION.cff` with the open-access 2014 SpringerPlus paper as preferred citation; Gitea and GitHub render a "Cite this repository" affordance from it
- `CMakePresets.json` with `release`, `debug`, and `sanitize` configure/build/test presets, each keeping its own warm cache under `build/<name>`
  - The `sanitize` test preset carries the `ASAN_OPTIONS`/`UBSAN_OPTIONS` halt-on-error contract previously duplicated inline in CI yaml

### Changed
- `VCP_WARNING_FLAGS` expanded with `-Wconversion -Wsign-conversion -Wshadow -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough` plus GCC-only `-Wlogical-op` and `-Wduplicated-cond`
  - Mechanical fallout (explicit casts, parameter renames) fixed across nine headers; no semantics altered
- `.gitea/workflows/ci.yml` invokes named CMake presets instead of inline flags, so CI reproduces a developer's local `cmake --preset=...` invocation byte-for-byte
- `.gitignore` simplified: the `build-*/` glob is removed — presets place all per-config trees under the already-ignored `build/`
- CI `build-and-test` matrix extended with Clang; both GCC and Clang now build everything and run the full ctest suite at Release and Debug
- Test and benchmark targets now compile with the same warning flags as the CLI tools, via a new shared `VCP_WARNING_FLAGS` CMake variable
- `std::make_pair(...)` replaced with the C++17 CTAD `std::pair{...}` form throughout headers, tools, and tests (~50 call sites)
- All four graph foundation classes declare defaulted noexcept move constructor/assignment, so `std::vector<graph>` reallocations move instead of copying the heap CSR arrays
- `[[nodiscard]]` added to every observer on the public API surface, including the nullable-iterator find-edge family
- `noexcept` added to one-line accessors on the four graph foundation classes, `square_matrix::size`, and the new defaulted moves
- Top-level `const` dropped from `generate_vector` return types (18 declaration/definition pairs), letting call sites move instead of copy the returned map
- The primary `vcp<n, r, d>` template now carries `requires (n >= 2 && r >= 1)`, turning nonsense instantiations into a single-line diagnostic at the use site
- **BREAKING**: CMake minimum requirement raised from 3.21 to 3.24; all current target distros ship >= 3.24 in their default repositories
- **BREAKING**: Minimum C++ standard raised from C++14 to C++20
- **BREAKING**: License changed from GPL-3.0-or-later to BSD 3-Clause
  - Versions before 2.0.0 remain available under their original GPL-3.0-or-later terms
- The 16-line GPL boilerplate in every public header and tool source (22 files) replaced with a two-line SPDX form; the full BSD 3-Clause text remains authoritative in `LICENSE`
- Documented thread-safety contract per `vcp<n, r, d>` specialization in README and in-header comments ("one `vcp` per thread" is always safe and recommended)
- `regression/run.sh` now encodes expected legacy-vs-current divergence for fixture classes where legacy has a known, since-fixed bug; full regression is 165/165 PASS
- `regression/validate_bug2_invariants.py`: dropped a redundant `strict=False` on a `zip()` call (3.10+ default).
- `tools/vcp_gen_er_fixture.cpp`: flush-and-check after both write loops, so a mid-write I/O failure raises instead of silently leaving a truncated fixture
- Public-header includes narrowed: `<iostream>` dropped or replaced with `<istream>`/`<ostream>` in the eight headers that over-included it
- Renamed top-level `script/` directory to `scripts/` for plural consistency with `tools/`; repo-layout change only
- Standardized include-guard convention to `VCP_<NAME>_HPP` across all 17 headers, matching sibling libraries `kdtree` and `mRMR`
- `std::make_unique<T[]>` replaces the explicit `new T[N]` form at every owning-array allocation site, removing every explicit `new`/`delete` from the public headers
- `vcp<4, 1, false>::generate_vector` hot path: stride-based V3→V4 re-encoding replaces two chained equality tests — a 14–18% wall-clock improvement on `vcp_generate 4 1 0`
- `VCP_SANITIZE` now enables ASan and UBSan on every built target with `-fno-sanitize-recover=all`; the new CI `sanitize` job runs this configuration on every PR
- **BREAKING**: Move headers from `inc/vcp/` to `include/vcp/` (use `#include <vcp/vcp.hpp>`)
- **BREAKING**: Move CLI tool sources from `src/` to `tools/` following header-only library conventions
- **BREAKING**: Minimum C++ standard raised from C++11 to C++14
- Replace vendored Boost 1.53 with FetchContent of a modern Boost.Multiprecision
- Replace vendored TCLAP with CLI11 for command-line argument parsing
- Rename `COPYING` to `LICENSE` (GPL-3.0 content unchanged)
- Rename `CHANGES` to `CHANGELOG.md` with Keep-a-Changelog formatting
- Migrated `std::copy`/`std::max_element` call sites to `std::ranges` equivalents and restructured two tests so clang-tidy 20.x runs clean with no behavior change

### Removed
- Redundant duplicated `vertex_id_t`/`edge_id_t`/iterator `using` declarations from the three derived graph headers; `graph.hpp` is now their single source of truth
- Dead template-metafunction primitives `vcp::TMP_power`, `vcp::TMP_min`, and `vcp::TMP_max` in `vcp_dynamic_mapper.hpp`; unreferenced, removed without replacement
- Redundant naked forward declarations of the free `vcp::edge` and `vcp::edge_value` overloads in `vcp.hpp`
- `build_local_gcc.sh` — obsolete gcc-4.8 bootstrap script from 2013
- `AUTHORS` — redundant with license header and source-file copyright notices
- Vendored `lib/boost_1_53_0/` (~23 MB, 2013 vintage)
- Vendored `lib/tclap-1.2.1/`
- `VCP_INSTRUMENT_K` instrumentation hook and its one-off k-probe harness; the recorded cardinality results remain under `benchmark/small_map_study/results/`
  - CMake gate for the remaining microbenchmark renamed from `VCP_BUILD_K_PROBE` to `VCP_BUILD_SMALL_MAP_STUDY`

### Fixed
- Removed the compile-time `MAX_NEIGHBORS` cap; `vcp<4, r, d>` now sizes the `v3Vertices` scratch buffer to an exact graph-derived bound, eliminating a heap overflow risk under `NDEBUG`
  - `detail/dense_or_sparse_map.hpp` and the new `detail/v3_buffer_bound.hpp` are now in `VCP_PUBLIC_HEADERS`, fixing a latent install gap
- `vcp_static_mapper` computes integer powers via bit shifts instead of `std::pow`, removing `double` round-trip inexactness; oversized shifts now throw `std::length_error`
- `square_matrix::size()` no longer recovers the side length via `std::sqrt`, removing a latent off-by-one for n² ≥ 2²⁵ and making `size()` constant-time on both specializations
- Removed redundant `typename` preceding non-type template parameters in `vcp.hpp` and `square_matrix.hpp`
- `directed_graph` out/in edge-iterator definitions now return `const_edge_iterator` to match their declarations; both aliases resolve identically, so behavior is unchanged
- `vcp<3, r, false>::generate_vector` reads the correct edge value for v3 candidates exclusive to v2; at r ≥ 2 the wrong value routed subgraphs to a wrong canonical bucket (present since 2012)
- Four latent clang-tidy errors, and the CI `lint` gap that excluded `include/vcp/detail/*.hpp` from scanning
- `detail::pair_hash` uses the 64-bit hash-combine constant on 64-bit platforms, closing a collision risk in the sparse tier
- Debug-only assertions guard the unsigned-arithmetic tail of `vcp<4, 1, d>::generate_vector` against silent underflow; compiled away under `NDEBUG`, exercised by the CI sanitize build
- Graph istream operators parse vertex ids with `std::from_chars` instead of `atol`; malformed non-empty tokens now throw `std::invalid_argument` instead of silently reading 0
- Merge loops in `vcp_4_1_1.hpp` and `vcp_4_r_1.hpp` use explicit exhaustion sentinels instead of an iterator-argument coincidence; byte-identical output
- `vcp_dynamic_mapper::canonical_subgraph_address` no longer double-counts pairs for undirected symmetric-matrix input, which made canonicalization a silent no-op (latent for current callers)
- `vcp_dynamic_mapper::element_structure`: corrected the per-slot shift (`>>= r`, not `>>= r_pset`) and an ill-formed const-ref mutation; the function had never been instantiated
- `multirelational_*_graph::relation_count`: removed a function-local `static` that bled the first call's result across instances, and guarded `std::max_element` against empty edges (UB)
- `square_matrix<T, 0>` fixed-to-dynamic conversion paths: the copy constructor dropped the last cell and copy assignment under-sized its resize (heap overflow); both latent until now
- `vcp<4, r, true>::generate_vector`: a wrong-slot diagonal write misrouted v3-v4-connected subgraphs, and a missing post-enumeration self-correction inflated counts for connected v1v2
- `vcp<4, r, false>::generate_vector` reads the correct edge value for v2-exclusive v3 candidates; two copy-paste defects misrouted subgraphs for all r ≥ 2 undirected graphs (present since 2012)
- Pad `dbug2_*` regression fixtures with trailing empty lines so the line-count-based graph parser reads the correct `num_vertices`
- `vcp<4, 1, true>::generate_vector` produces correct output on mutual-edge and shared-neighbor inputs, fixing three co-located counting/addressing defects
  - New regression Phase 6 (`validate_bug2_invariants.py`) asserts ground-truth invariants across six new `dbug2_*` fixtures
- `multirelational_*_graph::edge_value` returns 0 for the `edges_end()` sentinel instead of an out-of-bounds read, matching the simple graph siblings' contract
- `directed_to_undirected` no longer infinite-loops when out- and in-neighbor sets differ, and skips isolated vertices whose output block formed a pointer before `begin()`
- Graph `operator<<` overloads no longer form a pre-`begin()` pointer on isolated vertices (UB); rewritten with the separator-before-non-first idiom, byte-identical output

## [1.0.0] - 2012-05-01

### Added
- Initial commit of source files to the repository
- Header-only library for Vertex Collocation Profile computation
- Compressed sparse row (CSR) graph classes: `graph`, `directed_graph`, `multirelational_graph`, `multirelational_directed_graph`
- VCP algorithm template with partial specializations for (n=3,r=1) and (n=4,r=1), directed and undirected
- Static and dynamic subgraph-to-element mappers
- CLI tools: `vcp_generate`, `vcp_map`, `directed_to_undirected`, `ell_2_pairs`

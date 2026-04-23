# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `vcp_gen_er_fixture` CLI tool (`tools/vcp_gen_er_fixture.cpp`) — C++ generator for sparse Erdős–Rényi fixtures using a two-pass CSR-style sampling procedure (degree count, then RNG-replayed slot fill, then streamed row emission). Produces the undirected graph, a byte-identical hardlink for the bidirectional-directed variant, and a uniform-sample pairs file in a single invocation. Completes in ~6 s at n=10M vs the earlier Python prototype's ~1-2 min, with peak RSS ~400 MB instead of a risk of OOM on set-based adjacency. Same tool can be reused for any future sparse random-graph fixture need.
- Large-tier benchmark workload (`./benchmark/run.sh --large`) on a 10M-node avg-degree-5 sparse ER graph, exercising the four unirelational VCP procedures (3/4 × undirected/directed). Fixture generation is delegated to `vcp_gen_er_fixture`.
- `regression/run.sh` Phase L auto-detects the benchmark's large fixture and runs legacy-vs-current byte-diff consistency checks on the same 10M-node graph. Known divergence on `vcp 4 1 1` (legacy has an unfixed underflow) is encoded as an expected result; identical output on 3/4 × {undirected, directed} is the PASS criterion.
- Benchmark suite under `benchmark/` with two layers:
  - Cross-ref tool-level runner (`benchmark/run.sh`) that builds arbitrary git revisions in isolated worktrees, times the CLI tools against a curated set of seeded fixtures, and emits a markdown comparison table. Auto-detects legacy (Makefile) vs modern (CMake) build systems and prepends the 2012 baseline commit by default so current-branch performance can be measured directly against the original GitHub state.
  - Catch2 library-level micro-benchmarks (`benchmark/bench_vcp.cpp`) gated on the new `VCP_BUILD_BENCHMARKS` CMake option, mirroring the `bench_kdtree.cpp` / `bench_mrmr.cpp` pattern in sibling libraries.
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
- `VCP_SANITIZE` now enables AddressSanitizer **and** UndefinedBehaviorSanitizer (previously only ASan), applies to every built target (tools, tests, benchmarks — previously only CLI tools), and passes `-fno-sanitize-recover=all` so every sanitizer diagnostic is a hard error. The new CI `sanitize` job builds this configuration and runs the test suite on every PR.
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
- `vcp<4, 1, true>::generate_vector` now produces correct output on inputs containing mutual edges or shared-neighbor structure. The directed 4-vertex specialization had three co-located defects in the post-enumeration accounting:
  - **Counting (line 579)**: the asymmetric-v3v4 slot over-subtracted the v₁v₂ edge by 1 when v₁v₂ was mutual, causing unsigned underflow to 2⁶⁴−1. The original `static_cast<bool>(v1v2)` treated any v₁v₂ edge as asymmetric, but the v₁v₂ pair contributes to `amutualPairs` only when genuinely asymmetric. Fix: subtract 1 iff `v1v2` is OUT or IN.
  - **Counting (line 581)**: the mutual-v3v4 slot was missing the symmetric correction, inflating the count by 1 when v₁v₂ was mutual. The v₁v₂ mutual edge lives in `mutualPairs` but is absent from the `connections − amutuals` observation term, so the correction belongs here. Fix: subtract 1 iff `v1v2` is BOTH.
  - **Addressing (lines 538, 547)**: the re-encoding from V₁₃/V₂₃ (v₃ role) to V₁₄/V₂₄ (v₄ role) dropped the V₁₃ unit-conversion factor, inflating the V₁₄ slot value by V₁₃ = 4. In the worst case (both c₁₃ and c₂₃ = BOTH) the generated address reached 4239, overflowing the 4096-entry `element_address` static map (out-of-bounds read in `.rodata`, followed by a write to an arbitrary `counts[]` slot). Fix: exploit the stride equivalence V₁₄/V₁₃ = V₂₄/V₂₃ = 4 and compute `contrib = 4 * temp` directly; the explicit form is documented in-place for future-proofing against enum-layout changes.
  - The three defects coincide on K₄-bidirectional, which is now the minimal regression fixture. Added a new Phase 6 to `regression/run.sh` (`validate_bug2_invariants.py`) that asserts conservation (`sum(counts) == C(V−2, 2)`), non-underflow (`max(counts) <= C(V−2, 2)`), and slot localization (`counts[expected_element] == 1`) across six new directed fixtures (`dbug2_*`) that partition the (v₁v₂ type × shared-v₃ pattern) trigger space. Phase 6 reports absolute correctness independent of the byte-level legacy-vs-current diff used elsewhere.
- `multirelational_graph::edge_value` and `multirelational_directed_graph::edge_value` now return 0 (unconnected) when passed the `edges_end()` sentinel, matching the contract of the simple `graph` and `directed_graph` siblings. Callers in `vcp::edge_value()` pipe `edge()` results directly into `edge_value()` without a null check; the previous unconditional array lookup caused an out-of-bounds read on any multirelational graph with missing edges (e.g. paths).
- `directed_to_undirected` no longer infinite-loops on directed graphs whose out-neighbor and in-neighbor sets differ. The two tail merge loops were missing their iterator advances. Also skip the output block for isolated vertices (empty neighbor list), whose `neighbors.end() - 1` formed a pointer before `begin()`.
- `graph::operator<<`, `directed_graph::operator<<`, `multirelational_graph::operator<<`, and `multirelational_directed_graph::operator<<` no longer compute `neighbors_end(vIt) - 1` on an isolated vertex at the start of the shared edges array, which formed a pointer before the array base (undefined behavior per C++ [expr.add]/4, flagged by UBSan's pointer-overflow check). All four operators are rewritten using the separator-before-non-first idiom, producing byte-identical output on all non-pathological cases.

## [1.0.0] - 2012-05-01

### Added
- Initial commit of source files to the repository
- Header-only library for Vertex Collocation Profile computation
- Compressed sparse row (CSR) graph classes: `graph`, `directed_graph`, `multirelational_graph`, `multirelational_directed_graph`
- VCP algorithm template with partial specializations for (n=3,r=1) and (n=4,r=1), directed and undirected
- Static and dynamic subgraph-to-element mappers
- CLI tools: `vcp_generate`, `vcp_map`, `directed_to_undirected`, `ell_2_pairs`

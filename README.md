# vcp — Vertex Collocation Profiles

Author: Ryan N. Lichtenwalter (<rlichtenwalter@gmail.com>)

Header-only C++ library and companion CLI tools for computing Vertex Collocation
Profiles (VCPs) over graphs. The library provides processor- and memory-efficient
graph representations (compressed sparse row), VCP algorithm template classes
with hand-tuned specializations, and static/dynamic subgraph-to-element mappers.

## Contents

The code base has two parts:

1. **A header-only library** in `include/vcp/` with:
   - Efficient graph classes (`graph`, `directed_graph`, `multirelational_graph`,
     `multirelational_directed_graph`) backed by compressed sparse row storage.
   - VCP computation templates with specializations for common parameter
     combinations, plus naive fallbacks.
   - Static and dynamic subgraph-to-element mappers.
2. **Four CLI tools** (in `tools/`):
   - `vcp_generate` — generate VCP vectors for vertex pairs read from stdin.
   - `vcp_map` — output the subgraph-to-element mapping for a given VCP.
   - `directed_to_undirected` — convert a directed graph to undirected form.
   - `ell_2_pairs` — enumerate all vertex pairs at distance 2 from each other.

All tools accept `--help` and `--version`.

## Building

Requires a C++20 compiler, CMake ≥ 3.24, and Boost headers (for
`boost::multiprecision`, used to support multirelational VCPs whose address
space exceeds 64 bits). On Debian/Ubuntu:

```bash
sudo apt-get install libboost-dev
```

Configure and build:

```bash
cmake -B build
cmake --build build
```

Debug build (assertions and debug symbols):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

The tools will be placed in `build/` after a successful build.

To force a fresh configure (drop the cached CMake state and reconfigure
from scratch — useful after changing the C++ standard, toolchain, or
build options):

```bash
cmake -B build --fresh
```

The `MAX_NEIGHBORS` compile-time cap (default: 16384) can be adjusted at
configure time:

```bash
cmake -B build -DVCP_MAX_NEIGHBORS=32768
```

## Testing

```bash
ctest --test-dir build --output-on-failure
```

### Sanitized build

Set `-DVCP_SANITIZE=ON` with a Debug build to enable AddressSanitizer and
UndefinedBehaviorSanitizer across every target (library consumers, tools,
tests, and — when combined with `-DVCP_BUILD_BENCHMARKS=ON` — benchmarks):

```bash
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug -DVCP_SANITIZE=ON
cmake --build build-san
ASAN_OPTIONS=halt_on_error=1:detect_leaks=1:abort_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:abort_on_error=1:print_stacktrace=1 \
ctest --test-dir build-san --output-on-failure
```

The `sanitize` CI job runs this combination on every PR, so UB and memory
errors surface as failures rather than reaching `main`. The
`-fno-sanitize-recover=all` flag makes every sanitizer diagnostic a hard
error; Release builds are never affected.

## Library usage

Because the library is header-only, you do not need to build anything to use it
from your own code — just add `include/` to your include path and link against
Boost headers. For example:

```bash
g++ -std=c++20 -I /path/to/vcp/include myfile.cpp
```

In `myfile.cpp`:

```cpp
#include <vcp/graph.hpp>
#include <vcp/vcp.hpp>
```

### Available headers

All headers live under `include/vcp/`.

- `graph.hpp` — processor- and memory-efficient unirelational graph class in CSR
  format. Supports both undirected and directed graphs; directed graphs have
  efficient edge access in only one direction.
- `directed_graph.hpp` — bidirectional directed graph with fast edge access in
  either direction.
- `multirelational_graph.hpp` — as `graph.hpp`, plus a per-edge integral-encoded
  bitset tracking which relations are present.
- `multirelational_directed_graph.hpp` — as `directed_graph.hpp`, with
  multirelational support.
- `vcp.hpp` — top-level header including all VCP algorithms and specializations.
  Include this when using any VCP algorithm. Provides a naive algorithm that
  generalizes to any `(n, r, d)`, almost always overridden by a more efficient
  specialization.
- `vcp_<n>_<r>_<d>.hpp` — full or partial specializations of the VCP template
  for particular `(n, r, d)` values.
- `vcp_static_mapper.hpp` — general VCP subgraph-to-element mapper. The mapping
  table is generated at construction and must fit in memory; implementations
  are extremely fast but memory requirements become problematic for large
  combinations of `n` and `r`.
- `vcp_dynamic_mapper.hpp` — dynamic subgraph-to-element mapping with no upfront
  table. Requires compile-time `n`, `r`, `d` for maximal efficiency and
  minimal memory. Scales well for sparsely populated VCPs with large `n` and
  `r`.
- `square_matrix.hpp` — statically or dynamically allocated square matrix; with
  `n == 0` at compile time, the allocation is dynamic with runtime-determined
  size.

## Installation

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build
cmake --install build
```

After installation, downstream CMake projects can:

```cmake
find_package(vcp REQUIRED)
target_link_libraries(your_target PRIVATE vcp::vcp)
```

Or via pkg-config:

```bash
pkg-config --cflags vcp
```

## Design notes

Several design decisions that may not be obvious from the source:

- **Header-only.** The library is header-only rather than split into interface
  and implementation. Four reasons: (1) extensive templating requires
  implementations to be visible in the same translation unit as the headers,
  (2) header-only libraries compose easily with other projects without linking,
  (3) the project is open source and there is no reason to separate
  implementation from interface, and (4) benchmarks showed small but
  consistent runtime improvements over a non-header-only build. Interface
  clarity is preserved by defining member functions outside the class body
  within each header.

- **Multiplication/addition over bit shifting.** The classes `vcp_3_1_0`,
  `vcp_3_1_1`, `vcp_4_1_0`, and `vcp_4_1_1` use multiplications and additions
  instead of shifts and bitwise ORs when composing canonical subgraph
  addresses. On x86-64 with barrel shifters, shifts were consistently a
  fractional percentage slower in benchmarks — likely due to more ALU adders
  being available under deep superscalar pipelining.

- **Return-by-value of `std::array`.** The `generate_vector` members of the
  `vcp_3_1_0`, `vcp_3_1_1`, `vcp_4_1_0`, and `vcp_4_1_1` classes return a
  `std::array` by value. Benchmarks showed this was faster than, or equal to,
  writing to an iterator range, a `std::array` reference, or a `std::vector`
  reference — RVO keeps the actual storage on the caller's stack.

- **Thread safety.** The rule differs by specialization:

  - `vcp<3, r, d>` — **safe** for shared-instance concurrent calls. The
    class holds only a `const` reference to the graph; all per-call state
    in `generate_vector` is stack-local.
  - `vcp<4, r, d>` — **not safe** for shared-instance concurrent calls.
    The class holds a `v3Vertices` scratch buffer allocated in the
    constructor and reused across every call; concurrent calls would race
    on the scratch and produce incorrect output.

  Regardless of specialization, the underlying graph classes (`graph`,
  `directed_graph`, `multirelational_graph`, `multirelational_directed_graph`)
  are safe to share across threads for read-only access, provided no
  thread is concurrently calling a mutating `operator>>`.

  The **uniform pattern** — one `vcp<n, r, d>` instance per thread, shared
  graph — is always safe and is the recommended default. Construction is
  cheap; the only per-instance cost is the scratch allocation for the
  `vcp<4, ...>` specializations (256 KB to 2 MB depending on `r`).

  A future revision may expose a thread-safe scratch-handling option
  (caller-managed scratch, `thread_local`, or similar) for the `vcp<4, ...>`
  specializations if a concrete multithreaded use case emerges. The
  `vcp<3, ...>` specializations need no such change.

- **Container strategy for `temp_edge_types` and `edge_types`** in the
  `vcp<4, r, false>` and `vcp<4, r, true>` generic specializations.
  These two internal maps used to be `std::map<connectivity_address_type,
  unsigned long>` (and a pair-keyed variant for the directed case).
  They now route through `detail::dense_or_sparse_map`, a two-tier
  facade that picks the container at compile time from the key-space
  size:

  - **tier 1 (dense)** — `std::array<unsigned long, KeyCount>` backed
    by a heap allocation plus a small presence-tracking vector. Fires
    when `KeyCount * sizeof(unsigned long)` fits under an 8 MB byte
    budget. Covers `r ≤ 20` (undirected) and `r ≤ 10` (directed,
    where the key is a sorted pair and the effective bit count is
    `2r`). Construction is O(1) in real work — `new` reserves virtual
    memory without touching physical pages — and each call is zero
    per-entry allocations because the container is a class member
    reused across `generate_vector` calls, cleared in O(k) at the top
    of each call.

  - **tier 2 (sparse)** — `std::unordered_map` with a `pair_hash`
    combine for the directed-edge pair key. Fires above the byte
    budget (so `r ≥ 21` undirected or `r ≥ 11` directed). Avoids the
    per-node allocator cost `std::map` would incur at the cardinalities
    observed in practice (up to k ≈ 3000 on 30-bit diverse edge
    values; see `benchmark/small_map_study/results/summary_phase_B.md`).

  A three-tier design with a `std::vector<std::pair>` linear-scan
  middle tier was evaluated in
  `benchmark/small_map_study/results/summary_phase_C.md` and dropped:
  at every cardinality we measured, linear scan was beaten by either
  the dense tier (when the key space was small enough that dense was
  feasible) or by the hash tier (when k ≥ ~16). The linear-scan tier
  had no operating regime where it was strictly best, so carrying it
  would have been complexity without performance.

  `counts` — the per-call return map — is intentionally **not**
  refactored by this change. Its keys are `subgraph_address_type`,
  which becomes `boost::multiprecision::cpp_int` for `r ≥ 11`
  undirected (or `r ≥ 6` directed), making a dense tier infeasible
  and requiring a cpp_int-compatible hash for a sparse tier. Left for
  a follow-up once the larger-k cost distribution on representative
  workloads is characterized.

  Tool-level measurement: `vcp_gen 4 2 0 n=200` (the workload that
  stresses the refactored maps most in the default instantiation set)
  improves from 2.11s to 2.09s median of 5 runs, a modest 2–3%. The
  architectural benefit — zero per-call allocations for the two
  refactored maps at every r — is the more durable win and is the
  reason to prefer the new design even where the wall-clock delta is
  small. See `benchmark/small_map_study/` for the full study.

## Dependencies

- **Boost** (headers only) — for `boost::multiprecision::cpp_int`, used to
  represent multirelational address spaces that exceed 64 bits.
- **[CLI11](https://github.com/CLIUtils/CLI11)** — for command-line parsing in
  the CLI tools. Fetched automatically by CMake.
- **[Catch2](https://github.com/catchorg/Catch2)** v3 — for tests. Fetched
  automatically by CMake when tests are built.

## References

If you use this code base in academic work, please cite one or both of the
following papers, which introduce vertex collocation profiles and the
algorithms implemented here:

1. Lichtenwalter, R. N., & Chawla, N. V. (2012). *Vertex collocation profiles:
   subgraph counting for link analysis and prediction.* In Proceedings of the
   21st International Conference on World Wide Web (WWW '12), pp. 1019–1028.
   [doi:10.1145/2187836.2187973](https://doi.org/10.1145/2187836.2187973)

2. Lichtenwalter, R. N., & Chawla, N. V. (2014). *Vertex collocation profiles:
   theory, computation, and results.* SpringerPlus, 3:116. Open access
   (CC BY 2.0).
   [doi:10.1186/2193-1801-3-116](https://doi.org/10.1186/2193-1801-3-116)
   · [PMC open-access copy](https://pmc.ncbi.nlm.nih.gov/articles/PMC4212056/)

## License

BSD 3-Clause — see [LICENSE](LICENSE). Versions prior to 2.0.0 were
released under GPL-3.0-or-later; downstream users who received those
earlier releases retain their rights under that license.

# vcp regression harness

A golden-master characterization harness for the `vcp` CLI tools. Compares the
outputs of two builds — a "legacy" build (pre-modernization) and a "current"
build — on a common set of fixture graphs and paper-derived known-answer
scenarios.

The harness answers two questions:

1. **Correctness.** Does each build produce outputs consistent with the
   published paper (Lichtenwalter & Chawla, "Vertex collocation profiles:
   theory, computation, and results", SpringerPlus 2014, 3:116)?
2. **Regression.** Do the two builds produce identical outputs on the same
   inputs?

## Structure

```text
regression/
├── README.md                  # this file
├── build-legacy.sh            # build legacy binaries from a git revision
├── build-current.sh           # build current (CMake) binaries from a git revision
├── generate_fixtures.py       # deterministic fixture generator
├── run.sh                     # main comparison driver
├── fixtures/                  # committed generated graphs + pair files
│   ├── graphs_undirected/
│   ├── graphs_directed/
│   ├── graphs_multirelational_r2/
│   └── pairs/
├── expected/
│   └── cardinalities.txt      # paper Table 2 ground truth
├── bin/                       # (gitignored) built binaries
│   ├── legacy/
│   └── current/
├── worktree-legacy/           # (gitignored) legacy worktree
├── worktree-current/          # (gitignored) current worktree
└── results/                   # (gitignored) regression report + per-scenario outputs
```

## Usage

Build both versions, then run the harness:

```bash
./build-legacy.sh  ba1592d                          # pre-modernization commit
./build-current.sh feature/tooling-modernization    # any target revision
./run.sh
```

A markdown report is written to `results/report.md`. Per-scenario outputs land
under `results/paper/` and `results/characterization/` so any mismatches can be
inspected directly.

To re-target a future revision (say, a post-C++20 refactor), rebuild:

```bash
./build-current.sh feature/cpp20-upgrade
./run.sh
```

## Design notes

- The legacy build strips `-Werror` from the 2013-era Makefile. GCC 14 emits
  warnings (e.g. `-Walloc-size`) that GCC 4.8 didn't, which would otherwise
  abort the build. The warnings themselves are legitimate — they're exactly what
  the modernization pass fixed. Here we're building the legacy code to
  characterize its output, not to hold it to modern clean-compile standards.
- Both builds happen in separate `git worktree` checkouts so the main working
  tree HEAD is never moved. Worktrees are cleaned up on re-runs.
- Fixtures are committed (generated deterministically from seeded PRNG) so
  regressions are reproducible across machines without re-running the Python
  generator.
- Version-string and help-text output are deliberately excluded from the
  comparison — TCLAP (legacy) and CLI11 (current) format `--help` differently,
  and the version string changed from `1.0.0` to `2.0.0`. Neither represents a
  real behavior change.

## Known findings (2026-04-22 run)

These were uncovered by the harness on the 2026-04-22 run against `ba1592d`
(legacy) and `feature/tooling-modernization` HEAD (current). Source fixes are
explicitly out of scope for this harness — findings are flagged here for
follow-up work.

### Preexisting bugs present in BOTH versions

1. **`directed_to_undirected` infinite loop** on any directed graph whose
   in-degree and out-degree sets are not identical. Lines 56-64 of
   `src/directed_to_undirected.cpp` have two `while (it != end)` tail loops that
   push to `neighbors` without incrementing the iterator. OOM-killed on 7 of 10
   directed fixtures.
2. **`directed_to_undirected` segfault** on graphs with no edges
   (`disolated_5`). After the merge loops complete with an empty `neighbors`
   vector, `neighbors.end() - 1` is dereferenced — undefined behavior below
   `begin()`.
3. **`vcp_generate` on k4 directed with (n=4, r=1, d=1)** emits
   `18446744073709551615` (`SIZE_MAX`) as one element in the legacy output and
   segfaults in the current build — both are manifestations of what appears to
   be an uninitialized-memory read. Current's behavior is louder; legacy's is
   silent but wrong.

### Fixed incidentally by the modernization

- **`vcp_generate` on multirelational (r=2) path graphs** previously emitted
  huge pointer-like element addresses (e.g.
  `12948048434194688,0 12948048434194716,1`). The modernization's
  default-member-init / enum-underlying-type changes reset these to the expected
  clean values (`0,0 28,1`). The relative offsets within each VCP map were
  always correct; only the base was corrupt.

### Not yet investigated

- Performance: the harness diffs outputs, not execution time. A follow-up could
  add `time` wrappers and flag >10% regressions on the random graphs.

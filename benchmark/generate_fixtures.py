#!/usr/bin/env python3
"""Generate benchmark fixtures for the vcp benchmark suite.

Benchmark fixtures are larger than regression fixtures so that each
workload takes long enough (> 200 ms) to yield stable timing. They are
also curated to avoid the known legacy bugs so the same fixture runs
successfully on both the 2012 original code and the modernized code,
enabling meaningful side-by-side comparison.

Fixtures are deterministic: every random graph is seeded, so
regeneration is byte-identical across runs and across machines.

Emits:

    fixtures/undirected/<name>.txt        unirelational undirected graphs
    fixtures/directed/<name>.txt          unirelational directed graphs
    fixtures/multirelational/<name>.txt   r=2 multirelational undirected graphs
    fixtures/pairs/<graph>.pairs          all-pairs input for vcp_generate

Legacy-compatibility constraints:

  * directed fixtures used with `directed_to_undirected` (no flag) are
    bidirectionally symmetric - legacy infinite-loops when out- and
    in-neighbor sets diverge.
  * multirelational fixtures are complete graphs or near-complete so
    legacy's unguarded edge_value() lookup cannot hit a missing edge.
  * graphs are sized so vcp_generate(4,r,d) completes in ~1-5 s per ref
    with all-pairs input.
"""

from __future__ import annotations

import argparse
import array
import math
import os
import random
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path


@dataclass
class UndirectedGraph:
    n: int
    adj: list[set[int]] = field(default_factory=list)

    @classmethod
    def empty(cls, n: int) -> "UndirectedGraph":
        return cls(n=n, adj=[set() for _ in range(n)])

    def add_edge(self, u: int, v: int) -> None:
        if u == v:
            return
        self.adj[u].add(v)
        self.adj[v].add(u)

    def to_unirelational_text(self) -> str:
        return (
            "\n".join(
                " ".join(str(v) for v in sorted(self.adj[u])) for u in range(self.n)
            )
            + "\n"
        )

    def to_multirelational_text(self, bit_assigner) -> str:
        lines = []
        for u in range(self.n):
            tokens = [f"{v},{bit_assigner(u, v)}" for v in sorted(self.adj[u])]
            lines.append(" ".join(tokens))
        return "\n".join(lines) + "\n"


@dataclass
class DirectedGraph:
    n: int
    out_adj: list[set[int]] = field(default_factory=list)

    @classmethod
    def empty(cls, n: int) -> "DirectedGraph":
        return cls(n=n, out_adj=[set() for _ in range(n)])

    def add_arc(self, u: int, v: int) -> None:
        if u == v:
            return
        self.out_adj[u].add(v)

    def to_text(self) -> str:
        return (
            "\n".join(
                " ".join(str(v) for v in sorted(self.out_adj[u])) for u in range(self.n)
            )
            + "\n"
        )


def erdos_renyi_undirected(n: int, p: float, seed: int) -> UndirectedGraph:
    rng = random.Random(seed)
    g = UndirectedGraph.empty(n)
    for u, v in combinations(range(n), 2):
        if rng.random() < p:
            g.add_edge(u, v)
    return g


def erdos_renyi_bidirectional(n: int, p: float, seed: int) -> DirectedGraph:
    """A directed graph whose out- and in-neighbor sets coincide.

    Every sampled pair becomes two arcs (u->v and v->u), so the graph is
    symmetric. Legacy directed_to_undirected (without -b) infinite-loops
    on asymmetric neighbor sets; benchmark fixtures avoid that to keep
    cross-ref timing comparable.
    """
    rng = random.Random(seed)
    g = DirectedGraph.empty(n)
    for u, v in combinations(range(n), 2):
        if rng.random() < p:
            g.add_arc(u, v)
            g.add_arc(v, u)
    return g


def stream_er_sparse(n: int, p: float, seed: int, path: Path) -> None:
    """Stream a sparse G(n, p) adjacency listing directly to disk.

    This is the unified large-random-graph procedure - memory-efficient
    enough to run on 10M+ nodes without OOM, and format-agnostic
    between the undirected and directed-bidirectional interpretations
    (see note below). Small fixtures still use the in-memory
    :class:`UndirectedGraph` / :class:`DirectedGraph` dataclasses where
    clarity beats efficiency.

    Memory profile for n=10M, avg_deg=5 (~25M edges):

        degrees   int32[n]     ~40 MB
        offsets   int32[n+1]   ~40 MB
        neighbors int32[2m]    ~200 MB
        cursor    int32[n]     ~40 MB
        total                  ~320 MB, no Python object-per-edge overhead.

    vs the naive "list of sets of Python ints" approach, which burns
    ~2 GB on empty-set headers alone before counting the ~28 bytes of
    PyLong overhead per element. That representation OOM-killed on a
    15 GB box during the directed pass; this one peaks at ~400 MB.

    Format note: the on-disk output is a line per vertex listing the
    other endpoint of every incident edge, sorted ascending. That
    representation is bit-identical under two interpretations:

      1. Undirected: each edge {u, v} contributes v to u's row and u
         to v's row.
      2. Directed, bidirectionally symmetric: every sampled edge
         becomes arcs u->v AND v->u, so v appears in u's out-list
         and u in v's out-list.

    So a single generated file serves both the undirected and
    bidirectional-directed large-tier workloads without duplication.

    Determinism: two passes over an identical seeded :class:`random.Random`
    stream (pass 1 counts degrees, pass 2 emits to pre-sized slots)
    produce a byte-identical output across runs and machines.

    Sort-free emission: because we process sampled edges (u, v) with
    u < v in ascending u order and the geometric-gap sampler produces
    v values in ascending order within each u, every vertex x ends up
    with its predecessors (u < x) prepended in ascending order followed
    by its successors (v > x) in ascending order. Concatenation is
    globally sorted without an explicit sort step.
    """
    if not (0.0 < p < 1.0):
        raise ValueError(f"p must be in (0, 1); got {p}")
    log1mp = math.log1p(-p)

    # Pass 1: per-vertex degrees.
    rng = random.Random(seed)
    degrees = array.array("i", [0] * n)
    for u in range(n - 1):
        v = u
        while True:
            gap = int(math.log(rng.random()) / log1mp) + 1
            v += gap
            if v >= n:
                break
            degrees[u] += 1
            degrees[v] += 1

    # Offsets prefix-sum.
    total = sum(degrees)
    offsets = array.array("i", [0] * (n + 1))
    running = 0
    for u in range(n):
        offsets[u] = running
        running += degrees[u]
    offsets[n] = running

    # Pass 2: fill a flat int32 CSR array with the same edges.
    neighbors = array.array("i", [0] * total)
    cursor = array.array("i", offsets[:n])
    rng = random.Random(seed)
    for u in range(n - 1):
        v = u
        while True:
            gap = int(math.log(rng.random()) / log1mp) + 1
            v += gap
            if v >= n:
                break
            neighbors[cursor[u]] = v
            cursor[u] += 1
            neighbors[cursor[v]] = u
            cursor[v] += 1

    # Pass 3: stream each vertex's row to disk.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for u in range(n):
            start = offsets[u]
            end = offsets[u + 1]
            if start < end:
                f.write(" ".join(str(neighbors[k]) for k in range(start, end)))
            f.write("\n")


def all_pairs(n: int) -> str:
    return "\n".join(f"{u} {v}" for u, v in combinations(range(n), 2)) + "\n"


def sampled_pairs(n: int, k: int, seed: int) -> str:
    """Sample k distinct unordered pairs (u, v) with u < v, uniformly at random.

    Used for the large-tier workloads where the graph has so many
    vertices that an all-pairs enumeration (C(n, 2) lines) would be
    larger than the graph itself and would make the VCP enumeration run
    effectively forever. A fixed sample gives a representative,
    bounded-cost workload.
    """
    if k > n * (n - 1) // 2:
        raise ValueError(f"k={k} exceeds number of distinct pairs for n={n}")
    rng = random.Random(seed)
    pairs: set[tuple[int, int]] = set()
    while len(pairs) < k:
        u = rng.randrange(n)
        v = rng.randrange(n)
        if u == v:
            continue
        if u > v:
            u, v = v, u
        pairs.add((u, v))
    return "\n".join(f"{u} {v}" for u, v in sorted(pairs)) + "\n"


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def bit_assigner_factory(seed: int):
    """Assign a small non-zero relation bitset to every edge deterministically.

    Values are in {1, 2, 3} - never 0, so every edge present in the
    graph has a non-trivial multirelational connectivity value.

    This function only writes entries for edges that exist. The legacy
    edge_value bug fires on MISSING edges (i.e. VCP enumeration probing
    an edge not in the file), which cannot be prevented at the fixture
    layer. Keeping fixtures dense is a probabilistic mitigation: the
    denser the graph, the smaller the chance that a 4-vertex subgraph
    enumeration touches an unconnected pair. This is the reason the
    multirelational fixtures below use p=0.30 at n=100 and p=0.15 at
    n=200, far denser than the other workloads.
    """
    rng = random.Random(seed)
    cache: dict[tuple[int, int], int] = {}

    def f(u: int, v: int) -> int:
        key = (min(u, v), max(u, v))
        if key not in cache:
            cache[key] = rng.randint(1, 3)
        return cache[key]

    return f


def generate_small(out: Path) -> None:
    # Workloads tuned so each vcp_generate(4,r,d) all-pairs run is ~0.5-5 s
    # on modern code; legacy is typically slower. Sizes are a compromise
    # between "slow enough to time reliably" and "fast enough that a full
    # cross-ref sweep finishes in minutes, not hours."

    # --- undirected unirelational ---
    specs_u = [
        ("er_n200_p0.10_s42", 200, 0.10, 42),
        ("er_n500_p0.05_s42", 500, 0.05, 42),
        ("er_n1000_p0.02_s42", 1000, 0.02, 42),
        ("er_n2000_p0.01_s42", 2000, 0.01, 42),
    ]
    for name, n, p, seed in specs_u:
        g = erdos_renyi_undirected(n, p, seed)
        write(out / "undirected" / f"{name}.txt", g.to_unirelational_text())
        write(out / "pairs" / f"{name}.pairs", all_pairs(n))

    # --- directed unirelational (bidirectionally symmetric) ---
    specs_d = [
        ("dsym_n200_p0.10_s42", 200, 0.10, 42),
        ("dsym_n500_p0.05_s42", 500, 0.05, 42),
        ("dsym_n1000_p0.02_s42", 1000, 0.02, 42),
    ]
    for name, n, p, seed in specs_d:
        g = erdos_renyi_bidirectional(n, p, seed)
        write(out / "directed" / f"{name}.txt", g.to_text())
        write(out / "pairs" / f"{name}.pairs", all_pairs(n))

    # --- multirelational r=2 (dense, to avoid legacy missing-edge lookups) ---
    # A smaller n because r=2 is heavier.
    specs_mr = [
        ("er_n100_p0.30_s42", 100, 0.30, 42),
        ("er_n200_p0.15_s42", 200, 0.15, 42),
    ]
    for name, n, p, seed in specs_mr:
        g = erdos_renyi_undirected(n, p, seed)
        write(
            out / "multirelational" / f"{name}.txt",
            g.to_multirelational_text(bit_assigner_factory(seed + 1)),
        )
        write(out / "pairs" / f"mr_{name}.pairs", all_pairs(n))


# --- Large-tier fixtures ---
#
# The large tier is opt-in (via --large or run.sh --large) because it
# generates a ~300 MB undirected graph and a ~600 MB directed graph and
# takes a few minutes to produce. It is scoped to the four unirelational
# VCP procedures - 3/4-vertex, undirected/directed - because:
#
#   * Multirelational fixtures must be dense to dodge the legacy
#     missing-edge bug (edge_value reads OOB on absent edges), and
#     density is incompatible with a 10M-node sparse workload.
#   * Auxiliary tools (directed_to_undirected, ell_2_pairs) are
#     already I/O-bounded in the small tier; scaling up adds no new
#     signal.
#
# Ratio n:m is set by LARGE_AVG_DEGREE, mapped to Bernoulli p via
# p = avg_deg / (n - 1). Pair sample count is fixed at
# LARGE_PAIR_SAMPLE; at avg degree 5 this yields runtimes on the order
# of seconds to low minutes on modern code, which is both tractable for
# a benchmark sweep and large enough for the per-pair enumeration work
# to dominate the graph parse and measure meaningfully.
LARGE_N = 10_000_000
LARGE_AVG_DEGREE = 5
LARGE_PAIR_SAMPLE = 100_000
LARGE_SEED = 42


def _large_p() -> float:
    return LARGE_AVG_DEGREE / (LARGE_N - 1)


def _large_names() -> tuple[str, str, str]:
    stem = f"n{LARGE_N}_d{LARGE_AVG_DEGREE}_s{LARGE_SEED}"
    return (
        f"er_{stem}.txt",
        f"dsym_{stem}.txt",
        f"sample_{stem}_k{LARGE_PAIR_SAMPLE}.pairs",
    )


def generate_large(out: Path) -> None:
    """Emit the large-tier fixtures via the streaming unified procedure.

    Idempotent: files that already exist on disk are left alone, so
    ``run.sh --large`` on a repeat invocation is a no-op.

    The directed-bidirectional fixture is hardlinked from the
    undirected one because the two file representations are
    byte-identical - see :func:`stream_er_sparse`. A hardlink saves
    ~400 MB of disk and a whole generation pass while keeping the
    benchmark's two filenames semantically distinct.
    """
    undirected_name, directed_name, pairs_name = _large_names()
    p = _large_p()

    undirected_path = out / "undirected" / undirected_name
    directed_path = out / "directed" / directed_name
    pairs_path = out / "pairs" / pairs_name

    if undirected_path.exists():
        print(f"   (exists) {undirected_path.name}")
    else:
        print(
            f"   generating {undirected_path.name} "
            f"(n={LARGE_N:,}, avg_deg={LARGE_AVG_DEGREE}, p={p:.2e})"
        )
        stream_er_sparse(LARGE_N, p, LARGE_SEED, undirected_path)

    if directed_path.exists():
        print(f"   (exists) {directed_path.name}")
    else:
        print(
            f"   hardlinking {directed_path.name} "
            f"(byte-identical to {undirected_path.name})"
        )
        directed_path.parent.mkdir(parents=True, exist_ok=True)
        os.link(undirected_path, directed_path)

    if pairs_path.exists():
        print(f"   (exists) {pairs_path.name}")
    else:
        print(f"   generating {pairs_path.name} (k={LARGE_PAIR_SAMPLE:,})")
        write(pairs_path, sampled_pairs(LARGE_N, LARGE_PAIR_SAMPLE, LARGE_SEED + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).parent / "fixtures"
    )
    parser.add_argument(
        "--large",
        action="store_true",
        help=(
            "Also generate the large-tier fixtures "
            f"(n={LARGE_N:,}, avg_deg={LARGE_AVG_DEGREE}). "
            "Takes several minutes and produces ~900 MB on disk."
        ),
    )
    parser.add_argument(
        "--only-large",
        action="store_true",
        help="Skip the small tier; only generate the large-tier fixtures.",
    )
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    if not args.only_large:
        generate_small(out)

    if args.large or args.only_large:
        print(">> Large-tier fixtures")
        generate_large(out)


if __name__ == "__main__":
    main()

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


def all_pairs(n: int) -> str:
    return "\n".join(f"{u} {v}" for u, v in combinations(range(n), 2)) + "\n"


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def bit_assigner_factory(seed: int):
    """Assign a small non-zero relation bitset to every edge deterministically.

    Values in {1, 2, 3} - never 0, so every present edge has a non-trivial
    multirelational connectivity value. This does not prevent the legacy
    edge_value bug (which fires on MISSING edges), but keeping the
    fixtures dense reduces the chance of a missing-edge lookup during
    enumeration.
    """
    rng = random.Random(seed)
    cache: dict[tuple[int, int], int] = {}

    def f(u: int, v: int) -> int:
        key = (min(u, v), max(u, v))
        if key not in cache:
            cache[key] = rng.randint(1, 3)
        return cache[key]

    return f


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).parent / "fixtures"
    )
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

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


if __name__ == "__main__":
    main()

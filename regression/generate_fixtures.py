#!/usr/bin/env python3
"""Generate vcp regression-test fixture files deterministically.

The vcp CLI tools consume graphs in three closely related text formats
(one line per vertex, vertex IDs implicit by line number):

- Undirected unirelational (`vcp::graph`): space-separated neighbor IDs.
- Directed unirelational (`vcp::directed_graph`): space-separated OUT-neighbor
  IDs; in-neighbors are derived by the reader.
- Multirelational (`vcp::multirelational_graph<r>`): space-separated
  `neighbor_id,relation_bitset` pairs per line.

This script writes a curated set of structured and random graphs plus the
vertex-pair input lists for `vcp_generate`. All random graphs use a fixed
PRNG seed so regeneration is byte-identical across runs.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


@dataclass
class UndirectedGraph:
    """Undirected simple graph stored as a neighbor-set per vertex."""

    n: int
    adj: list[set[int]]

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
        """Emit the multirelational format.

        `bit_assigner(u, v)` returns a small positive integer encoding the
        relation bitset on the edge {u, v}. The same bitset is emitted on
        both endpoints (the underlying storage duplicates each edge).
        """
        lines = []
        for u in range(self.n):
            tokens = []
            for v in sorted(self.adj[u]):
                bits = bit_assigner(u, v)
                tokens.append(f"{v},{bits}")
            lines.append(" ".join(tokens))
        return "\n".join(lines) + "\n"


@dataclass
class DirectedGraph:
    """Directed graph stored as an out-neighbor-set per vertex."""

    n: int
    out_adj: list[set[int]]

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


def complete_graph(n: int) -> UndirectedGraph:
    g = UndirectedGraph.empty(n)
    for u, v in combinations(range(n), 2):
        g.add_edge(u, v)
    return g


def path_graph(n: int) -> UndirectedGraph:
    g = UndirectedGraph.empty(n)
    for i in range(n - 1):
        g.add_edge(i, i + 1)
    return g


def star_graph(n: int) -> UndirectedGraph:
    """Star with `n` total vertices: vertex 0 is the hub."""
    g = UndirectedGraph.empty(n)
    for i in range(1, n):
        g.add_edge(0, i)
    return g


def cycle_graph(n: int) -> UndirectedGraph:
    g = UndirectedGraph.empty(n)
    for i in range(n):
        g.add_edge(i, (i + 1) % n)
    return g


def disconnected_triangles(k: int) -> UndirectedGraph:
    """k disjoint triangles, 3k vertices total."""
    g = UndirectedGraph.empty(3 * k)
    for i in range(k):
        base = 3 * i
        g.add_edge(base, base + 1)
        g.add_edge(base + 1, base + 2)
        g.add_edge(base, base + 2)
    return g


def isolated_vertices(n: int) -> UndirectedGraph:
    return UndirectedGraph.empty(n)


def erdos_renyi(n: int, p: float, seed: int) -> UndirectedGraph:
    rng = random.Random(seed)
    g = UndirectedGraph.empty(n)
    for u, v in combinations(range(n), 2):
        if rng.random() < p:
            g.add_edge(u, v)
    return g


def directed_from(g: UndirectedGraph) -> DirectedGraph:
    """Random-orientation directed graph: each undirected edge becomes a single arc."""
    d = DirectedGraph.empty(g.n)
    for u in range(g.n):
        for v in g.adj[u]:
            if u < v:
                d.add_arc(u, v)  # always u -> v (smaller to larger)
    return d


def directed_bidirectional(g: UndirectedGraph) -> DirectedGraph:
    """Each undirected edge becomes two arcs (bidirectional)."""
    d = DirectedGraph.empty(g.n)
    for u in range(g.n):
        for v in g.adj[u]:
            d.add_arc(u, v)
    return d


def random_directed(n: int, p: float, seed: int) -> DirectedGraph:
    rng = random.Random(seed)
    d = DirectedGraph.empty(n)
    for u in range(n):
        for v in range(n):
            if u != v and rng.random() < p:
                d.add_arc(u, v)
    return d


def all_pairs(n: int) -> list[tuple[int, int]]:
    return [(u, v) for u, v in combinations(range(n), 2)]


# --- Fixture catalog ---

STRUCTURED_UNDIRECTED: dict[str, UndirectedGraph] = {
    "single_vertex": isolated_vertices(1),
    "two_vertices_no_edge": isolated_vertices(2),
    "two_vertices_edge": UndirectedGraph(n=2, adj=[{1}, {0}]),
    "p3_path": path_graph(3),
    "p4_path": path_graph(4),
    "k3_triangle": complete_graph(3),
    "k4_complete": complete_graph(4),
    "k5_complete": complete_graph(5),
    "k6_complete": complete_graph(6),
    "c4_cycle": cycle_graph(4),
    "c5_cycle": cycle_graph(5),
    "c6_cycle": cycle_graph(6),
    "s4_star": star_graph(5),  # hub + 4 leaves
    "s5_star": star_graph(6),
    "two_triangles_disconnected": disconnected_triangles(2),
    "three_triangles_disconnected": disconnected_triangles(3),
    "isolated_5": isolated_vertices(5),
}


def _edge_only_4() -> UndirectedGraph:
    """V=4, single edge 0-1, vertices 2 and 3 isolated."""
    g = UndirectedGraph.empty(4)
    g.add_edge(0, 1)
    return g


def _triangle_plus_isolate_4() -> UndirectedGraph:
    """V=4, triangle on {0, 1, 2}, vertex 3 isolated."""
    g = UndirectedGraph.empty(4)
    for u, v in [(0, 1), (1, 2), (0, 2)]:
        g.add_edge(u, v)
    return g


def _k22_4() -> UndirectedGraph:
    """V=4, complete bipartite K_{2,2}: parts {0,1} and {2,3}."""
    g = UndirectedGraph.empty(4)
    for u in (0, 1):
        for v in (2, 3):
            g.add_edge(u, v)
    return g


def _directed_asym_from(g: UndirectedGraph, orientation) -> DirectedGraph:
    """Directed orientation where each undirected edge is mapped to one arc
    determined by `orientation(u, v)` returning (src, dst)."""
    d = DirectedGraph.empty(g.n)
    for u in range(g.n):
        for v in g.adj[u]:
            if u < v:
                src, dst = orientation(u, v)
                d.add_arc(src, dst)
    return d


STRUCTURED_DIRECTED: dict[str, DirectedGraph] = {
    "dk3_oriented": directed_from(complete_graph(3)),
    "dk3_bidirectional": directed_bidirectional(complete_graph(3)),
    "dk4_oriented": directed_from(complete_graph(4)),
    "dk4_bidirectional": directed_bidirectional(complete_graph(4)),
    "dcycle4": directed_from(cycle_graph(4)),
    "dp4_forward": directed_from(path_graph(4)),
    "disolated_5": DirectedGraph.empty(5),
    # --- Bug 2 coverage fixtures (VCP_{4,1,1}) ---
    # Each graph isolates a specific combination of v1v2 connectivity and
    # Defect-C trigger (shared v3 with nonzero c13 AND c23). See
    # CHANGELOG.md / vcp_4_1_1.hpp commit message for the full taxonomy.
    # Mutual edge only: Defect-A/B trigger, no shared v3.
    "dbug2_mutual_edge": directed_bidirectional(_edge_only_4()),
    # Asymmetric edge only: no Defect triggers; control baseline.
    "dbug2_asym_edge": directed_from(_edge_only_4()),
    # Mutual triangle + isolate: Defect-A/B + Defect-C on (0,1), (0,2), (1,2).
    "dbug2_mutual_triangle": directed_bidirectional(_triangle_plus_isolate_4()),
    # Transitive tournament on {0,1,2} (arcs 0→1, 0→2, 1→2) + vertex 3 isolated.
    # Asymmetric shared-v3 exercise of Defect-C without triggering Defect-A/B
    # (v1v2 is OUT/IN, never BOTH).
    "dbug2_asym_triangle": _directed_asym_from(
        _triangle_plus_isolate_4(),
        orientation=lambda u, v: (u, v),  # smaller → larger (transitive)
    ),
    # K_{2,2} mutual: Defect-C via mutual shared-v3 without Defect-A/B
    # on within-part pairs (v1v2=0).
    "dbug2_k22_mutual": directed_bidirectional(_k22_4()),
    # K_{2,2} asymmetric (left-to-right): Defect-C via asymmetric shared-v3.
    "dbug2_k22_asym": _directed_asym_from(
        _k22_4(),
        orientation=lambda u, v: (u, v),
    ),
}

# Random graphs: fixed seed so regeneration is byte-stable.
RANDOM_UNDIRECTED: dict[str, UndirectedGraph] = {
    "er_n10_p0.3_s42": erdos_renyi(10, 0.3, 42),
    "er_n10_p0.5_s42": erdos_renyi(10, 0.5, 42),
    "er_n20_p0.2_s42": erdos_renyi(20, 0.2, 42),
    "er_n50_p0.1_s42": erdos_renyi(50, 0.1, 42),
    "er_n100_p0.05_s42": erdos_renyi(100, 0.05, 42),
}

RANDOM_DIRECTED: dict[str, DirectedGraph] = {
    "rd_n10_p0.3_s42": random_directed(10, 0.3, 42),
    "rd_n20_p0.15_s42": random_directed(20, 0.15, 42),
    "rd_n50_p0.05_s42": random_directed(50, 0.05, 42),
}


def varying_bit_assigner_r2(seed: int):
    """Deterministic per-edge bit pattern drawing from {1, 2, 3} for r=2."""
    rng = random.Random(seed)
    cache: dict[frozenset[int], int] = {}

    def f(u: int, v: int) -> int:
        key = frozenset({u, v})
        if key not in cache:
            cache[key] = rng.choice([1, 2, 3])
        return cache[key]

    return f


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=Path(__file__).parent / "fixtures",
        type=Path,
        help="Fixture output root (default: regression/fixtures)",
    )
    args = parser.parse_args()
    out: Path = args.output

    # Undirected unirelational
    for name, g in {**STRUCTURED_UNDIRECTED, **RANDOM_UNDIRECTED}.items():
        write(out / "graphs_undirected" / f"{name}.txt", g.to_unirelational_text())

    # Directed unirelational
    for name, d in {**STRUCTURED_DIRECTED, **RANDOM_DIRECTED}.items():
        write(out / "graphs_directed" / f"{name}.txt", d.to_text())

    # Multirelational r=2 — rebuild structured undirected graphs with
    # deterministic per-edge bit patterns in {1, 2, 3}.
    for name, g in STRUCTURED_UNDIRECTED.items():
        if g.n < 2 or all(not nbrs for nbrs in g.adj):
            # Skip trivial graphs (empty / single vertex) — no edges to annotate.
            continue
        assigner = varying_bit_assigner_r2(seed=hash(name) & 0xFFFFFFFF)
        write(
            out / "graphs_multirelational_r2" / f"{name}.txt",
            g.to_multirelational_text(assigner),
        )
    # Also a couple of random graphs with r=2.
    for name, g in RANDOM_UNDIRECTED.items():
        if g.n > 50:
            continue  # skip the largest; multirelational is O(r) heavier
        assigner = varying_bit_assigner_r2(seed=hash(name) & 0xFFFFFFFF)
        write(
            out / "graphs_multirelational_r2" / f"{name}.txt",
            g.to_multirelational_text(assigner),
        )

    # Pair files for vcp_generate. All-pairs for the small graphs, truncated
    # for the larger ones so test runtime stays bounded.
    pair_limits = {
        "single_vertex": 0,
        "two_vertices_no_edge": 1,
        "two_vertices_edge": 1,
        "p3_path": 3,
        "p4_path": 6,
        "k3_triangle": 3,
        "k4_complete": 6,
        "k5_complete": 10,
        "k6_complete": 15,
        "c4_cycle": 6,
        "c5_cycle": 10,
        "c6_cycle": 15,
        "s4_star": 10,
        "s5_star": 15,
        "two_triangles_disconnected": 15,
        "three_triangles_disconnected": 36,
        "isolated_5": 10,
        "er_n10_p0.3_s42": 45,
        "er_n10_p0.5_s42": 45,
        "er_n20_p0.2_s42": 50,
        "er_n50_p0.1_s42": 100,
        "er_n100_p0.05_s42": 200,
        # Bug-2 fixtures: all 6 pairs of V=4 graph.
        "dbug2_mutual_edge": 6,
        "dbug2_asym_edge": 6,
        "dbug2_mutual_triangle": 6,
        "dbug2_asym_triangle": 6,
        "dbug2_k22_mutual": 6,
        "dbug2_k22_asym": 6,
    }

    for name, limit in pair_limits.items():
        if limit == 0:
            continue
        fixture = (
            STRUCTURED_UNDIRECTED.get(name)
            or RANDOM_UNDIRECTED.get(name)
            or STRUCTURED_DIRECTED.get(name)
            or RANDOM_DIRECTED.get(name)
        )
        if fixture is None:
            raise KeyError(f"pair_limits references unknown fixture: {name}")
        pairs = all_pairs(fixture.n)[:limit]
        content = "\n".join(f"{u} {v}" for u, v in pairs) + "\n"
        write(out / "pairs" / f"{name}.txt", content)

    print(f"Fixtures written under {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

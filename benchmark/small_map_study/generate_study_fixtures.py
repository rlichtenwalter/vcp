#!/usr/bin/env python3
"""Generate multirelational fixtures for the small-map k-probe study.

Fixtures fall into three families:

  1. r=2 reuses the existing benchmark fixtures (no generation needed).
     These exercise the generic r>=2 code paths with a maximally-saturated
     key space (observed k saturates at 4 for edge_types).

  2. r=30 fixtures span more of the edge-value range (values uniformly
     sampled from [1, 2**30 - 1]) so that edge_types cardinality is
     driven by data diversity, not by the r=2 cap. The study needs
     these to separate "many distinct values" from "few, repeated".

  3. A hub fixture deliberately produces v3_count spikes: one high-degree
     vertex connected to every other, plus a sparse ER scaffold. When
     the hub participates as v1 or v2, v3_count can approach n, and
     temp_edge_types / counts grow along with it.

All fixtures are deterministic (seeded). Directed files are hardlinks
of the undirected files because the on-disk representation of a
bidirectionally-symmetric multirelational graph is byte-identical under
both interpretations (each pair (u, v, value) appears as u,value on
row v and as v,value on row u).

Emits to benchmark/fixtures/multirelational/ and benchmark/fixtures/pairs/.
"""

from __future__ import annotations

import argparse
import os
import random
from itertools import combinations
from pathlib import Path


def all_pairs(n: int) -> str:
    return "\n".join(f"{u} {v}" for u, v in combinations(range(n), 2)) + "\n"


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def er_multirelational(
    n: int, p: float, seed: int, value_range: tuple[int, int]
) -> str:
    """Erdős–Rényi multirelational graph, edge values sampled in value_range.

    value_range is closed [lo, hi]. The reader infers r from edge values;
    nothing is emitted to size r explicitly.
    """
    rng = random.Random(seed)
    val_rng = random.Random(seed + 1)
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    for u, v in combinations(range(n), 2):
        if rng.random() < p:
            val = val_rng.randint(*value_range)
            adj[u].append((v, val))
            adj[v].append((u, val))
    lines = []
    for u in range(n):
        tokens = [f"{v},{val}" for v, val in sorted(adj[u])]
        lines.append(" ".join(tokens))
    return "\n".join(lines) + "\n"


def hub_graph(
    n: int, hub: int, p_er: float, seed: int, value_range: tuple[int, int]
) -> str:
    """Graph with one hub vertex connected to everyone, plus ER scaffold.

    Stresses the v3_count expansion path. Every edge gets a random value
    from value_range.
    """
    rng = random.Random(seed)
    val_rng = random.Random(seed + 1)
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    # Hub edges
    for v in range(n):
        if v == hub:
            continue
        val = val_rng.randint(*value_range)
        adj[hub].append((v, val))
        adj[v].append((hub, val))
    # ER scaffold on non-hub vertices
    for u, v in combinations(range(n), 2):
        if u == hub or v == hub:
            continue
        if rng.random() < p_er:
            val = val_rng.randint(*value_range)
            adj[u].append((v, val))
            adj[v].append((u, val))
    lines = []
    for u in range(n):
        tokens = [f"{v},{val}" for v, val in sorted(adj[u])]
        lines.append(" ".join(tokens))
    return "\n".join(lines) + "\n"


def hardlink_as_directed(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    os.link(src, dst)


def generate(out: Path) -> None:
    # --- r=30, moderate diversity (4 bits used) ---
    # Edge values in [1, 15]: up to 15 distinct values feasible in any
    # temp_edge_types, but capped by the specific graph's diversity.
    specs_r30_diverse4 = [
        ("mr_r30_div4_n100_p0.30_s42", 100, 0.30, 42),
        ("mr_r30_div4_n200_p0.15_s42", 200, 0.15, 42),
    ]
    for name, n, p, seed in specs_r30_diverse4:
        content = er_multirelational(n, p, seed, value_range=(1, 15))
        ugp = out / "multirelational" / f"{name}.txt"
        write(ugp, content)
        write(out / "pairs" / f"mr_{name}.pairs", all_pairs(n))
        hardlink_as_directed(ugp, out / "multirelational" / f"{name}_dsym.txt")

    # --- r=30, high diversity (full 30-bit range) ---
    # Edge values in [1, 2**30 - 1]: extremely high potential diversity,
    # though in practice bounded by the number of distinct edges. This is
    # the stress case for "dense tier can't fit".
    specs_r30_full = [
        ("mr_r30_full_n100_p0.30_s42", 100, 0.30, 42),
        ("mr_r30_full_n200_p0.15_s42", 200, 0.15, 42),
    ]
    for name, n, p, seed in specs_r30_full:
        content = er_multirelational(n, p, seed, value_range=(1, (1 << 30) - 1))
        ugp = out / "multirelational" / f"{name}.txt"
        write(ugp, content)
        write(out / "pairs" / f"mr_{name}.pairs", all_pairs(n))
        hardlink_as_directed(ugp, out / "multirelational" / f"{name}_dsym.txt")

    # --- Hub graphs (v3_count spike path) ---
    hub_specs = [
        ("mr_r30_div4_hub_n200_s42", 200, 200 // 2, 0.05, 42, (1, 15)),
        ("mr_r30_full_hub_n200_s42", 200, 200 // 2, 0.05, 42, (1, (1 << 30) - 1)),
    ]
    for name, n, hub, p_er, seed, value_range in hub_specs:
        content = hub_graph(n, hub, p_er, seed, value_range)
        ugp = out / "multirelational" / f"{name}.txt"
        write(ugp, content)
        write(out / "pairs" / f"mr_{name}.pairs", all_pairs(n))
        hardlink_as_directed(ugp, out / "multirelational" / f"{name}_dsym.txt")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).parent.parent / "fixtures"
    )
    args = parser.parse_args()
    generate(args.output_dir)


if __name__ == "__main__":
    main()

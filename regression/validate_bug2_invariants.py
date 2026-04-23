#!/usr/bin/env python3
"""Validate VCP_{4,1,1} invariants on Bug-2 coverage fixtures.

For each fixture in regression/fixtures/graphs_directed/dbug2_*.txt, run
vcp_generate 4 1 1, parse the output vectors (one per pair), and assert:

  (1) Conservation: sum(counts) == C(V-2, 2) per pair. Catches Defect A
      (unsigned underflow to 2^64-1 → huge sum) and Defect B (off-by-one
      routing → sum = C(V-2, 2) + 1).

  (2) Non-underflow: max(counts) <= C(V-2, 2). Redundant with (1) on V=4
      graphs but explicit as a diagnostic signal.

  (3) Slot localization (V=4 only): for each pair, reconstruct the
      induced 4-vertex subgraph, compute its subgraph address via
      (V12, V13, V14, V23, V24, V34) = (1, 4, 16, 64, 256, 1024), look up
      the VCP element via `vcp_map 4 1 1`, and assert counts[element] == 1.
      Catches Defect C (address arithmetic off by the V1V3 stride → the
      single unit lands in the wrong slot, or out of bounds).

Exit code 0 on all fixtures pass, 1 otherwise.

Usage:
    validate_bug2_invariants.py --bin DIR [--fixtures DIR]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from itertools import permutations
from math import comb
from pathlib import Path

# Directedness values for VCP_{4,1,1} (see include/vcp/vcp_4_1_1.hpp).
OUT, IN, BOTH = 1, 2, 3
# Per-slot magnitudes for the 12-bit subgraph address.
V12, V13, V14, V23, V24, V34 = 1, 4, 16, 64, 256, 1024
EXPECTED_VECTOR_LEN = 2112


def load_element_address_map(vcp_map_bin: Path) -> list[int]:
    """Run `vcp_map 4 1 1` and return the 4096-entry subgraph→element map."""
    result = subprocess.run(
        [str(vcp_map_bin), "4", "1", "1"],
        capture_output=True,
        text=True,
        check=True,
    )
    mapping = [int(line) for line in result.stdout.strip().splitlines()]
    if len(mapping) != 4096:
        raise RuntimeError(f"vcp_map returned {len(mapping)} entries, expected 4096")
    return mapping


def read_directed_graph(path: Path) -> list[list[int]]:
    """Parse directed-graph text format: one line per vertex, space-separated
    out-neighbor IDs. In-neighbors are derived by the reader in the C++ code;
    here we store only out_adj and compute direction on demand."""
    out_adj: list[list[int]] = []
    with path.open() as f:
        for line in f:
            stripped = line.strip()
            out_adj.append([int(x) for x in stripped.split()] if stripped else [])
    return out_adj


def direction_value(out_adj: list[list[int]], u: int, v: int) -> int:
    """Return d(u, v) ∈ {0, OUT, IN, BOTH}."""
    u_to_v = v in out_adj[u]
    v_to_u = u in out_adj[v]
    return (OUT if u_to_v else 0) | (IN if v_to_u else 0)


def subgraph_address(
    out_adj: list[list[int]], v1: int, v2: int, v3: int, v4: int
) -> int:
    """Compute the 12-bit subgraph address for the ordered 4-tuple."""
    return (
        V12 * direction_value(out_adj, v1, v2)
        + V13 * direction_value(out_adj, v1, v3)
        + V14 * direction_value(out_adj, v1, v4)
        + V23 * direction_value(out_adj, v2, v3)
        + V24 * direction_value(out_adj, v2, v4)
        + V34 * direction_value(out_adj, v3, v4)
    )


def expected_element_v4(
    out_adj: list[list[int]], v1: int, v2: int, elem_map: list[int]
) -> int:
    """For V=4: the induced subgraph on {v1, v2, v3, v4} has a unique VCP
    element. v3/v4 are unlabeled, so both orderings must map to the same
    element — we assert this as a sanity check on the map."""
    assert len(out_adj) == 4, "expected_element_v4 requires V=4"
    others = [v for v in range(4) if v != v1 and v != v2]
    assert len(others) == 2
    # Check that both (v3, v4) orderings yield the same element index.
    elements = {
        elem_map[subgraph_address(out_adj, v1, v2, v3, v4)]
        for v3, v4 in permutations(others)
    }
    if len(elements) != 1:
        raise RuntimeError(
            f"pair ({v1},{v2}): v3/v4 orderings disagree on element: {elements}"
        )
    return elements.pop()


def validate_fixture(
    graph_path: Path,
    pair_path: Path,
    vcp_generate_bin: Path,
    elem_map: list[int],
) -> list[str]:
    """Run vcp_generate on one fixture and return a list of error strings
    (empty list if all invariants hold)."""
    out_adj = read_directed_graph(graph_path)
    v_count = len(out_adj)

    with pair_path.open() as f:
        pairs = [tuple(int(x) for x in line.split()) for line in f if line.strip()]

    proc = subprocess.run(
        [str(vcp_generate_bin), "4", "1", "1", str(graph_path)],
        input="\n".join(f"{a} {b}" for a, b in pairs) + "\n",
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return [
            f"vcp_generate crashed (exit={proc.returncode}): "
            f"{proc.stderr.strip()[:200]}"
        ]

    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    if len(lines) != len(pairs):
        return [f"expected {len(pairs)} output lines, got {len(lines)}"]

    expected_sum = comb(v_count - 2, 2)
    errors: list[str] = []
    for (v1, v2), line in zip(pairs, lines, strict=False):
        counts = [int(x) for x in line.split()]
        if len(counts) != EXPECTED_VECTOR_LEN:
            errors.append(
                f"pair ({v1},{v2}): vector length={len(counts)}, "
                f"expected {EXPECTED_VECTOR_LEN}"
            )
            continue
        total = sum(counts)
        peak = max(counts)
        if total != expected_sum:
            errors.append(
                f"pair ({v1},{v2}): sum={total}, expected {expected_sum} "
                f"[Defect A/B indicator]"
            )
        if peak > expected_sum:
            errors.append(
                f"pair ({v1},{v2}): max={peak} exceeds expected {expected_sum} "
                f"[underflow indicator]"
            )
        if v_count == 4:
            expected_idx = expected_element_v4(out_adj, v1, v2, elem_map)
            if counts[expected_idx] != 1:
                errors.append(
                    f"pair ({v1},{v2}): counts[{expected_idx}]="
                    f"{counts[expected_idx]}, expected 1 "
                    f"[Defect C indicator]"
                )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bin",
        type=Path,
        required=True,
        help="Directory containing vcp_generate and vcp_map binaries",
    )
    parser.add_argument(
        "--fixtures",
        type=Path,
        default=Path(__file__).parent / "fixtures",
        help="Fixtures root (default: regression/fixtures)",
    )
    args = parser.parse_args()

    vcp_generate = args.bin / "vcp_generate"
    vcp_map = args.bin / "vcp_map"
    if not vcp_generate.is_file():
        print(f"error: {vcp_generate} not found", file=sys.stderr)
        return 2
    if not vcp_map.is_file():
        print(f"error: {vcp_map} not found", file=sys.stderr)
        return 2

    graphs_dir = args.fixtures / "graphs_directed"
    pairs_dir = args.fixtures / "pairs"
    elem_map = load_element_address_map(vcp_map)

    fixtures = sorted(graphs_dir.glob("dbug2_*.txt"))
    if not fixtures:
        print(f"error: no dbug2_*.txt fixtures under {graphs_dir}", file=sys.stderr)
        return 2

    total_pass = 0
    total_fail = 0
    for graph_file in fixtures:
        name = graph_file.stem
        pair_file = pairs_dir / f"{name}.txt"
        if not pair_file.exists():
            print(f"  SKIP {name} (no pair file)")
            continue
        errors = validate_fixture(graph_file, pair_file, vcp_generate, elem_map)
        if errors:
            total_fail += 1
            print(f"  FAIL {name}")
            for err in errors:
                print(f"    - {err}")
        else:
            total_pass += 1
            print(f"  PASS {name}")

    print()
    print(f"Bug-2 invariant validation: {total_pass} pass, {total_fail} fail")
    return 0 if total_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

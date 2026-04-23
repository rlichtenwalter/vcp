#!/usr/bin/env bash
# Regression harness for vcp.
#
# Compares the outputs of two builds of the vcp CLI tools on a common set
# of fixtures and known-answer scenarios. Answers two questions:
#
#   1. Correctness: does each build produce outputs consistent with the
#      published paper (Lichtenwalter & Chawla 2014, SpringerPlus)?
#   2. Regression: do the two builds produce identical outputs on the
#      same inputs, modulo expected differences (version strings, help
#      text format)?
#
# Usage:
#   ./run.sh [--legacy-bin DIR] [--current-bin DIR] [--results-dir DIR]
#
# Defaults: bin/legacy, bin/current, results/
#
# Exit code: 0 if every comparison passes AND every paper check passes;
# 1 otherwise. Detailed results are written as markdown under results/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LEGACY_BIN="$SCRIPT_DIR/bin/legacy"
CURRENT_BIN="$SCRIPT_DIR/bin/current"
RESULTS_DIR="$SCRIPT_DIR/results"
FIXTURES_DIR="$SCRIPT_DIR/fixtures"
EXPECTED_DIR="$SCRIPT_DIR/expected"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --legacy-bin)   LEGACY_BIN="$2"; shift 2 ;;
        --current-bin)  CURRENT_BIN="$2"; shift 2 ;;
        --results-dir)  RESULTS_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$RESULTS_DIR"/{characterization,paper}
REPORT="$RESULTS_DIR/report.md"

# --- Counters ---
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
declare -a FAIL_NOTES=()

pass() { PASS_COUNT=$((PASS_COUNT + 1)); printf '  \033[32mPASS\033[0m %s\n' "$1"; }
fail() { FAIL_COUNT=$((FAIL_COUNT + 1)); printf '  \033[31mFAIL\033[0m %s\n' "$1"; FAIL_NOTES+=("$1"); }
skip() { SKIP_COUNT=$((SKIP_COUNT + 1)); printf '  \033[33mSKIP\033[0m %s\n' "$1"; }

# --- Sanity: both binary sets must be present ---
for dir in "$LEGACY_BIN" "$CURRENT_BIN"; do
    for tool in vcp_generate vcp_map directed_to_undirected ell_2_pairs; do
        if [ ! -x "$dir/$tool" ]; then
            echo "ERROR: missing binary $dir/$tool" >&2
            exit 2
        fi
    done
done

# --- Report header ---
{
    echo "# vcp regression harness report"
    echo
    echo "- Legacy binaries:  \`$LEGACY_BIN\`"
    echo "- Current binaries: \`$CURRENT_BIN\`"
    echo "- Run timestamp:    $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo
} > "$REPORT"

# =============================================================================
# Phase 1 — Paper correctness checks: vcp_map cardinalities
# =============================================================================
echo "=== Phase 1: paper cardinality checks ==="
{
    echo "## Phase 1 — Paper cardinality checks (vcp_map)"
    echo
    echo "Each row runs \`vcp_map n r d\` and counts distinct element IDs in"
    echo "stdout; that count should equal the cardinality \$|VCP_{n,r,d}|\$"
    echo "from the paper."
    echo
    echo "| Scenario | Expected | Legacy | Current | Result |"
    echo "|---|---|---|---|---|"
} >> "$REPORT"

# Skip very large cases to keep runtime bounded.
declare -A SKIP_LARGE=(
    ["3 5 0"]=1   # 32768 elements — fine
    ["3 4 0"]=1   # 4096 — fine
)

while read -r line; do
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    read -r n r d expected <<< "$line"
    scenario="vcp_map $n $r $d"

    legacy_out="$RESULTS_DIR/paper/legacy_vcp_map_${n}_${r}_${d}.txt"
    current_out="$RESULTS_DIR/paper/current_vcp_map_${n}_${r}_${d}.txt"

    "$LEGACY_BIN/vcp_map" "$n" "$r" "$d" > "$legacy_out" 2>/dev/null
    "$CURRENT_BIN/vcp_map" "$n" "$r" "$d" > "$current_out" 2>/dev/null

    legacy_card=$(sort -u "$legacy_out" | wc -l)
    current_card=$(sort -u "$current_out" | wc -l)

    status="PASS"
    if [ "$legacy_card" -ne "$expected" ] || [ "$current_card" -ne "$expected" ]; then
        status="FAIL"
        fail "cardinality $scenario (expected=$expected, legacy=$legacy_card, current=$current_card)"
    else
        pass "cardinality $scenario = $expected"
    fi
    echo "| \`$scenario\` | $expected | $legacy_card | $current_card | $status |" >> "$REPORT"
done < "$EXPECTED_DIR/cardinalities.txt"

# --- Paper Figure 3 spot check: subgraphs 1364 and 2388 → same element in VCP_{4,1,1} ---
echo
echo "=== Phase 1b: paper Figure 3 isomorphism spot check ==="
{
    echo
    echo "### Figure 3 spot check — VCP_{4,1,1}"
    echo
    echo "The paper shows that subgraphs 1364 and 2388 are isomorphic with"
    echo "respect to the fixed vertex pair and therefore map to the same"
    echo "VCP element."
    echo
    echo "| Subgraph | Legacy element | Current element |"
    echo "|---|---|---|"
} >> "$REPORT"

legacy_map_411="$RESULTS_DIR/paper/legacy_vcp_map_4_1_1.txt"
current_map_411="$RESULTS_DIR/paper/current_vcp_map_4_1_1.txt"
# Lines are 1-indexed in sed; subgraph IDs in the paper are 0-indexed.
legacy_1364=$(sed -n '1365p' "$legacy_map_411")
current_1364=$(sed -n '1365p' "$current_map_411")
legacy_2388=$(sed -n '2389p' "$legacy_map_411")
current_2388=$(sed -n '2389p' "$current_map_411")

{
    echo "| 1364 | $legacy_1364 | $current_1364 |"
    echo "| 2388 | $legacy_2388 | $current_2388 |"
} >> "$REPORT"

if [ "$legacy_1364" = "$legacy_2388" ] && [ "$current_1364" = "$current_2388" ]; then
    pass "Figure 3: subgraph 1364 and 2388 map to same element in both builds"
else
    fail "Figure 3: isomorphism check (legacy 1364=$legacy_1364, 2388=$legacy_2388; current 1364=$current_1364, 2388=$current_2388)"
fi

# =============================================================================
# Phase 2 — Characterization: legacy vs current on vcp_map full outputs
# =============================================================================
echo
echo "=== Phase 2: vcp_map characterization (legacy byte-diff current) ==="
{
    echo
    echo "## Phase 2 — vcp_map characterization (legacy vs current)"
    echo
    echo "Byte-for-byte diff of \`vcp_map\` output between the two builds."
    echo
    echo "| Scenario | Output lines | Diff bytes | Result |"
    echo "|---|---|---|---|"
} >> "$REPORT"

while read -r line; do
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    read -r n r d _expected <<< "$line"
    scenario="vcp_map $n $r $d"
    legacy_out="$RESULTS_DIR/paper/legacy_vcp_map_${n}_${r}_${d}.txt"
    current_out="$RESULTS_DIR/paper/current_vcp_map_${n}_${r}_${d}.txt"
    lines=$(wc -l < "$legacy_out")
    diff_bytes=$(diff "$legacy_out" "$current_out" 2>/dev/null | wc -c || echo "?")
    if cmp -s "$legacy_out" "$current_out"; then
        echo "| \`$scenario\` | $lines | 0 | PASS |" >> "$REPORT"
        pass "vcp_map $n $r $d (bytes match)"
    else
        echo "| \`$scenario\` | $lines | $diff_bytes | FAIL |" >> "$REPORT"
        fail "vcp_map $n $r $d: legacy and current differ"
    fi
done < "$EXPECTED_DIR/cardinalities.txt"

# =============================================================================
# Phase 3 — vcp_generate characterization
# =============================================================================
echo
echo "=== Phase 3: vcp_generate characterization ==="
{
    echo
    echo "## Phase 3 — vcp_generate characterization"
    echo
    echo "For each (graph fixture, n, r, d) combination, feed the fixture and"
    echo "its pair list to both builds and diff the VCP-vector outputs."
    echo
    echo "| Graph | n | r | d | Pairs | Result |"
    echo "|---|---|---|---|---|---|"
} >> "$REPORT"

# Only (n, r, d) combinations the CLI actually dispatches on:
VCP_GEN_PARAMS=(
    "3 1 0"
    "3 2 0"
    "3 1 1"
    "3 2 1"
    "4 1 0"
    "4 2 0"
    "4 1 1"
    "4 2 1"
)

for fixture_path in "$FIXTURES_DIR"/graphs_undirected/*.txt; do
    fixture_name=$(basename "$fixture_path" .txt)
    pair_file="$FIXTURES_DIR/pairs/${fixture_name}.txt"
    [ -f "$pair_file" ] || continue

    for params in "${VCP_GEN_PARAMS[@]}"; do
        read -r n r d <<< "$params"
        # d=0 requires unirelational or multirelational-r undirected graph.
        # For r=1 we pass the plain file. For r=2 we need the multirelational fixture.
        # For d=1 we need the directed fixture counterpart (if it exists).
        graph_for_tool="$fixture_path"
        if [ "$r" = "2" ]; then
            graph_for_tool="$FIXTURES_DIR/graphs_multirelational_r2/${fixture_name}.txt"
            [ -f "$graph_for_tool" ] || continue
        fi
        if [ "$d" = "1" ]; then
            # Map undirected name to directed counterpart where we have one.
            dname=""
            case "$fixture_name" in
                k3_triangle) dname="dk3_bidirectional" ;;
                k4_complete) dname="dk4_bidirectional" ;;
                p4_path)     dname="dp4_forward" ;;
                c4_cycle)    dname="dcycle4" ;;
                *) continue ;;
            esac
            graph_for_tool="$FIXTURES_DIR/graphs_directed/${dname}.txt"
            [ -f "$graph_for_tool" ] || continue
            if [ "$r" = "2" ]; then
                # No multirelational directed fixtures in this harness.
                continue
            fi
        fi

        tag="${fixture_name}_n${n}_r${r}_d${d}"
        legacy_out="$RESULTS_DIR/characterization/legacy_vcp_generate_${tag}.txt"
        current_out="$RESULTS_DIR/characterization/current_vcp_generate_${tag}.txt"

        "$LEGACY_BIN/vcp_generate" "$n" "$r" "$d" "$graph_for_tool" \
            < "$pair_file" > "$legacy_out" 2>/dev/null || true
        "$CURRENT_BIN/vcp_generate" "$n" "$r" "$d" "$graph_for_tool" \
            < "$pair_file" > "$current_out" 2>/dev/null || true

        pair_count=$(wc -l < "$pair_file")
        if cmp -s "$legacy_out" "$current_out"; then
            echo "| $fixture_name | $n | $r | $d | $pair_count | PASS |" >> "$REPORT"
            pass "vcp_generate $tag"
        else
            echo "| $fixture_name | $n | $r | $d | $pair_count | FAIL |" >> "$REPORT"
            fail "vcp_generate $tag"
        fi
    done
done

# =============================================================================
# Phase 4 — directed_to_undirected characterization
# =============================================================================
echo
echo "=== Phase 4: directed_to_undirected characterization ==="
{
    echo
    echo "## Phase 4 — directed_to_undirected characterization"
    echo
    echo "| Fixture | -b flag | Result |"
    echo "|---|---|---|"
} >> "$REPORT"

for fixture_path in "$FIXTURES_DIR"/graphs_directed/*.txt; do
    fixture_name=$(basename "$fixture_path" .txt)
    for flag in "" "-b"; do
        flag_label="${flag:-none}"
        tag="${fixture_name}_${flag_label}"
        legacy_out="$RESULTS_DIR/characterization/legacy_d2u_${tag}.txt"
        current_out="$RESULTS_DIR/characterization/current_d2u_${tag}.txt"

        "$LEGACY_BIN/directed_to_undirected" $flag < "$fixture_path" > "$legacy_out" 2>/dev/null || true
        "$CURRENT_BIN/directed_to_undirected" $flag < "$fixture_path" > "$current_out" 2>/dev/null || true

        if cmp -s "$legacy_out" "$current_out"; then
            echo "| $fixture_name | $flag_label | PASS |" >> "$REPORT"
            pass "directed_to_undirected $tag"
        else
            echo "| $fixture_name | $flag_label | FAIL |" >> "$REPORT"
            fail "directed_to_undirected $tag"
        fi
    done
done

# =============================================================================
# Phase 5 — ell_2_pairs characterization
# =============================================================================
echo
echo "=== Phase 5: ell_2_pairs characterization ==="
{
    echo
    echo "## Phase 5 — ell_2_pairs characterization"
    echo
    echo "| Fixture | Result |"
    echo "|---|---|"
} >> "$REPORT"

for fixture_path in "$FIXTURES_DIR"/graphs_undirected/*.txt; do
    fixture_name=$(basename "$fixture_path" .txt)
    legacy_out="$RESULTS_DIR/characterization/legacy_ell_${fixture_name}.txt"
    current_out="$RESULTS_DIR/characterization/current_ell_${fixture_name}.txt"

    "$LEGACY_BIN/ell_2_pairs" < "$fixture_path" > "$legacy_out" 2>/dev/null || true
    "$CURRENT_BIN/ell_2_pairs" < "$fixture_path" > "$current_out" 2>/dev/null || true

    if cmp -s "$legacy_out" "$current_out"; then
        echo "| $fixture_name | PASS |" >> "$REPORT"
        pass "ell_2_pairs $fixture_name"
    else
        echo "| $fixture_name | FAIL |" >> "$REPORT"
        fail "ell_2_pairs $fixture_name"
    fi
done

# =============================================================================
# Summary
# =============================================================================
TOTAL=$((PASS_COUNT + FAIL_COUNT + SKIP_COUNT))
{
    echo
    echo "## Summary"
    echo
    echo "- Total scenarios: $TOTAL"
    echo "- Passed:  $PASS_COUNT"
    echo "- Failed:  $FAIL_COUNT"
    echo "- Skipped: $SKIP_COUNT"
    echo
    if [ "$FAIL_COUNT" -ne 0 ]; then
        echo "### Failures"
        echo
        for note in "${FAIL_NOTES[@]}"; do
            echo "- $note"
        done
    fi
} >> "$REPORT"

echo
echo "=== Done ==="
echo "Total: $TOTAL, Pass: $PASS_COUNT, Fail: $FAIL_COUNT, Skip: $SKIP_COUNT"
echo "Report: $REPORT"

[ "$FAIL_COUNT" -eq 0 ]

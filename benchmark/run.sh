#!/usr/bin/env bash
# Cross-ref benchmark runner for the vcp CLI tools.
#
# Times a fixed set of realistic CLI workloads against one or more git
# revisions and emits a markdown comparison table. The workloads are
# chosen to exercise the major code paths - undirected and directed VCP
# enumeration at subgraph sizes 3 and 4, multirelational enumeration, and
# the I/O-dominated auxiliary tools.
#
# Fixtures are curated to avoid the known legacy bugs (infinite loops in
# directed_to_undirected on asymmetric neighbor sets, out-of-bounds reads
# in multirelational edge_value on missing edges), so the same inputs run
# on both the 2012 original code and the modernized code for
# apples-to-apples comparison. Any workload that still fails for a ref
# is reported as CRASH/TIMEOUT in the result table rather than aborting
# the whole sweep.
#
# Usage:
#   ./run.sh [--refs REVSPEC[,REVSPEC...]] [--runs N] [--timeout SEC]
#            [--skip-build] [--no-legacy]
#
# Defaults:
#   --refs "develop"      one ref (legacy added implicitly unless --no-legacy)
#   --runs 3              repetitions per workload per ref
#   --timeout 300         per-invocation wall-clock seconds
#
# Output:
#   results/<timestamp>/comparison.md    comparison table
#   results/<timestamp>/raw/             per-run timing files
#   results/latest -> <timestamp>        convenience symlink
#
# Exit code: 0 if all workloads produced a timing for at least one ref.

set -euo pipefail

# The timing primitive below uses $EPOCHREALTIME (bash 5.0+). Guarding
# here keeps the failure mode obvious - on older bash the variable
# silently expands to empty and every timing would become 0.000.
if [ "${BASH_VERSINFO[0]}" -lt 5 ]; then
    echo "ERROR: bash 5+ required (for \$EPOCHREALTIME); found ${BASH_VERSION}" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FIXTURES_DIR="$SCRIPT_DIR/fixtures"
RESULTS_ROOT="$SCRIPT_DIR/results"

# --- Arg parsing ---
REFS_ARG="develop"
RUNS=3
TIMEOUT=300
SKIP_BUILD=0
NO_LEGACY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --refs)       REFS_ARG="$2"; shift 2 ;;
        --runs)       RUNS="$2";     shift 2 ;;
        --timeout)    TIMEOUT="$2";  shift 2 ;;
        --skip-build) SKIP_BUILD=1;  shift ;;
        --no-legacy)  NO_LEGACY=1;   shift ;;
        -h|--help)
            sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

# Require odd $RUNS. The median-of-N computation below selects the
# ((N+1)/2)-th sorted value, which is the true median only when N is
# odd; for even N it returns the lower of the two middle values, which
# would be a silent misreport. Rather than averaging (ties in timing
# measurements are essentially impossible anyway), enforce the
# invariant.
if (( RUNS < 1 || RUNS % 2 == 0 )); then
    echo "ERROR: --runs must be a positive odd integer; got $RUNS" >&2
    exit 2
fi

# Build the ref list. Legacy is prepended unless opted out so the common
# case (compare current branch against the 2012 baseline) is one command.
IFS=',' read -ra REFS <<< "$REFS_ARG"
if [ "$NO_LEGACY" -eq 0 ]; then
    REFS=("ba1592d" "${REFS[@]}")
fi

# Generate fixtures if missing.
if [ ! -d "$FIXTURES_DIR/undirected" ]; then
    echo ">> Generating fixtures"
    python3 "$SCRIPT_DIR/generate_fixtures.py"
fi

# --- Slug derivation (must match build-ref.sh) ---
slug_for() {
    local rev="$1"
    if [ "$rev" = "ba1592d" ]; then
        echo "legacy"
    else
        echo "$rev" | tr '/' '_' | tr -cd 'A-Za-z0-9_.-'
    fi
}

# --- Build each ref ---
declare -a SLUGS
for ref in "${REFS[@]}"; do
    slug="$(slug_for "$ref")"
    SLUGS+=("$slug")
    if [ "$SKIP_BUILD" -eq 1 ] && [ -x "$SCRIPT_DIR/bin/$slug/vcp_generate" ]; then
        echo ">> Skipping build for $ref (slug: $slug) - binaries present"
        continue
    fi
    bash "$SCRIPT_DIR/build-ref.sh" "$ref" "$slug" >/dev/null
    echo ">> Built $ref (slug: $slug)"
done

# --- Results scaffolding ---
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RESULTS_DIR="$RESULTS_ROOT/$TIMESTAMP"
RAW_DIR="$RESULTS_DIR/raw"
mkdir -p "$RAW_DIR"

# Overwrite the "latest" symlink.
ln -sfn "$TIMESTAMP" "$RESULTS_ROOT/latest"

REPORT="$RESULTS_DIR/comparison.md"

# --- Timer primitive ---
# Uses bash 5's $EPOCHREALTIME (microsecond-precision wall clock) instead
# of /usr/bin/time (GNU time), which is not installed on every host. Each
# repetition's elapsed seconds are formatted to 3 decimal places via awk.
#
# _run_once <input_file> <cmd...>
#   Runs <cmd...> with <input_file> on stdin and /dev/null on stdout,
#   under the `timeout` guard. Echoes either the elapsed seconds or one
#   of: TIMEOUT / CRASH-<rc>. Does not read stderr from the tool (we
#   don't want tool diagnostics to contaminate the timing output).
_run_once() {
    local input="$1"; shift
    local start end rc
    start="$EPOCHREALTIME"
    # Note: bash `if C; then T; fi` returns 0 when C fails and there is no
    # else branch, masking the real exit code. Capture it inline with `|| rc=$?`
    # instead. `set -e` is also satisfied by this pattern.
    rc=0
    timeout "$TIMEOUT" "$@" < "$input" >/dev/null 2>/dev/null || rc=$?
    end="$EPOCHREALTIME"
    case "$rc" in
        0)   awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f\n", e - s }' ;;
        124) echo "TIMEOUT" ;;
        *)   echo "CRASH-$rc" ;;
    esac
    # The string output above IS the interface. Callers read the captured
    # stdout; they do not test $? of this function. Return 0 unconditionally
    # so a partial run does not trip set -e in the outer shell.
    return 0
}

# _collect_runs <label> <slug> <input_file> <cmd...>
#   Runs <cmd...> $RUNS times, writes per-run results to the raw log, and
#   echoes the median seconds (or TIMEOUT/CRASH if any run failed -
#   silent partial successes would hide a failure mode).
_collect_runs() {
    local label="$1"; shift
    local slug="$1";  shift
    local input="$1"; shift
    local raw="$RAW_DIR/${label}__${slug}.log"
    : > "$raw"

    local outcomes=()
    for _ in $(seq 1 "$RUNS"); do
        local result
        result="$(_run_once "$input" "$@")"
        echo "$result" >> "$raw"
        case "$result" in
            TIMEOUT*|CRASH*)
                echo "$result"
                return
                ;;
        esac
        outcomes+=("$result")
    done
    printf '%s\n' "${outcomes[@]}" | sort -n | awk -v n="$RUNS" 'NR==int((n+1)/2) {print; exit}'
}

# --- Pipe-style timing (for vcp_generate, which reads pairs from stdin) ---
# time_pipe <label> <slug> <tool> <pairs_file> <tool args...>
time_pipe() {
    local label="$1"; shift
    local slug="$1";  shift
    local tool="$1";  shift
    local pairs="$1"; shift
    local bin_dir="$SCRIPT_DIR/bin/$slug"

    if [ ! -x "$bin_dir/$tool" ]; then
        echo "MISSING"
        return
    fi
    _collect_runs "$label" "$slug" "$pairs" "$bin_dir/$tool" "$@"
}

# --- Stdin-redirect timing (for directed_to_undirected, ell_2_pairs) ---
# time_stdin <label> <slug> <tool> <input_file> [tool args...]
time_stdin() {
    local label="$1"; shift
    local slug="$1";  shift
    local tool="$1";  shift
    local input="$1"; shift
    local bin_dir="$SCRIPT_DIR/bin/$slug"

    if [ ! -x "$bin_dir/$tool" ]; then
        echo "MISSING"
        return
    fi
    _collect_runs "$label" "$slug" "$input" "$bin_dir/$tool" "$@"
}

# --- Report header ---
{
    echo "# vcp benchmark comparison"
    echo
    echo "Generated: $TIMESTAMP"
    echo "Runs per workload: $RUNS (median reported)"
    echo "Timeout: ${TIMEOUT}s"
    echo
    echo "## Refs"
    echo
    echo "| Slug | Revspec | Commit |"
    echo "|------|---------|--------|"
    for slug in "${SLUGS[@]}"; do
        local_ref="$(cat "$SCRIPT_DIR/bin/$slug/.revspec" 2>/dev/null || echo '?')"
        local_sha="$(cat "$SCRIPT_DIR/bin/$slug/.ref" 2>/dev/null || echo '?')"
        printf '| %s | %s | %s |\n' "$slug" "$local_ref" "${local_sha:0:10}"
    done
    echo
    echo "## Timings (seconds, median of $RUNS runs)"
    echo
    # Column header
    printf '| Workload |'
    for slug in "${SLUGS[@]}"; do
        printf ' %s |' "$slug"
    done
    echo
    printf '%s' '|----------|'
    for _ in "${SLUGS[@]}"; do
        printf '%s' '--------|'
    done
    echo
} > "$REPORT"

# --- Workload table ---
# Each workload is:
#   <label>                  a short human-readable tag used in the report
#   <kind>                   pipe|stdin
#   <tool> <args...>         what to run
#   <input>                  file fed on stdin (pairs for pipe, graph for stdin)
#
# The `pipe` kind passes tool args + a pairs file. The `stdin` kind passes
# tool args + a graph file.

run_workload() {
    local label="$1"; local kind="$2"; shift 2
    printf '| `%s` |' "$label" >> "$REPORT"
    for slug in "${SLUGS[@]}"; do
        echo ">>   [$slug] $label" >&2
        local result
        if [ "$kind" = "pipe" ]; then
            # pipe args: <tool> <pairs> <tool_args...>
            result="$(time_pipe "$label" "$slug" "$@")"
        else
            # stdin args: <tool> <input> <tool_args...>
            result="$(time_stdin "$label" "$slug" "$@")"
        fi
        printf ' %s |' "$result" >> "$REPORT"
    done
    echo >> "$REPORT"
}

# vcp_generate (the primary compute workload)
# Arguments: vcp_generate n r d graph_filename < pairs
echo ">> Running vcp_generate workloads" >&2

run_workload "vcp_gen 3 1 0 n=200"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n200_p0.10_s42.pairs" \
    3 1 0 "$FIXTURES_DIR/undirected/er_n200_p0.10_s42.txt"

run_workload "vcp_gen 3 1 0 n=500"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n500_p0.05_s42.pairs" \
    3 1 0 "$FIXTURES_DIR/undirected/er_n500_p0.05_s42.txt"

run_workload "vcp_gen 3 1 0 n=1000" pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n1000_p0.02_s42.pairs" \
    3 1 0 "$FIXTURES_DIR/undirected/er_n1000_p0.02_s42.txt"

run_workload "vcp_gen 4 1 0 n=200"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n200_p0.10_s42.pairs" \
    4 1 0 "$FIXTURES_DIR/undirected/er_n200_p0.10_s42.txt"

run_workload "vcp_gen 4 1 0 n=500"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n500_p0.05_s42.pairs" \
    4 1 0 "$FIXTURES_DIR/undirected/er_n500_p0.05_s42.txt"

run_workload "vcp_gen 4 1 0 n=1000" pipe vcp_generate \
    "$FIXTURES_DIR/pairs/er_n1000_p0.02_s42.pairs" \
    4 1 0 "$FIXTURES_DIR/undirected/er_n1000_p0.02_s42.txt"

# Directed: use bidirectionally-symmetric fixtures so legacy
# directed_to_undirected would work (legacy vcp_generate is fine on these
# regardless because the bug is in a different tool).
run_workload "vcp_gen 3 1 1 n=200"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/dsym_n200_p0.10_s42.pairs" \
    3 1 1 "$FIXTURES_DIR/directed/dsym_n200_p0.10_s42.txt"

run_workload "vcp_gen 3 1 1 n=500"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/dsym_n500_p0.05_s42.pairs" \
    3 1 1 "$FIXTURES_DIR/directed/dsym_n500_p0.05_s42.txt"

run_workload "vcp_gen 4 1 1 n=200"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/dsym_n200_p0.10_s42.pairs" \
    4 1 1 "$FIXTURES_DIR/directed/dsym_n200_p0.10_s42.txt"

run_workload "vcp_gen 4 1 1 n=500"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/dsym_n500_p0.05_s42.pairs" \
    4 1 1 "$FIXTURES_DIR/directed/dsym_n500_p0.05_s42.txt"

# Multirelational r=2 (dense so legacy does not hit a missing-edge read).
run_workload "vcp_gen 3 2 0 n=100"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/mr_er_n100_p0.30_s42.pairs" \
    3 2 0 "$FIXTURES_DIR/multirelational/er_n100_p0.30_s42.txt"

run_workload "vcp_gen 4 2 0 n=100"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/mr_er_n100_p0.30_s42.pairs" \
    4 2 0 "$FIXTURES_DIR/multirelational/er_n100_p0.30_s42.txt"

run_workload "vcp_gen 4 2 0 n=200"  pipe vcp_generate \
    "$FIXTURES_DIR/pairs/mr_er_n200_p0.15_s42.pairs" \
    4 2 0 "$FIXTURES_DIR/multirelational/er_n200_p0.15_s42.txt"

# directed_to_undirected: symmetric fixture so legacy's unfixed loop
# terminates. -b mode exercises the intersection branch only.
echo ">> Running auxiliary tools" >&2
run_workload "d2u -b n=1000"   stdin directed_to_undirected \
    "$FIXTURES_DIR/directed/dsym_n1000_p0.02_s42.txt" -b

run_workload "d2u none n=1000" stdin directed_to_undirected \
    "$FIXTURES_DIR/directed/dsym_n1000_p0.02_s42.txt"

# ell_2_pairs
run_workload "ell_2_pairs n=1000" stdin ell_2_pairs \
    "$FIXTURES_DIR/undirected/er_n1000_p0.02_s42.txt"

run_workload "ell_2_pairs n=2000" stdin ell_2_pairs \
    "$FIXTURES_DIR/undirected/er_n2000_p0.01_s42.txt"

# Footer
{
    echo
    echo "## Interpreting results"
    echo
    echo "- Cells are wall-clock seconds (smaller is better)."
    echo "- \`CRASH\` / \`TIMEOUT\` / \`MISSING\` replace the timing when a run failed."
    echo "- Raw per-run timings are under \`raw/<workload>__<slug>.log\`."
} >> "$REPORT"

echo
echo "=== Done ==="
echo "Report: $REPORT"
echo "Latest: $RESULTS_ROOT/latest"

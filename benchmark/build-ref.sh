#!/usr/bin/env bash
# Build the vcp CLI tools at a given git revision into benchmark/bin/<slug>/.
#
# Auto-detects the build system by checking for CMakeLists.txt at the root:
#   * modern (has CMakeLists.txt): cmake -B build -DCMAKE_BUILD_TYPE=Release
#   * legacy (no CMakeLists.txt):  make with -Werror stripped
#
# Each ref is built in its own detached-HEAD worktree so the main working
# tree HEAD is never touched.
#
# Usage:
#   ./build-ref.sh <revspec> [<slug>]
#
# If <slug> is omitted, it's derived from the revspec: branches/tags keep
# their name (with / replaced by _), commit SHAs use the short hash.
# <slug> is sanitized to [A-Za-z0-9_.-] to play nicely with filesystem paths.
#
# Examples:
#   ./build-ref.sh main
#   ./build-ref.sh ba1592d legacy
#   ./build-ref.sh feature/benchmark-suite

set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <revspec> [<slug>]" >&2
    exit 2
fi

REV="$1"
SLUG="${2:-}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"

# Derive slug if not provided.
if [ -z "$SLUG" ]; then
    SLUG="$REV"
fi
# Sanitize slug for filesystem use.
SLUG="$(echo "$SLUG" | tr '/' '_' | tr -cd 'A-Za-z0-9_.-')"

if [ -z "$SLUG" ]; then
    echo "ERROR: slug derivation yielded an empty string for revspec '$REV'" >&2
    exit 2
fi

WORKTREE_DIR="$SCRIPT_DIR/worktree-$SLUG"
BIN_OUT="$SCRIPT_DIR/bin/$SLUG"

echo ">> Building vcp at revision $REV (slug: $SLUG)"
echo "   repo:     $REPO_ROOT"
echo "   worktree: $WORKTREE_DIR"
echo "   output:   $BIN_OUT"

# Resolve revspec to a concrete SHA so the worktree uses a stable snapshot.
REV_SHA="$(git -C "$REPO_ROOT" rev-parse --verify "$REV^{commit}")"

# Idempotent: remove any stale worktree before creating a fresh one.
# If the directory exists but is not a tracked worktree, fail loudly
# rather than `rm -rf`ing an unrelated directory that happens to share
# the name - the slug is sanitized but a typo in `bin/` or a leftover
# from a manual experiment could otherwise be silently deleted.
if [ -d "$WORKTREE_DIR" ]; then
    echo ">> Removing stale worktree at $WORKTREE_DIR"
    if ! git -C "$REPO_ROOT" worktree remove --force "$WORKTREE_DIR"; then
        echo "ERROR: $WORKTREE_DIR exists but is not a git worktree." >&2
        echo "       Refusing to delete it automatically. Inspect and remove manually." >&2
        exit 1
    fi
fi

git -C "$REPO_ROOT" worktree add --detach "$WORKTREE_DIR" "$REV_SHA"

if [ -f "$WORKTREE_DIR/CMakeLists.txt" ]; then
    echo ">> Modern build: CMake (Release)"
    cmake -S "$WORKTREE_DIR" -B "$WORKTREE_DIR/build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DVCP_BUILD_TESTS=OFF \
        >/dev/null
    cmake --build "$WORKTREE_DIR/build" -j "$(nproc)" >/dev/null

    mkdir -p "$BIN_OUT"
    for tool in vcp_generate vcp_map directed_to_undirected ell_2_pairs; do
        cp "$WORKTREE_DIR/build/$tool" "$BIN_OUT/$tool"
    done
else
    echo ">> Legacy build: Makefile with -Werror stripped"
    # Matches regression/build-legacy.sh. The legacy Makefile includes
    # -Werror in its COMMON_FLAGS; GCC 14 emits warnings the 2012 GCC 4.8
    # didn't, so we override the flag set entirely.
    MAX_NEIGHBORS=16384
    CPP_FLAGS="-Wall -Wextra -Wno-unused-local-typedefs -std=c++11 -pedantic"
    CPP_FLAGS="$CPP_FLAGS -I ./inc -I ./lib/boost_1_53_0 -I ./lib/tclap-1.2.1/include"
    CPP_FLAGS="$CPP_FLAGS -D MAX_NEIGHBORS=${MAX_NEIGHBORS} -O3 -D NDEBUG"

    make -C "$WORKTREE_DIR" -j "$(nproc)" CPP_FLAGS="$CPP_FLAGS" >/dev/null

    mkdir -p "$BIN_OUT"
    for tool in vcp_generate vcp_map directed_to_undirected ell_2_pairs; do
        cp "$WORKTREE_DIR/bin/$tool" "$BIN_OUT/$tool"
    done
fi

# Drop a sidecar file noting what ref these binaries were built from, so
# results reports can cite the exact SHA even after the branch moves.
echo "$REV_SHA" > "$BIN_OUT/.ref"
echo "$REV"     > "$BIN_OUT/.revspec"

echo ">> Binaries installed to $BIN_OUT"
ls -la "$BIN_OUT"

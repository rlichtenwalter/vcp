#!/usr/bin/env bash
# Build legacy vcp binaries at a specified git revision into bin/legacy/.
#
# The legacy code uses the 2013-era Makefile + vendored Boost 1.53 + TCLAP
# 1.2.1 and is compiled as-is. Two concessions are made to let it build on
# a 2026 toolchain (GCC 14) without touching its source:
#
#   1. `-Werror` is stripped so GCC 14's newer warnings (which didn't exist
#      in GCC 4.8) don't block compilation. The warnings themselves are
#      legitimate — they're what the modernization pass fixed — but here
#      we're building the LEGACY code to characterize its output, not to
#      hold it to modern clean-compile standards.
#   2. `-std=c++11` is kept; the legacy code is C++11.
#
# Uses `git worktree add` so the main working tree HEAD is never moved.
# The worktree is placed under regression/worktree-legacy/ and cleaned up
# on re-runs.

set -euo pipefail

REV="${1:-ba1592d}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
WORKTREE_DIR="$SCRIPT_DIR/worktree-legacy"
BIN_OUT="$SCRIPT_DIR/bin/legacy"

echo ">> Building legacy vcp at revision $REV"
echo "   repo:     $REPO_ROOT"
echo "   worktree: $WORKTREE_DIR"
echo "   output:   $BIN_OUT"

# Remove any previous worktree before creating a new one so re-runs are idempotent.
if [ -d "$WORKTREE_DIR" ]; then
    echo ">> Removing stale worktree at $WORKTREE_DIR"
    git -C "$REPO_ROOT" worktree remove --force "$WORKTREE_DIR" 2>/dev/null || rm -rf "$WORKTREE_DIR"
fi

git -C "$REPO_ROOT" worktree add --detach "$WORKTREE_DIR" "$REV"

echo ">> Clean-compiling the legacy Makefile with -Werror stripped"
# Temporarily relax -Werror so GCC 14 warnings don't abort the build.
# Override CPP_FLAGS directly by invoking make with an overriding flag set.
# The legacy Makefile composes its flags in COMMON_FLAGS / CPP_FLAGS; we
# bypass those by passing CPP_FLAGS on the command line.
MAX_NEIGHBORS=16384
CPP_FLAGS="-Wall -Wextra -Wno-unused-local-typedefs -std=c++11 -pedantic -I ./inc -I ./lib/boost_1_53_0 -I ./lib/tclap-1.2.1/include -D MAX_NEIGHBORS=${MAX_NEIGHBORS} -O3 -D NDEBUG"

make -C "$WORKTREE_DIR" -j "$(nproc)" CPP_FLAGS="$CPP_FLAGS"

mkdir -p "$BIN_OUT"
for tool in vcp_generate vcp_map directed_to_undirected ell_2_pairs; do
    cp "$WORKTREE_DIR/bin/$tool" "$BIN_OUT/$tool"
done

echo ">> Legacy binaries installed to $BIN_OUT"
ls -la "$BIN_OUT"

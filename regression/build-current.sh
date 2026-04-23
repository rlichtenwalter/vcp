#!/usr/bin/env bash
# Build current vcp binaries at a specified git revision into bin/current/.
#
# Mirrors build-legacy.sh but for the modernized CMake build. Uses
# `git worktree add` so the main working tree HEAD is never moved.

set -euo pipefail

REV="${1:-feature/tooling-modernization}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
WORKTREE_DIR="$SCRIPT_DIR/worktree-current"
BIN_OUT="$SCRIPT_DIR/bin/current"

echo ">> Building current vcp at revision $REV"
echo "   repo:     $REPO_ROOT"
echo "   worktree: $WORKTREE_DIR"
echo "   output:   $BIN_OUT"

if [ -d "$WORKTREE_DIR" ]; then
    echo ">> Removing stale worktree at $WORKTREE_DIR"
    git -C "$REPO_ROOT" worktree remove --force "$WORKTREE_DIR" 2>/dev/null || rm -rf "$WORKTREE_DIR"
fi

git -C "$REPO_ROOT" worktree add --detach "$WORKTREE_DIR" "$REV"

echo ">> Configuring CMake (Release, compile_commands.json off)"
cmake -S "$WORKTREE_DIR" -B "$WORKTREE_DIR/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DVCP_BUILD_TESTS=OFF \
    >/dev/null

echo ">> Building"
cmake --build "$WORKTREE_DIR/build" -j "$(nproc)" >/dev/null

mkdir -p "$BIN_OUT"
for tool in vcp_generate vcp_map directed_to_undirected ell_2_pairs; do
    cp "$WORKTREE_DIR/build/$tool" "$BIN_OUT/$tool"
done

echo ">> Current binaries installed to $BIN_OUT"
ls -la "$BIN_OUT"

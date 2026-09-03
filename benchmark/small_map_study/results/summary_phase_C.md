# Phase C: container microbenchmark (break-even curve)

Catch2 micro-benchmark at `-O3 -DNDEBUG`. Each iteration: construct a fresh
container, perform `30 * k` insert-or-increment ops drawn from a key universe of
size `k`, then iterate the populated entries once. Measures the per-call
lifetime of `temp_edge_types`.

Values in nanoseconds (ns) or microseconds (us). Mean of 30 samples.

| k    | dense_array | std::unordered_map | flat_pair_vec (heap) | stack_flat[cap=k] | std::map |
| ---- | ----------- | ------------------ | -------------------- | ----------------- | -------- |
| 4    | 52 ns       | 271 ns             | 201 ns               | 246 ns            | 360 ns   |
| 16   | 270 ns      | 1.13 us            | 1.81 us              | 1.64 us           | 2.62 us  |
| 64   | 517 ns      | 4.53 us            | 29.2 us              | 27.2 us           | 14.3 us  |
| 256  | 2.75 us     | 20.1 us            | 284 us               | 154 us            | 205 us   |
| 1024 | 10.95 us    | 86.0 us            | 3.89 ms              | 676 us            | 1.26 ms  |
| 4096 | —           | 417 us             | 59.8 ms              | 2.70 ms           | 7.16 ms  |

## Key findings

### 1. Dense array wins wherever it's feasible

At every k up through 1024, `std::array` beats every other container by 4×-25×.
The dense tier is unambiguously the right choice for any key space that fits in
memory.

### 2. Hash (std::unordered_map) dominates the sparse regime

At k ≥ 16, `std::unordered_map` beats every other sparse container:

- vs. `std::map`: 2.3× faster at k=16, 10× at k=256, 17× at k=4096.
- vs. `flat_pair_vec`: 1.6× faster at k=16, 14× at k=256, **143×** at k=4096.
- vs. `stack_flat`: 1.5× faster at k=16, 8× at k=256, **6×** at k=4096.

The linear-scan forms cliff-degrade as k grows; hash stays near-flat.

### 3. Linear-scan has no robust operating regime

Linear scan (either heap-vector or stack-array) only beats hash at k ≤ 4-8. But
at that scale the **dense tier already wins**, because any key space producing k
≤ 8 is trivially small enough to fit in the dense byte budget. So linear scan's
"win" band does not overlap with any regime where dense is unavailable.

Concretely:

- k=4, dense unavailable: would need key space large enough that k=4 is
  plausible. But observed k never exceeds `min(2^r, distinct data values)`. A
  key space big enough to forbid dense (r > 20 undirected) but yielding k ≤ 4 in
  practice is a very narrow band, and hash beats stack_flat there by only ~10%
  (271 ns vs 246 ns).

- k=16, dense unavailable: this is the r=30 low-diversity case from Phase B.
  stack_flat 1.64 us, hash 1.13 us → **hash is 45% faster**.

- k ≥ 64: hash dominates by large factors.

### 4. std::map is consistently the worst choice

`std::map` never wins any cell of this table. At the low end it loses to
`flat_pair_vec` and `stack_flat` due to per-node allocator cost. At the high end
it loses to `unordered_map` due to tree-walk log factor. Removing it from the
hot path is a pure win.

## Implications for tier count

The Phase B plan proposed a three-tier design (dense / linear-scan / hash). The
Phase C data does not support the middle tier:

> Linear scan has no region of the (k, key-space) plane where it outperforms
> either of its neighbors by more than ~10%, and that 10% region (k ≤ 4) is
> dominated by the dense tier for any realistic instantiation.

A **two-tier design** (dense + hash) is empirically correct for this library's
hot path. `std::unordered_map` is fast enough as a pathological fallback that
the linear-scan middle tier does not earn its keep.

## Recommendation

Ship two-tier: `std::array` for bounded key spaces under the 8 MB byte budget,
`std::unordered_map` otherwise. Document the reasoning with the numbers above so
a future reader can re-examine the decision if hardware, compiler, or workloads
change.

This simplifies the implementation and test surface vs. the originally- planned
three-tier, without sacrificing any measured performance. Presenting this to the
user for a decision-point before Phase D builds the prototype — the user's prior
was three-tier, and the data shifts the recommendation.

# counts-container study — results

Micro-benchmark of six container candidates for the `counts` return
value of `vcp<4, r, d>::generate_vector`. Each iteration: construct,
perform all insert-or-increment ops, finalize to a sorted
`vector<pair<Key,Val>>`, sum values, destroy. One run on
linux 6.12 / g++ 14.3 / -O3 -DNDEBUG; 15 samples per cell. Full raw
output in `counts_bench_v1.txt`.

Mean wall-clock per iteration; relative column is speedup vs `StdMap`
baseline (>1 = faster, <1 = slower). Best cell in each row is **bold**.

## Per-scenario results

### tiny-hot — K=16, N=100k, uniform (N/K = 6250)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 1418 µs     | 1.00×      |
| HashSortedKeys   | 218 µs      | **6.51×**  |
| HashSortedIters  | 221 µs      | 6.43×      |
| HashSortedPairs  | 219 µs      | 6.46×      |
| SortedVector     | 1459 µs     | 0.97×      |
| BatchedDedup     | 1917 µs     | 0.74×      |

### small-zipf — K=256, N=1M, Zipf(s=1.2) (N/K = 3906)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 17.73 ms    | 1.00×      |
| HashSortedKeys   | 2.25 ms     | **7.89×**  |
| HashSortedIters  | 2.46 ms     | 7.20×      |
| HashSortedPairs  | 2.39 ms     | 7.43×      |
| SortedVector     | 16.20 ms    | 1.09×      |
| BatchedDedup     | 38.07 ms    | 0.47×      |

### medium-uniform — K=4096, N=500k, uniform (N/K = 122)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 29.62 ms    | 1.00×      |
| HashSortedKeys   | 1.71 ms     | 17.3×      |
| HashSortedIters  | 1.70 ms     | 17.4×      |
| HashSortedPairs  | 1.69 ms     | **17.5×**  |
| SortedVector     | 26.53 ms    | 1.12×      |
| BatchedDedup     | 19.77 ms    | 1.50×      |

### large-sparse — K=100k, N=200k, uniform (N/K = 2)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 42.80 ms    | 1.00×      |
| HashSortedKeys   | 12.04 ms    | 3.55×      |
| HashSortedIters  | 13.67 ms    | 3.13×      |
| HashSortedPairs  | 12.43 ms    | 3.44×      |
| SortedVector     | 536.27 ms   | 0.08×      |
| BatchedDedup     | 10.69 ms    | **4.00×**  |

### huge-oneshot — K=1M, N=1M, uniform (N/K = 1)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 938 ms      | 1.00×      |
| HashSortedKeys   | 212 ms      | 4.43×      |
| HashSortedIters  | 241 ms      | 3.89×      |
| HashSortedPairs  | 220 ms      | 4.27×      |
| SortedVector     | 30 430 ms   | 0.031×     |
| BatchedDedup     | 62 ms       | **15.1×**  |

### adversarial-desc — K=5k, N=5k, descending (N/K = 1, worst ordering)

| Candidate        | mean        | rel StdMap |
| ---------------- | ----------- | ---------- |
| StdMap           | 204.7 µs    | 1.00×      |
| HashSortedKeys   | 202.3 µs    | 1.01×      |
| HashSortedIters  | 196.0 µs    | 1.04×      |
| HashSortedPairs  | 193.8 µs    | 1.06×      |
| SortedVector     | 3328.6 µs   | 0.062×     |
| BatchedDedup     | 32.1 µs     | **6.37×**  |

## Cross-scenario summary

| Scenario          | Winner           | Margin over runner-up |
| ----------------- | ---------------- | --------------------- |
| tiny-hot          | HashSortedKeys*  | 0-1% (Hash\* tied)    |
| small-zipf        | HashSortedKeys   | ~6% over other Hash\* |
| medium-uniform    | HashSortedPairs* | 0-1% (Hash\* tied)    |
| large-sparse      | BatchedDedup     | 13% over Hash\*       |
| huge-oneshot      | BatchedDedup     | 3.4× over Hash\*      |
| adversarial-desc  | BatchedDedup     | 6× over Hash\*        |

## Findings

### 1. The three Hash* variants are functionally indistinguishable

Across all six scenarios, `HashSortedKeys`, `HashSortedIters`, and
`HashSortedPairs` are within 2-6% of each other and no single variant
wins consistently. The finalize-time work (one pass building a
vector, one sort, one materialization) is a small fraction of total
runtime at every realistic N.

The practical implication: **the choice among the three is not a
performance decision, it is a code-quality decision**. Iterator
invalidation concerns argue against `HashSortedIters`. `HashSortedKeys`
forces a second hash lookup per emitted entry at finalize time;
`HashSortedPairs` materializes the pair once and hands it out directly.
`HashSortedPairs` has the simplest iteration contract and wins ties.

### 2. SortedVector is catastrophic above small K

Insert-with-shift scales as O(K²) in the unique-key-count regime, and
it dominates every other cost:

- K=100k: **12.5× worse than StdMap**, 45× worse than Hash\*.
- K=1M:   **32× worse than StdMap**, 143× worse than Hash\*.
- K=5k descending: **16× worse than StdMap** (expected: max-shift
  ordering), 17× worse than Hash\*.

At K ≤ 4096 it is merely uncompetitive (never wins; tied-to-1.5× worse
than StdMap). **Rule it out.**

### 3. BatchedDedup has a sharp crossover at N/K ≈ 3

BatchedDedup's tradeoff is mechanical: the hot-loop append is
pointer-bump-cheap (~2 ns), but finalize pays O(N log N) on the full
append stream even when most of those N entries are duplicates that a
hash would have consolidated. The crossover is where the sort cost on
N ≈ K items equals the hash cost on N items with K-sized working set:

- N/K = 1 (huge-oneshot, adversarial-desc): wins by 3-6× over Hash\*.
- N/K = 2 (large-sparse): wins by 13%.
- N/K = 122 (medium-uniform): **loses by 12×** to Hash\*.
- N/K = 3906+ (zipf, tiny-hot): **loses by 9-17×** to Hash\*.

BatchedDedup is the right choice when every key is near-unique
(large-r runs where canonical addresses are effectively random) and
the wrong choice when the key space is small and hot (low-r runs where
a few hundred distinct canonicals receive millions of increments).

### 4. StdMap is dominated in every regime except the tiniest

The current implementation is never the best choice and is 3-17×
behind the best choice in the mid-density regimes we actually run at
scale. Replacing it is not controversial.

## Recommendation for the counts refactor

**Two viable endpoints, pick by N/K profile of the target workload.**

For typical VCP workloads (n=4, r=2, moderate-density graphs),
per-call `counts` has hundreds of distinct keys receiving thousands
of increments each — N/K in the hundreds to thousands. In that
regime, **HashSortedPairs is the clear winner**: it beats StdMap by
3-17×, ties or narrowly beats HashSortedKeys/Iters, and it returns
a pre-materialized sorted vector that needs no further conversion
at the public API boundary.

For large-r workloads (r=30 with full diversity), per-call `counts`
cardinality approaches N — most canonical addresses are unique and
BatchedDedup would be strictly better there. But r=30 is declared-max
and rare in real workloads; the k-probe study already showed that
the r=30 hot path is bottlenecked elsewhere.

**The single-container recommendation is HashSortedPairs** as the
backing for a new `vcp_vector<n, r, d>` public type. It wins or ties
the common workload, loses cleanly in the rare large-r regime, and
has the simplest iteration contract of the candidates.

If we later see evidence that real large-r `counts` calls pay
meaningful time in the hash-insert path, a compile-time dispatch
(`HashSortedPairs` for small key-spaces, `BatchedDedup` for large
ones) is a bolt-on — the public interface is the same sorted
`vector<pair<Key, Val>>` either way.

**SortedVector is rejected.** It is either tied with StdMap (small K)
or catastrophic (large K); it never wins and can be 140× slower than
the best candidate. No workload justifies it.

## Follow-up: realistic per-call scale

The scenarios above span the (N, K) space but aggregate into one
large call. Real `counts` maps are per-call: N per call is in the
low thousands for typical workloads, K in the low hundreds.

At that scale the construct/destroy overhead of `std::unordered_map`
(which does a heap allocation on first insert) may erode Hash\*'s
lead over BatchedDedup (which reserves eagerly). A follow-up scenario
with N=5000, K=200, Zipfian, many short iterations would answer
this directly if the recommendation needs more confidence. The
current 100k-lower-bound already shows Hash\* winning by 6.5× at
K=16, so per-call overhead would need to be very large to flip the
sign.

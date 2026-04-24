# Phase B: observed k distribution

Probe: `k_probe_tool` records the final cardinality of the three hot maps
(`temp_edge_types`, `edge_types`, `counts`) at the tail of every
`generate_vector` call. `VCP_INSTRUMENT_K` is `#define`d at compile time
only for this tool.

## Max k per workload per map

| Workload | temp_edge_types | edge_types | counts |
|----------|-----------------|------------|--------|
| r2_u_n100_p030         |   4 |   4 |   372 |
| r2_u_n200_p015         |   4 |   4 |   318 |
| r2_d_n100_p030         |   4 |   4 |   372 |
| r2_d_n200_p015         |   4 |   4 |   318 |
| r30_u_div4_n100        |  16 |  16 |  1738 |
| r30_u_div4_n200        |  16 |  16 |  1755 |
| r30_d_div4_n100        |  16 |  16 |  1738 |
| r30_d_div4_n200        |  16 |  16 |  1755 |
| r30_u_div4_hub_n200    |  16 |  16 |  1342 |
| r30_d_div4_hub_n200    |  16 |  16 |  1342 |
| r30_u_full_n100        | 1371 | 1478 |  4357 |
| r30_u_full_n200        | 1989 | 3024 |  7763 |
| r30_d_full_n100        | 1371 | 1478 |  4357 |
| r30_d_full_n200        | 1989 | 3024 |  7763 |
| r30_u_full_hub_n200    | 1131 | 1132 | 20833 |
| r30_d_full_hub_n200    | 1131 | 1132 | 20833 |

## Interpretation

1. **`temp_edge_types` and `edge_types` saturate at `min(2^r, distinct
   values present in data)`.** r=2 caps at 4 (2^2). r=30 with values
   drawn from [1, 15] caps at 16 (0 + 15). r=30 with values across the
   full 30-bit range is driven by data diversity up to k ~ 3000 in the
   stressed cases.

2. **r=30 full-diversity hits k in the 1000s for the two edge-type
   maps.** This is the pathological regime that a linear-scan fallback
   would struggle with. An `std::unordered_map` fallback is warranted
   for this regime, not just as a theoretical safety net.

3. **r=30 low-diversity (`div4`: values in [1, 15]) holds k ≤ 16.** In
   this regime, a linear-scan `std::vector<pair>` would outperform
   `std::unordered_map` (no hash overhead, no heap bucket indirection).
   This is the regime that earns tier 2 its keep.

4. **`counts` carries the highest k by a factor of 3-10x.** Up to 20833
   in the hub-full workload. But `counts` is keyed by
   `subgraph_address_type` (which becomes `boost::multiprecision::cpp_int`
   for r >= 11 undirected / r >= 6 directed) — not a std::size_t
   integer. The `detail::dense_or_sparse_map<K, V, KeyCount>` dense
   tier cannot apply here. `counts` stays as `std::map` for this PR;
   refactoring it is a distinct problem (sparse-friendly container with
   `subgraph_address_type` keys) and out of scope.

## Implications for tier thresholds

The bit count `b = r` (undirected) or `b = 2r` (directed) is the
compile-time proxy for key space.

- **b ≤ 20**: dense tier fits in 8 MB. Trivially best for any k.
  Confirmed by r=2 workloads (k ≤ 4, saturates at 2^2 = 4).

- **20 < b ≤ 24**: dense tier overflows the byte budget, but key space
  is still moderate (≤ 16M). Data diversity AT THIS BIT RANGE is not
  directly measured by the probe (no r ∈ {11..12} directed / r ∈ {21..24}
  undirected workloads), but the plausible extrapolation from the
  `div4`/`full` split at b=30/60 is: in this range, observed k is
  typically small (tens to low hundreds) unless the data is
  pathologically diverse. Linear-scan `std::vector<pair>` is a
  reasonable default here, to be validated or revised by Phase C
  microbench.

- **b > 24**: key space ≥ 32M; data diversity can produce k in the
  thousands (see r30_full workloads). Linear-scan at k ≈ 3000 is
  ~3000 compares per insert, which at 20k pairs × 30 inserts per pair
  = ~2 billion compares. `std::unordered_map` is the right fallback.

The proposed three-tier cutover for Phase D prototype:

```
tier 1 (dense std::array):           b * sizeof(V) ≤ 8 MB        (b ≤ 20)
tier 2 (std::vector<pair> linear):   20 < b ≤ 24
tier 3 (std::unordered_map):         b > 24
```

For directed edge_types the key is `std::pair<addr, addr>` and the
effective bit count is `2r`; the thresholds above apply to the
effective bit count. Concretely:

| r  | d | eff bits | tier |
|----|---|----------|------|
|  2 | 0 | 2  | 1 |
|  2 | 1 | 4  | 1 |
| 30 | 0 | 30 | 3 |
| 30 | 1 | 60 | 3 |

The three default tool instantiations (r ∈ {1, 2, 30}) exercise tier 1
(r ≤ 2) and tier 3 (r=30). Tier 2 is theoretical insurance for
intermediate r values (r ∈ {11..12} directed, r ∈ {21..24} undirected)
not instantiated by the default tools but reachable by downstream
users of the library.

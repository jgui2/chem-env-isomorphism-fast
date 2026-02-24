# Fast Chemical Environment Isomorphism Backend

Accelerated and memory-efficient graph isomorphism backend for chemical environment grouping workflows.

---

## Overview

This repository provides an optimized backend implementation for chemical graph isomorphism and environment deduplication used in a larger catalysis workflow.

The logic of Step2 remains structurally unchanged.  

Performance improvements come from replacing the underlying isomorphism and comparison backend.

## Context

Originally, environment grouping (Step2) relied on:

- NetworkX VF2-based isomorphism
- Nested pairwise comparison across configurations
- Graph object storage during grouping

This repository replaces this backend.

Step2 continues to call:

```python
unique_chem_envs(chem_envs_groups, metadata=paths)
```

Minor improvements were also introduced in Step2 to reduced temporary memory usage.

The overall Step2 logic remains unchanged, but is now faster and more efficient with a redesigned comparison pipeline.

---

# 1. Changes

The following components were altered:

- Introduced new class: `LuksIsomorphism`
- Modified `compare_chem_envs`
- Modified `unique_adsorbate`
- Reimplemented `unique_chem_envs`

---

# 2. Isomorphism Acceleration for `compare_chem_envs` and `unique_adsorbate`

## Cheap Invariants (Early Rejection)

For `compare_chem_envs` and `unique_adsorbate`, we added fast structural filters before exact matching:

- Node count
- Edge count
- Degree sequence
- Multiset of element labels

These invariants eliminate most non-isomorphic cases before invoking expensive matching.

## VF2 → bliss Replacement

Exact isomorphism now uses:

- NetworkX → igraph conversion
- `"element"` encoded as vertex color
- `igraph.Graph.isomorphic_bliss()`

This replaces NetworkX VF2 and significantly improves eifficiency for chemical graphs.

---

# 3. Redesign of `unique_chem_envs`

### Previous Logic

The original implementation used deeply nested loops:

- Double loop over configurations
- Inside each comparison, double loop over graphs within configurations
- Each graph pair compared via isomorphism

Additionally, full graph objects were cached in memory during grouping.

---

### New Logic

The new backend replaces pairwise comparison with canonical labeling and hashing.

For each graph:

1. Perform canonical permutation (via bliss)
2. Get canonical labeling form as a string
3. Hash the string to bytes (BLAKE2b)

For each configuration:

1. Collect all ego-graph hash values
2. Sort them
3. Stack into a tuple
4. Use this tuple as the configuration key

Grouping is then performed using:

```python
from collections import defaultdict
```

Each configuration is inserted into a dictionary based on its canonical key.

---

### Complexity and Memory Improvement

#### Previous Complexity

Let:

- C = number of configurations
- G = number of graphs per configuration

The previous implementation required:

- Pairwise comparison between configurations → O(C²)
- Inside each comparison, pairwise graph comparison → O(G²)

Total worst-case complexity:

O(C² × G²)

This leads to poor efficiency as both C and G grow.

#### New Complexity

For each configuration:

- Canonical labeling is computed once per graph
- Hash keys are constructed and sorted once per configuration
- Dictionary bucketing is performed once for each configuration

Overall complexity:

O(C × G)

After canonical permutation, grouping is linear in the number of configurations and adsorbates.

---

# 4. Usage

## Isomorphism Check For Two Graphs

```python
from chemical_environment_new import LuksIsomorphism

luks = LuksIsomorphism()
result = luks.is_isomorphic(G1, G2)
```

## Compare Two Configurations Directly

```python
from chemical_environment_new import compare_chem_envs

# config1 and config2 are lists of NetworkX graphs
result = compare_chem_envs(config1, config2)

print(result)  # True if configurations are isomorphic
```


## Environment Grouping for a List of Configurations

```python
from chemical_environment_new import unique_chem_envs

groups = unique_chem_envs(chem_envs_groups, metadata=paths)
```

---

# 5. Dependencies

- networkx
- numpy
- ase
- python-igraph

Recommended installation:

```
conda install -c conda-forge python-igraph
```

---

# Summary

This repository modernizes the backend of chemical environment comparison by:

- Introducing bliss-based colored isomorphism and canonical labeling concepts
- Adding cheap structural invariants
- Replacing nested pairwise comparison with canonical hash grouping
- Reducing complexity 
- Reducing memory usage by replacing graph object storage with byte hashes

The external workflow remains unchanged while achieving speed and memory improvements.

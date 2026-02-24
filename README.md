# chem-env-isomorphism-fast

Memory-efficient chemical environment isomorphism backend (igraph + bliss)

---

# Fast Chemical Environment Isomorphism Backend

Accelerated and memory-efficient graph isomorphism backend for chemical environment grouping workflows.

---

## Overview

This repository provides an optimized backend implementation for chemical graph isomorphism and environment deduplication used in a larger catalysis workflow.

The workflow-level Step2 logic remains structurally unchanged.  
Performance improvements come from replacing the underlying isomorphism and comparison backend.

Target graph type:

- Undirected graphs
- Nodes represent atoms
- Edges represent bonds
- Node attribute `"element"` defines atom type
- Bounded valence chemical graphs

---

## Context

Originally, environment grouping (Step2) relied on:

- NetworkX VF2-based isomorphism
- Nested pairwise comparison across configurations
- Graph object storage during grouping

This repository replaces that backend while keeping the external workflow structure intact.

Step2 continues to call:

```python
unique_chem_envs(chem_envs_groups, metadata=paths)
```

but now benefits from a redesigned comparison engine.

---

# 1. Backend Structural Changes

The following components were redesigned:

- Introduced new class: `LuksIsomorphism`
- Modified `compare_chem_envs`
- Modified `unique_adsorbate`
- Reimplemented backend of `unique_chem_envs`

---

# 2. Isomorphism Acceleration

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

This replaces NetworkX VF2 and significantly improves performance for bounded-valence chemical graphs.

---

# 3. Redesign of `unique_chem_envs`

### Previous Logic

The original implementation used deeply nested loops:

- Double loop over configurations
- Inside each comparison, double loop over graphs within configurations
- Each graph pair compared via isomorphism

This results in approximately:

O(N⁴) behavior in worst case  
(two configuration loops × two graph loops)

Additionally, full graph objects were retained in memory during grouping.

---

### New Logic

The redesigned backend replaces pairwise comparison with canonical labeling and hashing.

For each graph:

1. Perform canonical labeling (via bliss)
2. Represent canonical form as a string
3. Hash the string to fixed-length bytes (BLAKE2b)

For each configuration:

1. Collect all graph hash values
2. Sort them
3. Stack into a tuple
4. Use this tuple as the configuration key

Grouping is then performed using:

```python
from collections import defaultdict
```

Each configuration is inserted into a dictionary bucket based on its canonical key.

---

### Complexity and Memory Improvement

## Previous Complexity

Let:

- C = number of configurations
- G = average number of graphs per configuration

The previous implementation required:

- Pairwise comparison between configurations → O(C²)
- Inside each comparison, pairwise graph comparison → O(G²)

Total worst-case complexity:

O(C² × G² × T_iso)

where T_iso is the cost of a single graph isomorphism.

This leads to poor scalability as both C and G grow.

### New Complexity

For each configuration:

- Canonical labeling is computed once per graph
- Hash keys are constructed and sorted
- Dictionary bucketing is performed

Overall complexity:

O(C × G × T_label)

After canonical labeling, grouping is linear in the number of configurations.

---

# 4. Effect on Step2 Workflow

The Step2 workflow:

- Still constructs environment groups
- Still calls `unique_chem_envs`
- Keeps the same control structure

Performance gains come entirely from the redesigned backend.

Minor improvements were also introduced in Step2:

- Reduced temporary memory usage
- Stream-based file reading instead of full `readlines()` loading

However, the overall Step2 logic remains unchanged.

---

# 5. Usage

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

# 6. Assumptions

- Undirected graphs
- No self-loops
- No multi-edges
- Element-preserving mapping
- Bounded valence

---

# 7. Dependencies

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

- Introducing bliss-based colored isomorphism
- Adding cheap structural invariants
- Replacing nested pairwise comparison with canonical hash grouping
- Reducing complexity from ~O(N⁴) to ~O(N)
- Reducing memory usage by replacing graph object storage with byte hashes

The external workflow remains unchanged while achieving substantial speed and memory improvements.

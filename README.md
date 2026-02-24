# chem-env-isomorphism-fast
Memory-efficient chemical environment isomorphism checker (igraph + bliss)
# Fast Chemical Environment Isomorphism Backend

Accelerated and memory-efficient graph isomorphism backend for chemical environment grouping workflows.

---

## Overview

This repository provides an optimized backend implementation for chemical graph isomorphism checking used in environment grouping workflows.

The workflow-level Step2 code remains largely unchanged.  
Performance improvement comes from replacing the underlying isomorphism and environment comparison logic.

Target graph type:

- Undirected graphs
- Nodes represent atoms
- Edges represent bonds
- Node attribute `"element"` defines atom type
- Bounded valence chemical graphs

---
## Context

This repository contains an optimized backend extracted from an existing chemical environment grouping workflow.

Originally, environment grouping (Step2) was implemented inside a larger catalysis workflow and relied on:

- NetworkX-based isomorphism (VF2)
- Object-level pairwise environment comparison

This repository isolates and replaces that backend with:

- bliss-based colored graph isomorphism (via igraph)
- canonical hash-based environment comparison

The external Step2 workflow code remains structurally unchanged and continues to call:

```python
unique_chem_envs(chem_envs_groups, metadata=paths)
```

Thus, this repository serves as a drop-in backend replacement for accelerating environment comparison without modifying workflow-level logic.

---


# 1. Isomorphism Backend

### Class

`LuksIsomorphism`

### Purpose

Provide fast element-preserving graph isomorphism checking.

### Pipeline

#### Stage 1 — Cheap Invariants

- Node count
- Edge count
- Degree sequence
- Multiset of element labels

Early rejection avoids expensive checks.

#### Stage 2 — Exact Colored Isomorphism

- Convert NetworkX → igraph
- Encode `"element"` as vertex color
- Use `igraph.Graph.isomorphic_bliss()`

This replaces NetworkX VF2-based isomorphism.

---

# 2. Canonical Hash-Based Environment Comparison

The main structural improvement is replacing object-based comparison with canonical key comparison.

Instead of storing graph objects and performing repeated pairwise comparisons:

- Each graph is converted to a canonical representation
- Canonical representation is hashed (BLAKE2b)
- Environment groups are represented by sorted tuples of graph hash keys
- Grouping is performed using dictionary bucketing

---

# 3. Effect on Step2 Workflow

The external Step2 workflow code:

- Still calls `unique_chem_envs`
- Still groups chemical environments
- Logic structure unchanged

Performance improvement comes from:

- Faster isomorphism backend
- Elimination of repeated pairwise comparisons
- Reduced memory usage (no graph objects retained)

Thus Step2 runtime and memory footprint are significantly reduced without altering workflow-level logic.

---

# 4. Usage

## Isomorphism Check

```python
from chemical_environment_new import LuksIsomorphism

luks = LuksIsomorphism()
result = luks.is_isomorphic(G1, G2)
```

## Environment Grouping (Called by Step2)

```python
from chemical_environment_new import unique_chem_envs

groups = unique_chem_envs(chem_envs_groups, metadata=paths)
```

---

# 5. Assumptions

- Undirected graphs
- No self-loops
- No multi-edges
- Element-preserving mapping
- Bounded valence

---

# 6. Dependencies

- networkx
- numpy
- ase
- python-igraph

Recommended installation:

```
conda install -c conda-forge python-igraph
```

---

# 7. Changes compared to the original workflow

- Replaced NetworkX VF2 isomorphism with igraph/bliss colored isomorphism.
- Added cheap invariants for early rejection.
- Added canonical digest keys (BLAKE2b) to represent graphs/environments.
- Minor Step2 memory tweaks.

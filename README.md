# chem-env-isomorphism-fast
Memory-efficient chemical environment isomorphism checker (igraph + bliss)
# Fast Chemical Environment Isomorphism & Deduplication

Memory-efficient and accelerated graph isomorphism and environment grouping for heterogeneous catalysis workflows.

---

## Overview

This repository provides a scalable implementation for:

- Fast chemical graph isomorphism checking  
- Environment grouping (Step2 in workflow)  
- Reduced memory usage for large configuration sets  

Target graph type:

- Undirected graphs  
- Nodes represent atoms  
- Edges represent bonds  
- Each node has `"element"` attribute (e.g. `"Pt"`, `"C"`, `"O"`, `"H"`)  
- Bounded valence (chemical graphs)

---

# 1. Isomorphism Checker

### Class

`LuksIsomorphism`

### Input

Two NetworkX graphs that:

- Are undirected  
- Contain node attribute `"element"`  

### Output

Boolean value indicating element-preserving isomorphism.

---

## Algorithm

The checker uses a two-stage pipeline:

### Stage 1 — Cheap Invariants (Early Rejection)

Before running exact isomorphism:

- Compare number of nodes  
- Compare number of edges  
- Compare sorted degree sequence  
- Compare sorted multiset of element labels  

If any mismatch occurs, graphs are rejected immediately.

---

### Stage 2 — Exact Colored Isomorphism (bliss backend)

If invariants match:

1. Convert NetworkX → igraph  
2. Encode element types as vertex colors  
3. Call:

```python
igraph.Graph.isomorphic_bliss()
```

This approach is significantly faster than NetworkX VF2 for bounded-valence chemical graphs.

---

# 2. Step2 Optimization (Environment Deduplication)

## Previous Implementation

Old `unique_chem_envs`:

- Stored full graph objects  
- Compared new group against all stored groups  
- Used pairwise isomorphism checks  

Result:

- O(N²) comparisons  
- High memory usage  
- Poor scalability  

---

## New Implementation — Canonical Key Bucketing

Instead of storing graph objects, we compute canonical hash keys.

### For each graph:

1. Convert NetworkX → igraph  
2. Compute canonical permutation via bliss  
3. Generate canonical representation (colors + edges)  
4. Hash using BLAKE2b  

### For each environment group:

- Build sorted tuple of graph hash keys  
- Use tuple as dictionary key  
- Bucket indices directly  

---

## Why This Is Faster

- No graph objects stored after key generation  
- No O(N²) pairwise comparisons  
- Dictionary-based grouping (~O(N))  
- Fixed-size byte hashes reduce memory footprint  

---

# 3. Usage

## A. Isomorphism Check

```python
from chemical_environment_new import LuksIsomorphism

luks = LuksIsomorphism()
result = luks.is_isomorphic(G1, G2)

print(result)
```

`G1` and `G2` must be NetworkX graphs with `"element"` node attributes.

---

## B. Environment Deduplication (Step2)

```python
from chemical_environment_new import unique_chem_envs

groups = unique_chem_envs(chem_envs_groups, metadata=paths)

for bucket in groups:
    print("Duplicate group:")
    print(bucket)
```

Return value:

- List of buckets  
- Each bucket contains metadata entries for isomorphic groups  

---

# 4. Assumptions

- Undirected graphs  
- No self-loops  
- No multi-edges  
- Element-preserving mapping  
- Bounded valence  

---

# 5. Dependencies

Python ≥ 3.9

Required packages:

```
networkx
numpy
ase
python-igraph
```

Recommended igraph installation:

```
conda install -c conda-forge python-igraph
```

---

# Summary

This implementation improves:

- Isomorphism speed via bliss backend  
- Step2 scalability via canonical hashing  
- Memory efficiency by avoiding graph object storage  

Designed for large-scale adsorbate configuration enumeration workflows.

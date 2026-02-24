# chem-env-isomorphism-fast
Memory-efficient chemical environment isomorphism checker (igraph + bliss)
Fast Chemical Environment Isomorphism & Deduplication
This repository provides a memory-efficient and accelerated implementation for chemical environment graph isomorphism checking and environment deduplication in heterogeneous catalysis workflows.
The implementation targets bounded-valence, undirected chemical graphs where:
Nodes represent atoms
Edges represent bonds
Each node carries an "element" attribute (e.g., "Pt", "C", "O", "H")
The goal is to accelerate:
Graph isomorphism checking
Environment grouping (Step2 in the workflow)
Reduce memory footprint during large-scale configuration enumeration
1. Isomorphism Checker
The LuksIsomorphism class provides a fast colored graph isomorphism check.
Input
Two NetworkX graphs that:
Are undirected
Contain node attribute "element"
Output
Boolean value indicating whether the two graphs are isomorphic under element-preserving mapping.
Algorithm Overview
The checker uses a two-stage pipeline:
Stage 1 — Cheap Invariants (Early Rejection)
Before running an exact isomorphism routine, the following invariants are compared:
Number of nodes
Number of edges
Sorted degree sequence
Sorted multiset of node element labels
If any invariant differs, the graphs are rejected immediately.
Stage 2 — Exact Colored Isomorphism (bliss backend)
If invariants match:
Convert NetworkX graphs to igraph format
Encode atom element types as vertex colors
Use igraph.Graph.isomorphic_bliss() for exact colored graph isomorphism
This approach is significantly faster and more scalable than NetworkX’s VF2-based is_isomorphic for bounded-valence chemical graphs.
2. Step2 Optimization (Environment Deduplication)
Previous Implementation
The previous unique_chem_envs implementation:
Stored full graph objects for each unique environment group
Compared each new group against all stored groups
Used pairwise isomorphism checking
This resulted in:
O(N²) comparisons in worst case
High memory usage
Poor scalability for large configuration sets
New Implementation (Canonical Key Bucketing)
The optimized version replaces object-based comparison with canonical hashing.
For each graph:
Convert NetworkX → igraph
Compute canonical permutation using bliss
Generate canonical representation (colors + edge list)
Hash the canonical representation using BLAKE2b
For each environment group:
Compute a sorted tuple of graph hash keys
Use this tuple as dictionary key
Bucket indices directly
Benefits
No graph objects stored after key computation
No O(N²) pairwise comparisons
Dictionary-based grouping (~O(N))
Fixed-size hash keys reduce memory usage
3. Usage
A. Isomorphism Check
from chemical_environment_new import LuksIsomorphism

luks = LuksIsomorphism()
result = luks.is_isomorphic(G1, G2)

print(result)
G1 and G2 must be NetworkX graphs with "element" node attributes.
B. Environment Deduplication (Step2)
from chemical_environment_new import unique_chem_envs

# chem_envs_groups: list of environment groups
# each group is a list of NetworkX graphs
# metadata: optional list aligned with chem_envs_groups (e.g., file paths)

groups = unique_chem_envs(chem_envs_groups, metadata=paths)

for bucket in groups:
    print("Duplicate group:")
    print(bucket)
Return value:
List of buckets
Each bucket contains metadata entries corresponding to isomorphic environment groups
4. Assumptions
Undirected graphs
No self-loops
No multi-edges
Bounded valence (chemical graphs)
Element-preserving isomorphism
5. Dependencies
Python ≥ 3.9
Required packages:
networkx
numpy
ase
python-igraph
Recommended igraph installation:
conda install -c conda-forge python-igraph
Summary
This implementation improves:
Isomorphism check performance via bliss backend
Step2 scalability via canonical hashing
Memory efficiency by avoiding storage of graph objects
It is designed for large-scale adsorbate configuration enumeration workflows.

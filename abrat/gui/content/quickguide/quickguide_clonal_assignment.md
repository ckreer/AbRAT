### Clonal Assignment

{{AbRAT}} features a modular clonal assignment framework. Users can choose which chain(s) (e.g., heavy and/or light) 
to include in the clustering process and select which V(D)J gene segment information to use for pre-grouping. 
Each selected chain is clustered independently, and a global cluster is formed by combining the resulting subclusters. 
Rows with missing or non-computable values are assigned to subcluster 0.

The framework currently offers three clustering algorithms for CDR3 sequences (at the amino acid level):
- **Iterative Greedy CDR3 Clustering (Original Approach)**
- **Matrix-Based Greedy CDR3 Clustering (Global Greedy Approach)**
- **Hierarchical CDR3 Clustering**

For all methods, common preprocessing steps are applied: filtering out missing data, applying a configurable threshold 
on sequence length differences, and using a normalized Levenshtein distance threshold. Moreover, each algorithm can be run 
multiple times with different random initializations to minimize the number of unassigned sequences ("singles"), optimizing 
the final clustering result.
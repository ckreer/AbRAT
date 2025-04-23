.. AbRAT Clustering algorithm documentation file
.. include:: shared.rst

Clustering Algorithms
=====================

Iterative Greedy CDR3 Clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the original algorithm that was used in the cited publications in the Quick Guide.

**How It Works**

- Valid CDR3 sequences are processed sequentially:
  - The first sequence in the list starts a new cluster.
  - Each subsequent candidate is compared first with the cluster’s representative using a length filter (ensuring the difference in length does not exceed a specified threshold) and then against **all** members already in the cluster using the normalized Levenshtein distance.
  - Only if the candidate meets both criteria for every cluster member is it added to the cluster.
- The process repeats with the remaining unclustered sequences.
- Multiple iterations with different random orderings can be performed to minimize the number of singles.

**Pros**

- **High Homogeneity:** Every candidate is verified against all existing cluster members, ensuring very homogeneous clusters.
- **Historical Consistency:** This is the original approach used in previous versions.

**Cons**

- **Performance:** Sequential pairwise comparisons can be slow for large datasets.
- **Order Sensitivity:** The final clustering result may depend on the initial sequence order.

----

Matrix-Based Greedy CDR3 Clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**How It Works**

- A full pairwise similarity matrix is computed for all valid CDR3 sequences using both the normalized Levenshtein distance and a length filter.
- A global greedy strategy is applied:
  - In each iteration, the sequence with the highest overall similarity (i.e., the one with the most matches according to the filters) is chosen as the cluster representative.
  - All sequences similar to this representative (meeting both length and distance criteria) are grouped together.
  - These sequences are then removed from further consideration, and the process repeats.
- As with the iterative method, multiple iterations can be performed to reduce the number of singles.

**Pros**

- **Efficiency:** Vectorized computation of the similarity matrix makes it faster for larger datasets.
- **Global Perspective:** Quickly identifies clusters with many similar sequences based on overall similarity counts.

**Cons**

- **Potential Heterogeneity:** Because clustering is based on overall counts rather than verifying each candidate against every cluster member, some pairs within a cluster might slightly exceed the threshold.
- **Less Fine-Grained:** New candidates may not be compared as rigorously on a pairwise basis compared to the iterative approach.

----

Hierarchical CDR3 Clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**How It Works**

- For valid sequences, a full pairwise distance matrix is computed that combines the normalized Levenshtein distance with a length filter—this ensures that only sequences with acceptable length differences are compared.
- The resulting distance information is then condensed to focus on the unique distances between sequence pairs.
- Using this condensed distance data, a hierarchical clustering procedure is applied, where the most similar sequences are iteratively merged into clusters to form a dendrogram representing the overall structure.
- Finally, a threshold is applied to cut the dendrogram into flat clusters.
- Multiple iterations can also be run to optimize the result, though hierarchical clustering is deterministic and less order-sensitive than the greedy approaches.

**Pros**

- **Flexible Clustering Structure:** Hierarchical clustering produces a dendrogram, which can be "cut" at different levels, providing a flexible view of the data structure.
- **No Need to Predefine Cluster Number:** A threshold is used to form clusters, which can be particularly useful for heterogeneous datasets.

**Cons**

- **Computational Complexity:** Calculating the full distance matrix and linkage can be computationally intensive for very large datasets.
- **Sensitivity to Cutoff:** The choice of the cutoff threshold is critical and may require fine-tuning.

----

Summary of Differences
^^^^^^^^^^^^^^^^^^^^^^

- The **Iterative Greedy CDR3 Clustering (Original Approach)** builds clusters by sequentially verifying each candidate against all current cluster members. It produces very homogeneous clusters but is slower and more order-sensitive.
- The **Matrix-Based Greedy CDR3 Clustering (Global Greedy Approach)** uses a full similarity matrix to quickly form clusters based on overall similarity counts, trading off some pairwise strictness for efficiency.
- The **Hierarchical CDR3 Clustering** method computes a full distance matrix and uses a dendrogram-based approach to form clusters. It offers flexibility in selecting the clustering threshold but can be computationally heavy and sensitive to the chosen cutoff.

All methods share common preprocessing steps—such as filtering out missing data, applying a configurable length
difference threshold, and using a normalized Levenshtein distance threshold—and they all allow for multiple
iterations with random initializations to minimize the number of unassigned sequences ("singles").

Note that iterations are generally critical for the original (iterative) greedy algorithm because its outcome
is highly dependent on the order in which sequences are processed, leading to different clustering outcomes when the
order is changed. Running multiple iterations with different random orderings helps overcome local optima and
reduces the chance that a suboptimal ordering leads to too many unassigned sequences ("singles"). In contrast, the
matrix-based approach computes a full pairwise similarity matrix and uses a global view to choose the representative
with the highest overall similarity. Although randomization can still affect the tie-breaking or the grouping in
borderline cases, the global nature of the matrix reduces the sensitivity to ordering. Therefore, iterations might
have a less pronounced effect on the matrix-based algorithm. Hierarchical clustering is deterministic and largely
insensitive to input order; however, in edge cases with many nearly identical distances leading to tie-breaking
ambiguities, additional iterations with random variations may help.
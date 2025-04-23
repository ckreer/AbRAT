### JACCARD INDEX

# in app
# TODO: Calculate overlap indices and save in page state st.session_state.page_states[page_name]:
jaccard_results, jaccard_matrix = compute_pairwise_jaccard(data_dict=series_dict,
                                                           iterations = 20,
                                                           random_state = 42)
st.session_state.page_states[page_name]['results'][chain][mode]['jaccard_results'] = jaccard_results
st.session_state.page_states[page_name]['results'][chain][mode]['jaccard_matrix'] = jaccard_matrix


with st.container():
    st.markdown('''#### CDR3 Similarities''')
    index_dict = {'jaccard': 'Jaccard Index'}
    diversity_index = st.radio("Overlap Index",
                               ['jaccard'],
                               format_func=lambda x: index_dict.get(x),
                               horizontal=True
                               )
    cols = st.columns([1, 1, 1])
    for i, cdr3_chain in enumerate(selected_chain_ordered):
        with cols[i]:
            df = st.session_state.page_states[page_name]['results'][cdr3_chain][cdr3_mode][
                diversity_index + '_matrix'].copy()
            if any([col in alias_dict for col in df.columns]):
                df.columns = [alias_dict[col] for col in df.columns]
            if any([row in alias_dict for row in df.index]):
                df.index = [alias_dict[row] for row in df.index]
            st.write(f"**{chain_dict[cdr3_chain]}**")
            jaccard_fig = create_heatmap_with_annotations(df=df,
                                                          colorscale='blues',
                                                          margin_color='black',
                                                          margin_width=0.5,
                                                          color_bar_title=index_dict[diversity_index]
                                                          )
            st.plotly_chart(jaccard_fig, use_container_width=False, key=cdr3_chain + "jaccardPlot")

st.markdown(r"""
    #### CDR3 similarity:
    :blue[_AbRAT_] currently implements two functions to assess the similarity between different repertoires. 

    The **Jaccard Index** is defined as:
    $$
    J(A, B) = \frac{|A \cap B|}{|A \cup B|}
    $$
    where $A$ and $B$ are the sets of unique CDR3 sequences from two repertoires. In our workflow, each dataset is first 
    reduced to eliminate duplicate sequences, ensuring that each sequence is only counted once. Randomization is then 
    applied by subsampling each dataset to the size of the smallest repertoire (sampling without replacement) over 
    multiple iterations. This process minimizes the bias introduced by varying dataset sizes and provides robust average 
    similarity values. 

    The final similarity matrix contains values ranging from $0$ to $1$. A value of **$0$** indicates no similarity (no 
    shared or similar CDR3 sequences), while a value of **$1$** indicates complete overlap between the repertoires.

    _**Note:** This measures only accounts for **exact** CDR3 matches between repertoires._     
    """)


### ANLEITUNG für weighted Jaccard
st.markdown(r"""
    #### CDR3 similarity:
    :blue[_AbRAT_] currently implements two functions to assess the similarity between different repertoires.

    The **Jaccard Index** is defined as:
    $$
    J(A, B) = \frac{|A \cap B|}{|A \cup B|}
    $$
    where $A$ and $B$ are the sets of unique CDR3 sequences from two repertoires. In our workflow, each dataset is first
    reduced to eliminate duplicate sequences, ensuring that each sequence is only counted once. Randomization is then
    applied by subsampling each dataset to the size of the smallest repertoire (sampling without replacement) over
    multiple iterations. This process minimizes the bias introduced by varying dataset sizes and provides robust average
    similarity values.


    The **Weighted Jaccard Index** extends the classical Jaccard Index by considering partial matches between sequences.
    It is computed as:
    $$
    J_w(A, B) = \frac{\text{soft\_intersection}(A, B)}{|A| + |B| - \text{soft\_intersection}(A, B)}
    $$
    where **soft\_intersection** is calculated as the sum over each sequence $a$ in set $A$ of the maximum similarity
    score to any sequence $b$ in set $B$. The similarity between two sequences is typically measured with Python’s `SequenceMatcher`
     that returns a score between $0$ (no similarity) and $1$ (identical).

    **_Weighted Jaccard Index Example:_**

    Consider two example CDR3 sequences:
    - **Seq1:** `"CASSLGTDTQYF"`
    - **Seq2:** `"CASRLGNDTQYF"`

    The algorithm identifies the following matching blocks:
    - `"DTQYF"` with a length of $5$,
    - `"CAS"` with a length of $3$,
    - `"LG"` with a length of $2$.

    The total matching length is $5 + 3 + 2 = 10$. With both sequences having a length of $12$, the weighted similarity
    ratio is calculated as:
    $$
    \text{ratio} = \frac{2 \times 10}{12 + 12} \approx 0.83
    $$
    This score represents the degree of similarity between the two sequences, taking into account both exact and partial
    matches.

    The final similarity matrix contains values ranging from $0$ to $1$. A value of **$0$** indicates no similarity (no
    shared or similar CDR3 sequences), while a value of **$1$** indicates complete overlap between the repertoires.

    """)


# Anleitung Clustering Algorithmen:
st.markdown("""
---
    #### Iterative Greedy CDR3 Clustering (Original Approach)

    **How It Works:**  
    - Valid CDR3 sequences are processed sequentially:  
      - The first sequence in the list starts a new cluster.  
      - Each subsequent candidate is compared first with the cluster’s representative using a length filter (ensuring 
      the difference in length does not exceed a specified threshold) and then against **all** members already in the 
      cluster using the normalized Levenshtein distance.  
      - Only if the candidate meets both criteria for every cluster member is it added to the cluster.  
    - The process repeats with the remaining unclustered sequences.
    - Multiple iterations with different random orderings can be performed to minimize the number of singles.

    **Pros:**  
    - **High Homogeneity:** Every candidate is verified against all existing cluster members, ensuring very homogeneous 
    clusters.  
    - **Historical Consistency:** This is the original approach used in previous versions.

    **Cons:**  
    - **Performance:** Sequential pairwise comparisons can be slow for large datasets.  
    - **Order Sensitivity:** The final clustering result may depend on the initial sequence order.

    ---
    #### Matrix-Based Greedy CDR3 Clustering (Global Greedy Approach)

    **How It Works:**  
    - A full pairwise similarity matrix is computed for all valid CDR3 sequences using both the normalized Levenshtein 
    distance and a length filter.  
    - A global greedy strategy is applied:  
      - In each iteration, the sequence with the highest overall similarity (i.e., the one with the most matches 
      according to the filters) is chosen as the cluster representative.  
      - All sequences similar to this representative (meeting both length and distance criteria) are grouped together.  
      - These sequences are then removed from further consideration, and the process repeats.
    - As with the iterative method, multiple iterations can be performed to reduce the number of singles.

    **Pros:**  
    - **Efficiency:** Vectorized computation of the similarity matrix makes it faster for larger datasets.  
    - **Global Perspective:** Quickly identifies clusters with many similar sequences based on overall similarity counts.

    **Cons:**  
    - **Potential Heterogeneity:** Because clustering is based on overall counts rather than verifying each candidate 
    against every cluster member, some pairs within a cluster might slightly exceed the threshold.  
    - **Less Fine-Grained:** New candidates may not be compared as rigorously on a pairwise basis compared to the 
    iterative approach.

    ---

    #### Hierarchical CDR3 Clustering

    **How It Works:**  
    - For valid sequences, a full pairwise distance matrix is computed that combines the normalized Levenshtein distance 
    with a length filter—this ensures that only sequences with acceptable length differences are compared. 
    - The resulting distance information is then condensed to focus on the unique distances between sequence pairs. 
    - Using this condensed distance data, a hierarchical clustering procedure is applied, where the most similar 
    sequences are iteratively merged into clusters to form a dendrogram representing the overall structure. 
    - Finally, a threshold is applied to cut the dendrogram into flat clusters. 
    - Multiple iterations can also be run to optimize the result, though hierarchical clustering is deterministic and 
    less order-sensitive than the greedy approaches.

    **Pros:**  
    - **Flexible Clustering Structure:** Hierarchical clustering produces a dendrogram, which can be "cut" at different 
    levels, providing a flexible view of the data structure.  
    - **No Need to Predefine Cluster Number:** A threshold is used to form clusters, which can be particularly useful 
    for heterogeneous datasets.

    **Cons:**  
    - **Computational Complexity:** Calculating the full distance matrix and linkage can be computationally intensive 
    for very large datasets.  
    - **Sensitivity to Cutoff:** The choice of the cutoff threshold is critical and may require fine-tuning.

    ---
    #### Summary of Differences

    - The **Iterative Greedy CDR3 Clustering (Original Approach)** builds clusters by sequentially verifying each 
    candidate against all current cluster members. It produces very homogeneous clusters but is slower and more order-sensitive.  
    - The **Matrix-Based Greedy CDR3 Clustering (Global Greedy Approach)** uses a full similarity matrix to quickly 
    form clusters based on overall similarity counts, trading off some pairwise strictness for efficiency.  
    - The **Hierarchical CDR3 Clustering** method computes a full distance matrix and uses a dendrogram-based approach 
    to form clusters. It offers flexibility in selecting the clustering threshold but can be computationally heavy and 
    sensitive to the chosen cutoff.

    All methods share common preprocessing steps—such as filtering out missing data, applying a configurable length 
    difference threshold, and using a normalized Levenshtein distance threshold—and they all allow for multiple 
    iterations with random initializations to minimize the number of unassigned sequences ("singles"). 

    Note that iterations are generally critical for the original (iterative) greedy algorithm because its outcome 
    is highly dependent on the order in which sequences are processed, leading to different clustering outcomes when the 
    order is changed. Running multiple iterations with different random orderings helps overcome local optima and 
    reduces the chance that a suboptimal ordering leads to too many unassigned sequences ("singles").
    In contrast, the matrix-based approach computes a full pairwise similarity matrix and uses a global view to choose 
    the representative with the highest overall similarity. Although randomization can still affect the tie-breaking 
    or the grouping in borderline cases, the global nature of the matrix reduces the sensitivity to ordering. 
    Therefore, iterations might have a less pronounced effect on the matrix-based algorithm. Hierarchical clustering 
    is deterministic and largely insensitive to input order; however, in edge cases with many nearly identical 
    distances leading to tie-breaking ambiguities, additional iterations with random variations may help.

    """)

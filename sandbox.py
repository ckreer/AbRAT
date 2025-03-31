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
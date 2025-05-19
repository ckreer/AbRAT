#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib as mpl
import random
import re
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from rapidfuzz import distance

from abrat.core.utils import plot_pie_chart

# ============
# Column names
# =============
# can be adapted, if changed later, e.g. to airr style
sample_col = 'SAMPLE_INFORMATION'
hc_col = 'HEAVY_CHAIN'
lc_col = 'LIGHT_CHAIN'

v_gene_subcol = 'V_GENE'
d_gene_subcol = 'D_GENE'
j_gene_subcol = 'J_GENE'

top_v_gene_subcol = 'TOP_V'
top_d_gene_subcol = 'TOP_V'
top_j_gene_subcol = 'TOP_J'

cdr3_aa_subcol = 'CDR3_AA'
cdr3_nt_subcol = 'CDR3_NT'

hc_cluster_subcol = 'HC_SUBCLUSTER'
lc_cluster_subcol = 'LC_SUBCLUSTER'
representative_bcr_col = 'CLUSTER_REPRESENTATIVE'
hclc_cluster_subcol = 'HC-LC_CLUSTER'
cluster_size_subcol = 'CLUSTER_SIZE'
clonal_subcol = 'IS_CLONAL'
clone_subcol = 'CLONE'
clone_color_subcol = 'CLONE_COLOR'

clone_color_dict = {'palette': 'GnBu',
                    'non_clonal': '#D3D3D3',  # light gray
                    'undefined': '#FFFFFF'    # white
                   }

# ==========
# Functions
# ==========

def parse_clone_number(clone_label):
    """
    Parses the clone label and extracts the numeric clone number.

    Parameters:
        clone_label (str): The clone label in the format "Clone-<number>".

    Returns:
        int or float: The numeric clone number if found; otherwise, infinity.
    """
    match = re.match(r'^Clone-(\d+)$', clone_label)
    if match:
        return int(match.group(1))
    else:
        return float('inf')


def combine_clusters(row):
    """
    Generates a combined cluster ID based on the availability of heavy and light chain sub-clusters.

    Parameters:
        row (pd.Series): A DataFrame row containing the heavy chain and light chain cluster values.

    Returns:
        str or None: The combined cluster ID in the format "HC-LC" if available; otherwise, None.
    """
    hc_val = row[(hc_col, hc_cluster_subcol)]
    lc_val = row[(lc_col, lc_cluster_subcol)]
    # Here, hc_val and lc_val are already integers.
    if hc_val != 0 and lc_val != 0:
        return f"{hc_val}-{lc_val}"
    elif hc_val != 0:
        return f"{hc_val}-0"
    elif lc_val != 0:
        return f"0-{lc_val}"
    else:
        return None


def process_bcr_list(df, settings):
    """
    Processes a BCR list DataFrame to perform grouping, clustering, and clonal assignment based on provided settings.

    This function standardizes missing values, performs clustering on heavy and light chains (if selected),
    combines the clusters into a composite cluster ID, determines cluster sizes and clonality, assigns clone names,
    and optionally applies coloring to clones.

    Parameters:
        df (pd.DataFrame): The BCR list DataFrame.
        settings (dict): Dictionary with processing settings (e.g., heavy_chain, light_chain, autocolor).

    Returns:
        pd.DataFrame: The processed DataFrame with cluster IDs, representative BCRs, clone names, and additional clone information.
    """
    # Standardize missing values
    df = standardize_missing_values(df)

    # Heavy-chain clustering if selected
    if settings.get('heavy_chain', False) and hc_col in df.columns.levels[0]:
        heavy_df = df[hc_col]
        clusters, reps = cluster_bcr_chains(heavy_df, settings.get('heavy_chain', {}))
        df[(hc_col, hc_cluster_subcol)] = clusters
        df[(hc_col, representative_bcr_col)] = reps
    else:
        df[(hc_col, hc_cluster_subcol)] = 0
        df[(hc_col, representative_bcr_col)] = None

    # Light-chain clustering if selected
    if settings.get('light_chain', False) and lc_col in df.columns.levels[0]:
        light_df = df[lc_col]
        clusters, reps = cluster_bcr_chains(light_df, settings.get('light_chain', {}))
        df[(lc_col, lc_cluster_subcol)] = clusters
        df[(lc_col, representative_bcr_col)] = reps
    else:
        df[(lc_col, lc_cluster_subcol)] = 0
        df[(lc_col, representative_bcr_col)] = None

    # Convert cluster IDs to integer (fill NaNs with 0)
    df[(hc_col, hc_cluster_subcol)] = df[(hc_col, hc_cluster_subcol)].fillna(0).astype(int)
    df[(lc_col, lc_cluster_subcol)] = df[(lc_col, lc_cluster_subcol)].fillna(0).astype(int)

    # Combine heavy and light chain clusters into a composite cluster
    df[(sample_col, hclc_cluster_subcol)] = df.apply(combine_clusters, axis=1)
    # Determine cluster sizes and clonality
    df[(sample_col, cluster_size_subcol)] = df.groupby([(sample_col, hclc_cluster_subcol)])[[(sample_col, hclc_cluster_subcol)]].transform('count')
    df[(sample_col, clonal_subcol)] = df[(sample_col, cluster_size_subcol)] > 1

    df = assign_clone_names(df)
    if settings.get('autocolor', True):
        df = color_clones(df, clone_color_dict)
    return df


def cluster_bcr_chains(chain_df, settings):
    """
    Clusters a chain DataFrame based on provided settings.

    If settings['group_by'] is specified, the DataFrame is split into subgroups.
    Within each subgroup, clustering is applied using the algorithm specified in settings['algorithm']
    and its corresponding parameters in settings['params']. To ensure globally unique cluster IDs
    (with 0 reserved for "no cluster"), an offset starting at 1 is used.

    Parameters:
        chain_df (pd.DataFrame): DataFrame for a chain (e.g., heavy or light).
        settings (dict): Dictionary of clustering settings, including:
                         - 'group_by': column name(s) to group by (optional).
                         - 'algorithm': clustering algorithm to use (e.g., "Iterative CDR3 similarity",
                                        "Matrix CDR3 similarity", "Hierarchical CDR3 Clustering").
                         - 'params': dictionary of algorithm-specific parameters.

    Returns:
        tuple: (cluster_series, representative_series)
            - cluster_series (pd.Series): Unique numeric cluster IDs for each row.
            - representative_series (pd.Series): Representative BCR ID (i.e., the index of the representative) for each row.
    """
    # Determine the clustering function and algorithm-specific parameters.
    algorithm = settings.get('algorithm', 'Iterative CDR3 similarity')
    params = settings.get('params', {})  # Specific parameters for the algorithm
    iterations = params.get('iterations', 1)

    if algorithm == 'Iterative CDR3 similarity':
        clustering_function = iterative_greedy_cdr3_clustering  # Modified to return (clusters, reps)
    elif algorithm == 'Matrix CDR3 similarity':
        clustering_function = matrix_greedy_cdr3_clustering  # Modified to return (clusters, reps)
    elif algorithm == 'Hierarchical CDR3 Clustering':
        clustering_function = hierarchical_cdr3_clustering  # Modified to return (clusters, reps)
    else:
        raise ValueError(f"Algorithm '{algorithm}' is not implemented.")

    best_result = None
    best_reps = None
    best_singles = np.inf
    np.random.seed(42)  # Ensure reproducibility for the first iteration
    all_num_singles = []

    # Iterative search for the best clustering (fewest singles)
    for i in range(iterations):
        random_state = 42 if iterations == 1 else np.random.randint(0, int(1e6))
        permuted_df = chain_df.sample(frac=1, random_state=random_state).copy()

        # Split the DataFrame into subgroups if grouping is specified.
        if settings.get('group_by'):
            group_cols = settings['group_by']
            group_df = permuted_df.copy()
            if isinstance(group_cols, list):
                for col in group_cols:
                    group_df[col] = group_df[col].fillna("missing")
            else:
                group_df[group_cols] = group_df[group_cols].fillna("missing")
            sub_dfs = [group for _, group in group_df.groupby(group_cols)]
        else:
            sub_dfs = [permuted_df]

        cluster_ids_list = []
        rep_series_list = []
        offset = 1  # Starts at 1 because 0 is reserved for "no cluster"

        # Iterate over subgroups
        for sub_df in sub_dfs:
            # Call the chosen clustering method.
            subgroup_clusters, subgroup_reps = clustering_function(sub_df, params)

            # Adjust cluster IDs by adding the current offset.
            subgroup_clusters = subgroup_clusters + offset
            # rep_series contains the original BCR IDs (no offset needed).

            if not subgroup_clusters.empty:
                offset = subgroup_clusters.max() + 1

            cluster_ids_list.append(subgroup_clusters)
            rep_series_list.append(subgroup_reps)

        concatenated_clusters = pd.concat(cluster_ids_list).sort_index()
        concatenated_reps = pd.concat(rep_series_list).sort_index()

        # Evaluate: Count "singles" (clusters with only 1 element)
        cluster_counts = concatenated_clusters.value_counts()
        num_singles = (cluster_counts == 1).sum()
        all_num_singles.append(num_singles)

        # Select the iteration with the fewest singles.
        if num_singles < best_singles:
            best_singles = num_singles
            best_result = concatenated_clusters
            best_reps = concatenated_reps

    # print("Number of singles per iteration:", all_num_singles)
    return best_result, best_reps


def normalized_levenshtein(seq1, seq2):
    """
    Computes the normalized Levenshtein distance between two sequences.

    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        float: Normalized Levenshtein distance (0 if both sequences are empty).
    """
    lev = distance.Levenshtein.distance(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    return 0 if max_len == 0 else lev / max_len


def standardize_missing_values(df, missing_indicators=["N/A", "nan", "NaN", "", "n/a"]):
    """
    Replaces common missing value indicators in the DataFrame with np.nan.

    Parameters:
        df (pd.DataFrame): The DataFrame to process.
        missing_indicators (list): List of values to be considered as missing.

    Returns:
        pd.DataFrame: The DataFrame with missing indicators replaced by np.nan.
    """
    return df.replace(missing_indicators, np.nan)


def assign_clone_names(df):
    """
    Assigns clone names based on cluster sizes (e.g., Clone-1, Clone-2, ...).

    Clusters with only one member are labeled "Non-clonal". If the combined cluster (HC-LC_CLUSTER)
    is "0-0", it is labeled as "Undefined" and the 'IS_CLONAL' flag is set to False.
    Finally, the DataFrame is sorted by cluster size (descending), with "0-0" entries placed at the very end.

    Parameters:
        df (pd.DataFrame): DataFrame containing BCR clustering information.

    Returns:
        pd.DataFrame: DataFrame with clone names assigned in the ('SAMPLE_INFORMATION', 'CLONE') column,
                      and cluster sizes in the ('SAMPLE_INFORMATION', 'CLUSTER_SIZE') column.
    """
    cluster_col = (sample_col, hclc_cluster_subcol)

    # Compute cluster sizes based on the combined cluster column.
    cluster_sizes = df.groupby(cluster_col).size()
    cluster_size_map = cluster_sizes.to_dict()

    # Create mapping for clone names.
    clone_name_map = {}

    # For clusters with more than one member, assign clone names (sorted by size descending).
    multi_member_clusters = {cid: size for cid, size in cluster_size_map.items() if size > 1}
    sorted_clusters = sorted(multi_member_clusters.items(), key=lambda x: x[1], reverse=True)
    for rank, (cid, size) in enumerate(sorted_clusters, start=1):
        clone_name_map[cid] = f"Clone-{rank}"

    # For clusters with only one member, assign "Non-clonal".
    for cid, size in cluster_size_map.items():
        if size == 1:
            clone_name_map[cid] = "Non-clonal"

    # Override: if the combined cluster is "0-0", label it as "Undefined".
    if "0-0" in clone_name_map:
        clone_name_map["0-0"] = "Undefined"

    # Map the clone names to a new column ('SAMPLE_INFORMATION', 'CLONE').
    df[(sample_col, clone_subcol)] = df[cluster_col].map(clone_name_map)
    df[(sample_col, clone_subcol)] = df[(sample_col, clone_subcol)].fillna("Non-clonal")

    # Map the cluster sizes into a new column for sorting and display.
    df[(sample_col, cluster_size_subcol)] = df[cluster_col].map(cluster_size_map)

    # Create a temporary column to flag undefined clusters ("0-0").
    df['_undefined'] = (df[cluster_col] == "0-0").astype(int)

    # Extract clone rank from the clone name for final sorting.
    df['clone_rank'] = df[(sample_col, clone_subcol)].apply(parse_clone_number)

    # Sort the DataFrame: first by _undefined ascending, then by cluster size descending, then by clone rank.
    df = df.sort_values(
        by=['_undefined', (sample_col, cluster_size_subcol), 'clone_rank'],
        ascending=[True, False, True]
    ).drop(columns=['clone_rank'])

    # Drop the temporary flag column.
    df = df.drop(columns=['_undefined'])

    # Ensure that rows with combined cluster "0-0" have IS_CLONAL set to False.
    df.loc[df[cluster_col] == "0-0", (sample_col, clonal_subcol)] = False

    return df


def color_clones(df, config):
    """
    Assigns colors to clones based on a discrete colormap and stores the corresponding hex codes in a column.

    Clones labeled "Non-clonal" or "Undefined" receive default colors as defined in the configuration.
    The resulting colors are stored in the ('SAMPLE_INFORMATION', 'CLONE_COLOR') column.

    Parameters:
        df (pd.DataFrame): DataFrame containing clone labels in the ('SAMPLE_INFORMATION', 'CLONE') column.
        config (dict): Dictionary containing color configuration, e.g.:
                       {
                           'palette': 'GnBu_d',
                           'non_clonal': '#D3D3D3',  # light gray
                           'undefined': '#FFFFFF'    # white
                       }

    Returns:
        pd.DataFrame: DataFrame with an additional column ('SAMPLE_INFORMATION', 'CLONE_COLOR')
                      containing the assigned hex color for each clone.
    """
    # Retrieve configuration parameters.
    palette_name = config.get('palette', 'GnBu_d')
    default_non_clonal = config.get('non_clonal', "#D3D3D3")
    default_undefined = config.get('undefined', "#FFFFFF")

    # Get unique clone labels from the CLONE column.
    unique_clones = df[(sample_col, clone_subcol)].unique()
    unique_clones = [clone for clone in unique_clones if pd.notna(clone)]

    # Determine which clones receive a random color.
    # "Non-clonal" and "Undefined" will be assigned default colors.
    clones_to_color = [clone for clone in unique_clones if clone not in ["Non-clonal", "Undefined"]]
    n_colors = len(clones_to_color)

    # Retrieve the discrete colormap from Matplotlib.
    cmap = mpl.cm.get_cmap(palette_name, n_colors)
    # Convert colormap colors to hex strings.
    colors = [mpl.colors.rgb2hex(cmap(i)) for i in range(n_colors)]
    random.shuffle(colors)

    # Build a dictionary mapping clone labels to colors.
    clone_color_map = {}
    for clone in unique_clones:
        if clone == "Non-clonal":
            clone_color_map[clone] = default_non_clonal
        elif clone == "Undefined":
            clone_color_map[clone] = default_undefined
        elif clone in clones_to_color:
            clone_color_map[clone] = colors.pop()
        else:
            clone_color_map[clone] = default_non_clonal  # Fallback

    # Map the clone colors to the ('SAMPLE_INFORMATION', 'CLONE_COLOR') column.
    df[(sample_col, clone_color_subcol)] = df[(sample_col, clone_subcol)].map(clone_color_map)

    return df


def save_excel_with_row_colors(df, path_to_file, color_col=(sample_col, clone_color_subcol), sheet_name='Clustered'):
    """
    Saves the DataFrame as an Excel file and applies row coloring based on hex color codes.

    The DataFrame is exported to an Excel file with its index. Cell background colors (excluding the index column)
    are applied based on the hex color values from the specified color column. If an index name exists, the extra
    header row is taken into account.

    Parameters:
        df (pd.DataFrame): The DataFrame to save.
        path_to_file (str or Path): Output file path.
        color_col (str or tuple): The column (or MultiIndex tuple) in df containing hex color codes.
        sheet_name (str): Worksheet name in the Excel file.

    Returns:
        None
    """

    # Convert <NA> (pd.NA) to np.nan to ensure they export as real empty Excel cells
    df = df.astype(object).where(pd.notnull(df), np.nan)

    with pd.ExcelWriter(path_to_file, engine='xlsxwriter') as writer:
        # Write the DataFrame with its index.
        df.to_excel(writer, index=True, sheet_name=sheet_name)
        workbook = writer.book
        worksheet = writer.sheets[sheet_name]

        # Determine the number of header rows (e.g., 2 for a MultiIndex on columns).
        header_rows = df.columns.nlevels if isinstance(df.columns, pd.MultiIndex) else 1

        # Account for an extra row if an index name exists.
        if df.index.name is not None:
            header_rows += 1

        # Total number of columns = DataFrame columns + 1 (for the index).
        total_cols = df.shape[1] + 1

        # Iterate over each row of the DataFrame.
        for i in range(df.shape[0]):
            # Excel row index = header_rows + i.
            excel_row = header_rows + i
            # Retrieve the hex color for the row.
            cell_color = df.iloc[i][color_col]
            if pd.isna(cell_color) or not isinstance(cell_color, str):
                continue  # Skip rows without a valid color.

            cell_format = workbook.add_format({'bg_color': cell_color})

            # Column 0 is the index; start from column 1.
            for col_idx in range(1, total_cols):
                value = df.iloc[i, col_idx - 1]
                if not np.isscalar(value):
                    # If value is a list or array, convert it to a string.
                    worksheet.write(excel_row, col_idx, str(value), cell_format)
                elif pd.isna(value) or (isinstance(value, (int, float)) and not np.isfinite(value)):
                    worksheet.write(excel_row, col_idx, '', cell_format)
                else:
                    worksheet.write(excel_row, col_idx, value, cell_format)


def plot_clonality_donut_chart(df, clone_col, color_col):
    """
    Plots a clonality donut chart based on clone counts.

    Clones are sorted by size (largest first) and, among clones with equal counts,
    by their first occurrence in the DataFrame. Non-clonal entries ("Non-clonal" and "Undefined")
    are appended at the end regardless of size. The function uses clone counts and assigns colors
    based on the provided hex color codes.

    Parameters:
        df (pd.DataFrame): DataFrame containing clone information.
        clone_col (str): Column name with clone labels (e.g., "CLONE").
        color_col (str): Column name with corresponding hex color codes.

    Returns:
        matplotlib.figure.Figure: A Matplotlib Figure object representing the donut chart.
    """
    total = len(df)

    # Calculate clone counts (unsorted).
    clone_counts = df[clone_col].value_counts()

    # Record the first occurrence of each clone based on DataFrame order.
    first_occ = {}
    for pos, clone in enumerate(df[clone_col]):
        if clone not in first_occ:
            first_occ[clone] = pos

    # Separate clonal from non-clonal labels.
    clonal_labels = [label for label in clone_counts.index if label not in ["Non-clonal", "Undefined"]]
    non_clonal_labels = [label for label in clone_counts.index if label in ["Non-clonal", "Undefined"]]

    # Sort clonal labels: by count descending, then by first occurrence ascending.
    clonal_order = sorted(clonal_labels, key=lambda x: (-clone_counts[x], first_occ.get(x, float('inf'))))
    # Sort non-clonal labels by first occurrence.
    non_clonal_order = sorted(non_clonal_labels, key=lambda x: first_occ.get(x, float('inf')))
    final_labels = clonal_order + non_clonal_order

    # Build the color mapping: extract colors in the order of final_labels.
    colors_dict = df[[clone_col, color_col]].drop_duplicates().set_index(clone_col)[color_col].to_dict()
    colors = [colors_dict.get(label, "#000000") for label in final_labels]

    sizes = [clone_counts[label] for label in final_labels]

    # Create legend labels with clone counts and percentages.
    legend_labels = []
    for label in final_labels:
        count = clone_counts[label]
        percentage = (count / total * 100) if total > 0 else 0
        legend_labels.append(f"{label}: {count} ({round(percentage, 1)}%)")

    legend_dict = {
        'labels': legend_labels,
        'loc': 'center left',
        'title': 'Legend'
    }

    fig = plot_pie_chart(
        pie_sizes=sizes,
        pie_colors=colors,
        center_text=str(total),
        center_fontsize=14,
        counterclockwise=False,
        edges={'color': 'black', 'ewidth': 0.5},
        bbox_to_anchor_tuple=(0.8, 0.5),
        legend=legend_dict,
        figsize=(4, 2),
        legend_font_size=8,
        max_legend_items_per_col=7,
        fixed_pie_position=[0.0, 0.1, 0.7, 0.8],
        wedge_width=0.5,
        adjust_kwargs={"left": 0.0, "bottom": 0.3}
    )

    return fig


# ======================
# Clustering algorithms
# ======================

# Original cluster approach
def iterative_greedy_cdr3_clustering(sub_df, params):
    """
    Iteratively performs greedy clustering of CDR3 sequences based on normalized Levenshtein distance
    and an optional length difference filter.

    This function returns two pandas Series:
      - cluster_series (pd.Series): Cluster IDs (int) for each row.
      - rep_series (pd.Series): Representative BCR IDs (i.e., the index of the first candidate used
        to initialize the cluster) for each row.

    Rows with missing CDR3_AA values are assigned cluster 0 and have no representative.

    Parameters:
        sub_df (pd.DataFrame): DataFrame containing CDR3 sequences (in the column specified by `cdr3_aa_subcol`)
                               and other related information.
        params (dict): Dictionary of clustering parameters. Expected keys include:
                       - 'lev_threshold' (float): Normalized Levenshtein distance threshold (default: 0.25).
                       - 'length_threshold' (int or None): Maximum allowed length difference (default: 1).

    Returns:
        tuple: A tuple (cluster_series, rep_series) where:
               - cluster_series (pd.Series): Cluster IDs for each row.
               - rep_series (pd.Series): Representative BCR IDs for each row.
    """
    # Separate rows with valid and missing CDR3_AA sequences.
    valid_df = sub_df[sub_df[cdr3_aa_subcol].notna()].copy()
    missing_df = sub_df[sub_df[cdr3_aa_subcol].isna()].copy()

    # If no valid sequences exist, assign cluster 0 and no representative to all rows.
    if valid_df.empty:
        cluster_series = pd.Series(0, index=sub_df.index)
        rep_series = pd.Series([None] * len(sub_df), index=sub_df.index)
        return cluster_series, rep_series

    # List of sequences and corresponding indices for valid rows.
    sequences = valid_df[cdr3_aa_subcol].tolist()
    valid_indices = list(valid_df.index)
    n = len(sequences)

    # Retrieve filtering parameters.
    lev_threshold = params.get('lev_threshold', 0.25)
    length_threshold = params.get('length_threshold', 1)

    # Function to decide if two sequences pass the length difference filter.
    def passes_length_filter(seq1, seq2):
        diff = abs(len(seq1) - len(seq2))
        if length_threshold in [False, None]:
            return True
        elif length_threshold == 0:
            return diff == 0
        else:
            return diff <= length_threshold

    # Initialize cluster labels for valid sequences.
    cluster_labels = np.zeros(n, dtype=int)
    current_cluster = 1  # Start numbering clusters at 1 (0 is reserved for missing/unclusterable).
    unclustered = list(range(n))  # List of indices for unclustered sequences.
    rep_dict = {}  # Dictionary to store the representative for each cluster.

    # Iterative greedy clustering loop.
    while unclustered:
        # Start a new cluster with the first candidate in the unclustered list.
        rep_idx = unclustered[0]
        rep_dict[current_cluster] = valid_indices[rep_idx]  # The first selected sequence is the representative.
        current_cluster_members = [rep_idx]

        # Remove the representative from the unclustered list.
        unclustered.remove(rep_idx)

        # Iterate over a copy of the unclustered list to safely remove elements.
        for candidate in unclustered.copy():
            candidate_seq = sequences[candidate]
            valid_candidate = True
            # Check candidate against all sequences already in the cluster.
            for member in current_cluster_members:
                member_seq = sequences[member]
                # First, check the length difference filter.
                if not passes_length_filter(candidate_seq, member_seq):
                    valid_candidate = False
                    break
                # Then, compute normalized Levenshtein distance.
                if normalized_levenshtein(candidate_seq, member_seq) > lev_threshold:
                    valid_candidate = False
                    break
            if valid_candidate:
                current_cluster_members.append(candidate)
                unclustered.remove(candidate)

        # Assign the current cluster number to all members in the cluster.
        for idx in current_cluster_members:
            cluster_labels[idx] = current_cluster
        current_cluster += 1

    rep_series_valid = pd.Series({idx: rep_dict[cluster] for idx, cluster in zip(valid_indices, cluster_labels)})
    rep_series_missing = pd.Series([None] * len(missing_df), index=missing_df.index)

    cluster_series = pd.concat(
        [pd.Series(cluster_labels, index=valid_indices), pd.Series(0, index=missing_df.index)]
    ).sort_index()
    rep_series = pd.concat([rep_series_valid, rep_series_missing]).sort_index()

    return cluster_series, rep_series


def matrix_greedy_cdr3_clustering(sub_df, params):
    """
    Performs deterministic greedy clustering of CDR3 sequences using a full pairwise similarity matrix.

    This function returns two pandas Series:
      - cluster_series (pd.Series): Cluster IDs (int) for each row.
      - rep_series (pd.Series): Representative BCR IDs for each row, based on the candidate with maximal overlap.

    Rows with missing CDR3_AA values are assigned cluster 0 and have no representative.

    Parameters:
        sub_df (pd.DataFrame): DataFrame containing CDR3 sequences (in the column specified by `cdr3_aa_subcol`)
                               and other related information.
        params (dict): Dictionary of clustering parameters. Expected keys include:
                       - 'lev_threshold' (float): Normalized Levenshtein distance threshold (default: 0.25).
                       - 'length_threshold' (int or None): Maximum allowed length difference.

    Returns:
        tuple: A tuple (cluster_series, rep_series) where:
               - cluster_series (pd.Series): Cluster IDs for each row.
               - rep_series (pd.Series): Representative BCR IDs for each row.
    """
    valid_df = sub_df[sub_df[cdr3_aa_subcol].notna()].copy()
    missing_df = sub_df[sub_df[cdr3_aa_subcol].isna()].copy()

    if valid_df.empty:
        cluster_series = pd.Series(0, index=sub_df.index)
        rep_series = pd.Series([None] * len(sub_df), index=sub_df.index)
        return cluster_series, rep_series

    sequences = valid_df[cdr3_aa_subcol].tolist()
    b_cell_ids = list(valid_df.index)
    n = len(sequences)

    length_threshold = params.get('length_threshold', None)
    lengths = np.array([len(seq) for seq in sequences])
    if length_threshold in [False, None]:
        length_filter = np.ones((n, n), dtype=bool)
    elif length_threshold == 0:
        length_filter = (np.abs(np.subtract.outer(lengths, lengths)) == 0)
    else:
        length_filter = (np.abs(np.subtract.outer(lengths, lengths)) <= length_threshold)

    # Build the full pairwise normalized Levenshtein distance matrix.
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            d = normalized_levenshtein(sequences[i], sequences[j])
            D[i, j] = d
            D[j, i] = d

    lev_threshold = params.get('lev_threshold', 0.25)
    lev_filter = (D <= lev_threshold)
    valid_filter = np.logical_and(length_filter, lev_filter)

    cluster_labels = np.zeros(n, dtype=int)
    current_cluster = 0  # Start at 0; an offset can be added later if needed.
    rep_dict = {}
    unclustered = set(range(n))

    while unclustered:
        # For each candidate, count how many sequences it can cluster with.
        counts = {i: np.sum([valid_filter[i, j] for j in unclustered]) for i in unclustered}
        rep = max(counts, key=counts.get)
        rep_dict[current_cluster] = b_cell_ids[rep]
        cluster_members = {j for j in unclustered if valid_filter[rep, j]}
        for j in cluster_members:
            cluster_labels[j] = current_cluster
        current_cluster += 1
        unclustered -= cluster_members

    rep_series_valid = pd.Series({idx: rep_dict[cluster] for idx, cluster in zip(b_cell_ids, cluster_labels)})
    rep_series_missing = pd.Series([None] * len(missing_df), index=missing_df.index)

    cluster_series = pd.concat(
        [pd.Series(cluster_labels, index=b_cell_ids), pd.Series(0, index=missing_df.index)]
    ).sort_index()
    rep_series = pd.concat([rep_series_valid, rep_series_missing]).sort_index()

    return cluster_series, rep_series


def hierarchical_cdr3_clustering(sub_df, params):
    """
    Performs hierarchical clustering on valid CDR3 sequences using normalized Levenshtein distance
    and an optional length difference filter.

    This function returns two pandas Series:
      - cluster_series (pd.Series): Cluster IDs (int) for each row, determined using fcluster.
      - rep_series (pd.Series): Representative BCR IDs for each row, determined as the medoid of the cluster
        (i.e., the BCR with the lowest mean distance to all other BCRs in the cluster).

    Rows with missing CDR3_AA values are assigned cluster 0 and have no representative.

    Parameters:
        sub_df (pd.DataFrame): DataFrame containing CDR3 sequences (in the column specified by `cdr3_aa_subcol`)
                               and other related information.
        params (dict): Dictionary of clustering parameters. Expected keys include:
                       - 'lev_threshold' (float): Normalized Levenshtein distance threshold (default: 0.25).
                       - 'length_threshold' (int or None): Maximum allowed length difference (default: 1).
                       - 'hc_method' (str): Linkage method to use in hierarchical clustering (default: 'average').

    Returns:
        tuple: A tuple (cluster_series, rep_series) where:
               - cluster_series (pd.Series): Cluster IDs for each row.
               - rep_series (pd.Series): Representative BCR IDs for each row.
    """
    valid_df = sub_df[sub_df[cdr3_aa_subcol].notna()].copy()
    missing_df = sub_df[sub_df[cdr3_aa_subcol].isna()].copy()

    if valid_df.empty:
        cluster_series = pd.Series(0, index=sub_df.index)
        rep_series = pd.Series([None] * len(sub_df), index=sub_df.index)
        return cluster_series, rep_series

    sequences = valid_df[cdr3_aa_subcol].tolist()
    b_cell_ids = list(valid_df.index)
    n = len(sequences)

    if n == 1:
        result_valid = pd.Series([1], index=b_cell_ids)
        rep_series_valid = pd.Series([b_cell_ids[0]], index=b_cell_ids)
        result_missing = pd.Series(0, index=missing_df.index)
        rep_series_missing = pd.Series([None] * len(missing_df), index=missing_df.index)
        cluster_series = pd.concat([result_valid, result_missing]).sort_index()
        rep_series = pd.concat([rep_series_valid, rep_series_missing]).sort_index()
        return cluster_series, rep_series

    lev_threshold = params.get('lev_threshold', 0.25)
    length_threshold = params.get('length_threshold', 1)
    hc_method = params.get('hc_method', 'average')

    # Build the pairwise distance matrix D.
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if length_threshold not in [False, None]:
                diff = abs(len(sequences[i]) - len(sequences[j]))
                if (length_threshold == 0 and diff != 0) or (length_threshold > 0 and diff > length_threshold):
                    D[i, j] = 1.1
                    D[j, i] = 1.1
                    continue
            d = normalized_levenshtein(sequences[i], sequences[j])
            D[i, j] = d
            D[j, i] = d

    condensed_D = squareform(D)
    Z = linkage(condensed_D, method=hc_method)
    cluster_labels = fcluster(Z, t=lev_threshold, criterion='distance')

    # Determine the medoid for each cluster.
    rep_dict = {}
    unique_clusters = np.unique(cluster_labels)
    for cluster in unique_clusters:
        indices = [i for i, lab in enumerate(cluster_labels) if lab == cluster]
        if len(indices) == 1:
            rep_dict[cluster] = b_cell_ids[indices[0]]
        else:
            sub_D = D[np.ix_(indices, indices)]
            sum_dist = sub_D.sum(axis=1)  # Sum of distances per candidate.
            medoid_local_idx = indices[np.argmin(sum_dist)]
            rep_dict[cluster] = b_cell_ids[medoid_local_idx]

    rep_series_valid = pd.Series({idx: rep_dict[cluster] for idx, cluster in zip(b_cell_ids, cluster_labels)})
    rep_series_missing = pd.Series([None] * len(missing_df), index=missing_df.index)

    cluster_series = pd.concat(
        [pd.Series(cluster_labels, index=b_cell_ids), pd.Series(0, index=missing_df.index)]
    ).sort_index()
    rep_series = pd.concat([rep_series_valid, rep_series_missing]).sort_index()

    return cluster_series, rep_series


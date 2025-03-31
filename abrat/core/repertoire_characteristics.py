import pandas as pd
import numpy as np
from difflib import SequenceMatcher

import plotly.graph_objects as go
import peptides
import streamlit as st

def get_sample_information(df):
    """Extracts Sample information from dataframe and returns a dictionary"""
    infos = {'Cohort': ', '.join(df[('SAMPLE_INFORMATION', 'COHORT')].unique()),
             'Donor': ', '.join(df[('SAMPLE_INFORMATION', 'SUBJECT')].unique()),
             'Time points': ', '.join(df[('SAMPLE_INFORMATION', 'TIME_POINT')].unique()),
             'Sample': ', '.join(df[('SAMPLE_INFORMATION', 'TISSUE')].unique()),
             'Subset': ', '.join(df[('SAMPLE_INFORMATION', 'SUBSET')].unique()), 'Total BCRs': df.shape[0],
             'Heavy chains': int(df[('HEAVY_CHAIN', 'HEAVY_FOUND')].sum()),
             'Light chains': int(df[('LIGHT_CHAIN', 'LIGHT_FOUND')].sum()),
             'Kappa chains': int((df[('LIGHT_CHAIN', 'PCR_ISOTYPE')].dropna()=="KC").sum()),
             'Lambda chains': int((df[('LIGHT_CHAIN', 'PCR_ISOTYPE')].dropna()=="LC").sum())
             }
    # infos['Paired chains'] = df[('HEAVY_CHAIN','HEAVY_FOUND')].sum()
    return infos

def get_hydrophobicity_values(sequence, scale='kyte-doolittle'):
    """
    Returns a list of hydrophobicity values based on the given amino acid sequence
    and the selected scale.

    :param sequence: Amino acid sequence (single-letter code).
    :param scale: 'kyte-doolittle' or 'eisenberg' (default: 'kyte-doolittle').

    :return: List[float]: A list of the corresponding hydrophobicity values.

    :raises: ValueError: If an unknown amino acid symbol is encountered or an unsupported scale is specified.
    """

    # Kyte-Doolittle scale: Hydrophobicity values for each amino acid residue
    kyte_doolittle = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
        'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    # Eisenberg scale: Consensus scale (values according to Eisenberg et al.)
    eisenberg = {
        'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
        'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
        'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
        'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
    }

    # Select the scale based on the parameter
    scale = scale.lower()
    if scale == 'kyte-doolittle':
        scale_dict = kyte_doolittle
    elif scale == 'eisenberg':
        scale_dict = eisenberg
    else:
        raise ValueError("Unsupported scale. Please choose 'kyte-doolittle' or 'eisenberg'.")

    # Generate a list of hydrophobicity values
    values = []
    for aa in sequence.upper():
        if aa in scale_dict:
            values.append(scale_dict[aa])
        else:
            raise ValueError(f"Unknown amino acid '{aa}' in sequence.")

    return values


def compute_gravy_scores(sequences: pd.Series, scale: str = 'kyte-doolittle') -> pd.Series:
    """
    Computes the GRAVY (Grand Average of Hydropathy) scores for a Pandas Series of
    CDR3 amino acid sequences using the specified hydrophobicity scale.

    Parameters:
      swquences (pd.Series): Series of CDR3 amino acid sequences (using single-letter codes).
      scale (str): Hydrophobicity scale to use ('kyte-doolittle' or 'eisenberg').

    :returns: pd.Series: Series of GRAVY scores calculated for each sequence.
    """

    def compute_gravy(seq: str) -> float:
        # Get the hydrophobicity values for the given sequence based on the chosen scale.
        values = get_hydrophobicity_values(seq, scale=scale)
        # Check if sequence is not empty to avoid division by zero.
        if len(values) == 0:
            return None
        # Calculate the GRAVY score as the average of the hydrophobicity values.
        return sum(values) / len(values)

    # Apply the computation to each sequence in the series.
    return sequences.apply(compute_gravy)


def compute_histogram(series: pd.Series, bin_size: float = 1, min_val: float=None, max_val: float=None) -> pd.Series:
    """
    Computes a histogram from a Pandas Series using the specified bin size.

    Parameters:
      series (pd.Series): The input data as a numeric Pandas Series.
      bin_size (float): The size of the bins. Default is 1.
      min_val (float): min bin value
      max_val (float): max bin value

    Returns:
      pd.Series: A Pandas Series where the index represents the start of each bin
                 and the values are the counts of elements in each bin.
    """
    # Ensure the series contains numeric data
    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError("The input series must be numeric.")

    # Determine the minimum and maximum values in the series
    if min_val is None:
        min_val = series.min()
    if max_val is None:
        max_val = series.max()

    # Create bins from the minimum to the maximum value, ensuring the maximum is included
    bins = np.arange(min_val, max_val + bin_size, bin_size)

    # Compute the histogram counts and bin edges using numpy.histogram
    counts, bin_edges = np.histogram(series, bins=bins)

    # Calculate the midpoint for each bin as the average of the bin edges
    index = bin_edges[:-1]
    # [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]

    # Return the histogram as a Pandas Series with the index as index
    hist_series = pd.Series(counts, index=index)
    return hist_series

def compute_diversity_index(sequences: pd.Series, index_type: str = "shannon") -> float:
    """
    Computes a diversity index for a Series of CDR3 amino acid sequences.

    The function calculates the frequency distribution of unique sequences and then:
      - For the Shannon index: H = -∑ p_i * ln(p_i)
      - For the Inverse Simpson index: D = 1 / ∑ p_i^2

    Parameters:
        sequences (pd.Series): A Pandas Series containing CDR3 sequences (strings).
        index_type (str): The type of diversity index to compute. Use "shannon"
                          for the Shannon index or "inverse_simpson" for the Inverse Simpson index.

    Returns:
        float: The computed diversity index.

    Raises:
        ValueError: If an unknown index_type is provided.
    """
    # Get counts for each unique sequence
    counts = sequences.value_counts()
    total = counts.sum()

    # Calculate the proportion (p_i) for each unique sequence
    probs = counts / total

    if index_type.lower() == "shannon":
        # Shannon index: H = -sum(p_i * ln(p_i))
        shannon_index = -np.sum(probs * np.log(probs))
        return shannon_index
    elif index_type.lower() == "inverse_simpson":
        # Inverse Simpson index: D = 1 / sum(p_i^2)
        inverse_simpson = 1.0 / np.sum(probs ** 2)
        return inverse_simpson
    else:
        raise ValueError("index_type must be either 'shannon' or 'inverse_simpson'.")


def subsample_series_dict(series_dict: dict) -> tuple:
    """
    Determines the smallest Series (by length) from a dictionary of Series and returns:
    - the original dictionary (unchanged)
    - the minimal size,
    - and the key of the smallest Series.

    Parameters:
        series_dict (dict): Dictionary where keys are dataset names and values are Pandas Series.
        random_state (int): Seed for reproducibility.

    Returns:
        tuple: (series_dict, min_size, smallest_key)
    """
    sizes = {key: series.shape[0] for key, series in series_dict.items()}
    min_size = min(sizes.values())
    smallest_key = min(sizes, key=sizes.get)
    return series_dict, min_size, smallest_key


def compute_diversity_indices_for_subsampling(
        series_dict: dict,
        index_type: str = "shannon",
        iterations: int = 20,
        random_state: int = 42,
        seq_column: str = "CDR3_AA"
) -> dict:
    """
    For each Pandas Series in series_dict (which should contain CDR3 sequences),
    performs multiple iterations of subsampling (for those larger than the smallest)
    and computes the diversity index (e.g., Shannon) for each iteration.

    For the smallest Series, the index is computed directly.

    Parameters:
        series_dict (dict): Dictionary where keys are dataset names and values are Pandas Series.
        index_type (str): "shannon" or "inverse_simpson".
        iterations (int): Number of subsampling iterations.
        random_state (int): Seed for reproducibility.
        seq_column (str): The name to assign to the Series (for clarity; not used internally).

    Returns:
        dict: A dictionary with keys corresponding to dataset names and values being a dict of (n, mean, std)
              for the computed diversity index.
    """
    # Bestimme die Länge jeder Series und den kleinsten Datensatz
    series_dict, min_size, smallest_key = subsample_series_dict(series_dict)

    results = {}
    for key, series in series_dict.items():
        # Optional: falls die Series einen anderen Namen hat, benenne sie um:
        series = series.rename(seq_column)

        if series.shape[0] == min_size:
            # Für den kleinsten Datensatz direkt berechnen.
            diversity_value = compute_diversity_index(series, index_type=index_type)
            results[key] = {'n':min_size, 'mean':diversity_value, 'std':0.0}
        else:
            # Für größere Series mehrfach subsamplen.
            indices = []
            for i in range(iterations):
                subsample = series.sample(n=min_size, random_state=random_state + i)
                index_value = compute_diversity_index(subsample, index_type=index_type)
                indices.append(index_value)
            mean_val = np.mean(indices)
            std_val = np.std(indices)
            results[key] = {'n':min_size, 'mean':mean_val, 'std':std_val}
    return results


# CURRENTLY NOT IMPLEMENTED
def compute_pairwise_jaccard(data_dict: dict, iterations: int = 20, random_state: int = 42):
    """
    Computes the pairwise Jaccard index for multiple sequence datasets (e.g. CDR3 repertoires).

    For each dataset (key in data_dict), the function:
      - Extracts the unique sequences.
      - Determines the minimal set size across all datasets.
      - Performs 'iterations' rounds of subsampling (without replacement) from each dataset,
        sampling 'min_size' sequences.
      - Computes the Jaccard index for each pair of datasets:
          J(A,B) = |A ∩ B| / |A ∪ B|

    Returns:
      - results: A dictionary with keys as tuples (set1, set2) and values as dictionaries:
                 {'n': min_size, 'mean': mean_jaccard, 'std': std_jaccard}
      - matrix: A symmetric Pandas DataFrame with dataset names as rows and columns containing the mean Jaccard values.

    Parameters:
        data_dict (dict): Dictionary where keys are dataset names and values are lists or pandas Series of sequences (str, e.g. CDR3s).
        iterations (int): Number of subsampling iterations (default 20).
        random_state (int): Seed for reproducibility.
    """
    # Get unique sequences only
    unique_dict = {}
    for dataset, seqs in data_dict.items():
        if isinstance(seqs, pd.Series):
            seq_list = seqs.dropna().unique().tolist()
        else:
            seq_list = list(set(seqs))
        unique_dict[dataset] = seq_list

    # Determine the minimal set size
    min_size = min(len(seq_list) for seq_list in unique_dict.values())

    # Sort dataset names to ensure consistent key ordering
    dataset_names = sorted(list(unique_dict.keys()))

    # Create dictionary for pairwise results for all combinations.
    # Keys are tuples (dataset1, dataset2) with dataset1 <= dataset2 (alphabetically).
    pairwise_results = {}
    for i in range(len(dataset_names)):
        for j in range(i, len(dataset_names)):
            pairwise_results[(dataset_names[i], dataset_names[j])] = []

    rng = np.random.RandomState(random_state)

    # Subsampling and calculation of Jaccard indices
    for _ in range(iterations):
        subsamples = {}
        for dataset, seq_list in unique_dict.items():
            # Randomly choose without replacement
            subsample = rng.choice(seq_list, size=min_size, replace=False)
            subsamples[dataset] = set(subsample)
        # Calculate Jaccard index for each combination
        for i in range(len(dataset_names)):
            for j in range(i, len(dataset_names)):
                setA = subsamples[dataset_names[i]]
                setB = subsamples[dataset_names[j]]
                union = setA | setB
                intersection = setA & setB
                jacc = len(intersection) / len(union) if len(union) > 0 else 0
                pairwise_results[(dataset_names[i], dataset_names[j])].append(jacc)

    # Determine means and standard deviation for each pair
    results = {}
    for pair, values in pairwise_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        results[pair] = {'n': min_size, 'mean': mean_val, 'std': std_val}

    # Build a symmetric matrix of the results.
    matrix = pd.DataFrame(index=dataset_names, columns=dataset_names, dtype=float)
    for i in range(len(dataset_names)):
        for j in range(len(dataset_names)):
            # Ensure key is ordered alphabetically
            key = tuple(sorted((dataset_names[i], dataset_names[j])))
            matrix.iloc[i, j] = results[key]['mean']

    return results, matrix

def sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculates the similarity between two sequences using difflib's SequenceMatcher.
    Returns a value between 0 (no similarity) and 1 (identical).
    """
    return SequenceMatcher(None, seq1, seq2).ratio()

# Not included in the basic repertoire characteristics
def compute_pairwise_weighted_jaccard(data_dict: dict, iterations: int = 20, random_state: int = 42):
    """
    Computes the weighted pairwise Jaccard Index between multiple sequence datasets,
    taking into account not only exact matches but also similar CDR3s.

    For each dataset (key in data_dict):
      - Unique sequences are extracted.
      - The minimum dataset size is determined.
      - 'iterations' rounds of subsampling (without replacement) are performed,
        sampling 'min_size' sequences.
      - For each pair, the weighted Jaccard Index is computed:
          weighted J(A,B) = soft_intersection / soft_union,
      where:
          soft_intersection = sum over a in A of max_{b in B}(similarity(a, b))
          soft_union = |A| + |B| - soft_intersection

    Returns:
      - results: A dictionary with keys as tuples (set1, set2) and values as dictionaries:
                 {'n': min_size, 'mean': mean_weighted_jaccard, 'std': std_weighted_jaccard}
      - matrix: A symmetric Pandas DataFrame with dataset names as rows and columns,
                containing the mean weighted Jaccard values.
    """
    # Extract unique sequences
    unique_dict = {}
    for dataset, seqs in data_dict.items():
        if isinstance(seqs, pd.Series):
            seq_list = seqs.dropna().unique().tolist()
        else:
            seq_list = list(set(seqs))
        unique_dict[dataset] = seq_list

    # Determine the minimal size
    min_size = min(len(seq_list) for seq_list in unique_dict.values())

    # Sort dataset names for consistent key ordering
    dataset_names = sorted(list(unique_dict.keys()))

    # Prepare dictionary for pairwise results for all combinations.
    # Keys are tuples (dataset1, dataset2) with dataset1 <= dataset2 (alphabetically).
    pairwise_results = {}
    for i in range(len(dataset_names)):
        for j in range(i, len(dataset_names)):
            pairwise_results[(dataset_names[i], dataset_names[j])] = []

    rng = np.random.RandomState(random_state)

    # Subsampling and calculation of weighted Jaccard indices
    for _ in range(iterations):
        subsamples = {}
        for dataset, seq_list in unique_dict.items():
            # Randomly choose without replacement
            subsample = rng.choice(seq_list, size=min_size, replace=False)
            subsamples[dataset] = set(subsample)
        # Calculate the weighted Jaccard Index for each pair
        for i in range(len(dataset_names)):
            for j in range(i, len(dataset_names)):
                setA = subsamples[dataset_names[i]]
                setB = subsamples[dataset_names[j]]
                # Soft Intersection: For each a in setA, take the maximum similarity score to any b in setB
                soft_intersection = 0.0
                for a in setA:
                    # It is assumed that sequence_similarity() returns a value between 0 and 1
                    best_sim = max(sequence_similarity(a, b) for b in setB)
                    soft_intersection += best_sim
                # Define soft union analogous to the classical union
                union_score = len(setA) + len(setB) - soft_intersection
                jacc = soft_intersection / union_score if union_score > 0 else 0
                pairwise_results[(dataset_names[i], dataset_names[j])].append(jacc)

    # Compute means and standard deviations
    results = {}
    for pair, values in pairwise_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        results[pair] = {'n': min_size, 'mean': mean_val, 'std': std_val}

    # Represent results in a symmetric matrix
    matrix = pd.DataFrame(index=dataset_names, columns=dataset_names, dtype=float)
    for i in range(len(dataset_names)):
        for j in range(len(dataset_names)):
            # Retrieve the result using a key with consistent ordering
            key = tuple(sorted((dataset_names[i], dataset_names[j])))
            matrix.iloc[i, j] = results[key]['mean']

    return results, matrix

#=====================
# PLOTTING FUNCTIONS
#=====================

def create_heatmap_with_annotations(
        df: pd.DataFrame,
        colorscale: str = 'reds',
        margin_color = 'black',
        margin_width = 1,
        color_bar_title ="Colorbar"
) -> go.Figure:
    """
    Creates a heatmap with custom annotations for columns and rows.

    Parameters:
      - df: DataFrame containing the values for the heatmap.
            The columns and rows are used for the labels.
      - color_dict: Dictionary that specifies the text color for the column and row names.
      - colorscale: Colorscale for the heatmap (default: 'Viridis').
      - col_annotation_y: y offset (in paper coordinates) for column labels.
      - row_annotation_x: x offset (in paper coordinates) for row labels.

    Returns:
      - Plotly Figure (go.Figure) with the created heatmap and annotations.
    """
    n_rows, n_cols = df.shape

    # Zellenmittelpunkte: x=0.5..(n_cols-0.5), y=0.5..(n_rows-0.5)
    xvals = np.arange(n_cols) + 0.5
    yvals = np.arange(n_rows) + 0.5

    # Heatmap anlegen (numeric axes), horizontale Colorbar
    fig = go.Figure(data=go.Heatmap(
        x=xvals,
        y=yvals,
        z=df.values,
        colorscale=colorscale,
        colorbar=dict(
            orientation='h',  # horizontal
            x=0,  # linksbündig
            y=-0.1,  # näher an der Heatmap
            xanchor='left',
            yanchor='top',
            thickness=10,  # dünnere Leiste
            len=1,  # so breit wie die Heatmap-Domain
            outlinecolor = margin_color,  # Farbe der Colorbar-Umrandung
            outlinewidth = margin_width, # Dicke der Umrandung
            title=dict(  # Titel für die Colorbar
                text=color_bar_title,  # Text des Titels
                side='bottom',  # Titel unterhalb der Colorbar anzeigen
                #font=dict(color='black', size=12)
            )
        )
    ))

    # X-Achse (Spalten) oben
    fig.update_xaxes(
        side='top',
        tickmode='array',
        tickvals=xvals,
        ticktext=[str(col) for col in df.columns],
        range=[0, n_cols],  # erlaubt Zeichnen von Linien bis n_cols
        showgrid=False,
        zeroline=False,
        showline=False,
        scaleanchor='y',  # 1:1-Seitenverhältnis
        scaleratio=1
    )

    # Y-Achse (Zeilen) links, von oben nach unten
    fig.update_yaxes(
        tickmode='array',
        tickvals=yvals,
        ticktext=[str(idx) for idx in df.index],
        range=[0, n_rows],
        showgrid=False,
        zeroline=False,
        showline=False,
        autorange='reversed'  # oberste Zeile = df.index[0]
    )

    # Dünne Linien um jede Zelle (Shapes)
    # => Grid aus horizontalen und vertikalen Linien
    for i in range(n_rows + 1):
        fig.add_shape(
            type='line',
            x0=0, y0=i,
            x1=n_cols, y1=i,
            line=dict(color=margin_color, width=margin_width),
            xref='x', yref='y',
            layer='above'  # Linie über der Heatmap
        )
    for j in range(n_cols + 1):
        fig.add_shape(
            type='line',
            x0=j, y0=0,
            x1=j, y1=n_rows,
            line=dict(color=margin_color, width=margin_width),
            xref='x', yref='y',
            layer='above'
        )

    # Kompakte Abmessungen:
    cell_size_px = 70  # Pixel pro Zelle
    top_margin = 80
    bottom_margin = 50
    left_margin = 90
    right_margin = 30
    total_width = n_cols * cell_size_px + left_margin + right_margin
    total_height = n_rows * cell_size_px + top_margin + bottom_margin

    fig.update_layout(
        autosize=False,
        width=total_width,
        height=total_height,
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=bottom_margin),
        #plot_bgcolor='white'
    )

    return fig

#!/usr/bin/env python3

import pandas as pd
import numpy as np
from difflib import SequenceMatcher
import plotly.graph_objects as go


def get_sample_information(df):
    """
    Extracts sample information from the DataFrame and returns it as a dictionary.

    Parameters:
        df (pd.DataFrame): DataFrame containing sample information with MultiIndex columns
                           (e.g., ('SAMPLE_INFORMATION', 'COHORT'), ('SAMPLE_INFORMATION', 'SUBJECT'), etc.).

    Returns:
        dict: A dictionary with sample information including Cohort, Donor, Time points, Sample, Subset, Total BCRs,
              Heavy chains, Light chains, Kappa chains, and Lambda chains.
    """
    infos = {
        'Cohort': ', '.join(df[('SAMPLE_INFORMATION', 'COHORT')].unique()),
        'Donor': ', '.join(df[('SAMPLE_INFORMATION', 'SUBJECT')].unique()),
        'Time points': ', '.join(df[('SAMPLE_INFORMATION', 'TIME_POINT')].unique()),
        'Sample': ', '.join(df[('SAMPLE_INFORMATION', 'TISSUE')].unique()),
        'Subset': ', '.join(df[('SAMPLE_INFORMATION', 'SUBSET')].unique()),
        'Total BCRs': df.shape[0],
        'Heavy chains': int(df[('HEAVY_CHAIN', 'HEAVY_FOUND')].sum()),
        'Light chains': int(df[('LIGHT_CHAIN', 'LIGHT_FOUND')].sum()),
        'Kappa chains': int((df[('LIGHT_CHAIN', 'PCR_ISOTYPE')].dropna() == "KC").sum()),
        'Lambda chains': int((df[('LIGHT_CHAIN', 'PCR_ISOTYPE')].dropna() == "LC").sum())
    }
    return infos


def get_hydrophobicity_values(sequence, scale='kyte-doolittle', ignore_unknown=True, warn=True):
    """
    Returns a list of hydrophobicity values for the given amino acid sequence based on the selected scale.
    Unknown amino acids (e.g., 'X') are optionally ignored.

    Parameters:
        sequence (str): Amino acid sequence (using single-letter codes).
        scale (str): Hydrophobicity scale to use ('kyte-doolittle' or 'eisenberg'; default is 'kyte-doolittle').
        ignore_unknown (bool): If True, unknown residues are skipped. If False, raises a ValueError.
        warn (bool): If True, prints a warning when unknown residues are ignored.

    Returns:
        list: A list of hydrophobicity values corresponding to each known residue in the sequence.

    Raises:
        ValueError: If ignore_unknown=False and an unknown residue is encountered.
    """
    kyte_doolittle = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
        'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    eisenberg = {
        'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
        'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
        'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
        'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
    }

    scale_dict = {
        'kyte-doolittle': kyte_doolittle,
        'eisenberg': eisenberg
    }.get(scale.lower())

    if scale_dict is None:
        raise ValueError("Unsupported scale. Please choose 'kyte-doolittle' or 'eisenberg'.")

    values = []
    for aa in sequence.upper():
        if aa in scale_dict:
            values.append(scale_dict[aa])
        elif ignore_unknown:
            if warn:
                print(f"Warning: Unknown amino acid '{aa}' ignored.")
            continue
        else:
            raise ValueError(f"Unknown amino acid '{aa}' in sequence.")

    return values


def compute_gravy_scores(sequences: pd.Series, scale: str = 'kyte-doolittle') -> pd.Series:
    """
    Computes the GRAVY (Grand Average of Hydropathy) scores for a Series of CDR3 amino acid sequences.

    The GRAVY score is calculated as the average of the hydrophobicity values for the amino acids in the sequence.

    Parameters:
        sequences (pd.Series): Series of CDR3 amino acid sequences (using single-letter codes).
        scale (str): Hydrophobicity scale to use ('kyte-doolittle' or 'eisenberg'; default is 'kyte-doolittle').

    Returns:
        pd.Series: Series containing the GRAVY scores for each sequence.
    """
    def compute_gravy(seq: str) -> float:
        values = get_hydrophobicity_values(seq, scale=scale)
        if len(values) == 0:
            return None
        return sum(values) / len(values)
    return sequences.apply(compute_gravy)


def compute_histogram(series: pd.Series, bin_size: float = 1, min_val: float = None, max_val: float = None) -> pd.Series:
    """
    Computes a histogram from a numeric Pandas Series using the specified bin size.

    Parameters:
        series (pd.Series): Numeric data as a Pandas Series.
        bin_size (float): Size of each histogram bin (default is 1).
        min_val (float): Minimum value for the bins (if None, uses the series minimum).
        max_val (float): Maximum value for the bins (if None, uses the series maximum).

    Returns:
        pd.Series: A Pandas Series where the index represents the start of each bin and the values are the counts.
    """
    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError("The input series must be numeric.")
    if min_val is None:
        min_val = series.min()
    if max_val is None:
        max_val = series.max()
    bins = np.arange(min_val, max_val + bin_size, bin_size)
    counts, bin_edges = np.histogram(series, bins=bins)
    index = bin_edges[:-1]
    hist_series = pd.Series(counts, index=index)
    return hist_series


def compute_diversity_index(sequences: pd.Series, index_type: str = "shannon") -> float:
    """
    Computes a diversity index for a Series of CDR3 amino acid sequences.

    The diversity index is computed from the frequency distribution of unique sequences:
      - Shannon index: H = -∑ (p_i * ln(p_i))
      - Inverse Simpson index: D = 1 / ∑ (p_i^2)

    Parameters:
        sequences (pd.Series): Series containing CDR3 sequences.
        index_type (str): The type of diversity index to compute ("shannon" or "inverse_simpson").

    Returns:
        float: The computed diversity index.

    Raises:
        ValueError: If an unsupported index_type is provided.
    """
    counts = sequences.value_counts()
    total = counts.sum()
    probs = counts / total
    if index_type.lower() == "shannon":
        shannon_index = -np.sum(probs * np.log(probs))
        return shannon_index
    elif index_type.lower() == "inverse_simpson":
        inverse_simpson = 1.0 / np.sum(probs ** 2)
        return inverse_simpson
    else:
        raise ValueError("index_type must be either 'shannon' or 'inverse_simpson'.")


def subsample_series_dict(series_dict: dict) -> tuple:
    """
    Determines the smallest Series (by length) from a dictionary of Series and returns:
      - the original dictionary,
      - the minimal size,
      - and the key corresponding to the smallest Series.

    Parameters:
        series_dict (dict): Dictionary with dataset names as keys and Pandas Series as values.

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
    Computes diversity indices for each dataset in a dictionary of Series via subsampling.

    For each dataset in the dictionary, if the Series is larger than the smallest Series, multiple subsamples
    (of size equal to the smallest Series) are drawn and the diversity index is computed for each subsample.
    For the smallest Series, the diversity index is computed directly.

    Parameters:
        series_dict (dict): Dictionary where keys are dataset names and values are Pandas Series.
        index_type (str): Type of diversity index to compute ("shannon" or "inverse_simpson").
        iterations (int): Number of subsampling iterations (default is 20).
        random_state (int): Seed for reproducibility (default is 42).
        seq_column (str): The name to assign to the Series for clarity (not used internally).

    Returns:
        dict: A dictionary with dataset names as keys and a dict as value containing:
              'n' (minimal sample size), 'mean' (mean diversity index), and 'std' (standard deviation).
    """
    series_dict, min_size, smallest_key = subsample_series_dict(series_dict)
    results = {}
    for key, series in series_dict.items():
        series = series.rename(seq_column)
        if series.shape[0] == min_size:
            diversity_value = compute_diversity_index(series, index_type=index_type)
            results[key] = {'n': min_size, 'mean': diversity_value, 'std': 0.0}
        else:
            indices = []
            for i in range(iterations):
                subsample = series.sample(n=min_size, random_state=random_state + i)
                index_value = compute_diversity_index(subsample, index_type=index_type)
                indices.append(index_value)
            mean_val = np.mean(indices)
            std_val = np.std(indices)
            results[key] = {'n': min_size, 'mean': mean_val, 'std': std_val}
    return results


def compute_pairwise_jaccard(data_dict: dict, iterations: int = 20, random_state: int = 42):
    """
    Computes the pairwise Jaccard index for multiple sequence datasets via subsampling.

    For each dataset in the dictionary, unique sequences are extracted and the minimal set size is determined.
    For a specified number of iterations, subsamples of size equal to the minimal set size are drawn from each dataset.
    The Jaccard index for each pair of datasets is computed as:
        J(A, B) = |A ∩ B| / |A ∪ B|

    Parameters:
        data_dict (dict): Dictionary where keys are dataset names and values are lists or Pandas Series of sequences.
        iterations (int): Number of subsampling iterations (default is 20).
        random_state (int): Seed for reproducibility (default is 42).

    Returns:
        tuple: A tuple (results, matrix) where:
               - results (dict): Dictionary with keys as tuples (dataset1, dataset2) and values as a dict with:
                                 'n' (minimal set size), 'mean' (mean Jaccard index), and 'std' (standard deviation).
               - matrix (pd.DataFrame): A symmetric DataFrame with dataset names as rows and columns containing the mean Jaccard index.
    """
    unique_dict = {}
    for dataset, seqs in data_dict.items():
        if isinstance(seqs, pd.Series):
            seq_list = seqs.dropna().unique().tolist()
        else:
            seq_list = list(set(seqs))
        unique_dict[dataset] = seq_list

    min_size = min(len(seq_list) for seq_list in unique_dict.values())
    dataset_names = sorted(list(unique_dict.keys()))
    pairwise_results = {}
    for i in range(len(dataset_names)):
        for j in range(i, len(dataset_names)):
            pairwise_results[(dataset_names[i], dataset_names[j])] = []

    rng = np.random.RandomState(random_state)
    for _ in range(iterations):
        subsamples = {}
        for dataset, seq_list in unique_dict.items():
            subsample = rng.choice(seq_list, size=min_size, replace=False)
            subsamples[dataset] = set(subsample)
        for i in range(len(dataset_names)):
            for j in range(i, len(dataset_names)):
                setA = subsamples[dataset_names[i]]
                setB = subsamples[dataset_names[j]]
                union = setA | setB
                intersection = setA & setB
                jacc = len(intersection) / len(union) if len(union) > 0 else 0
                pairwise_results[(dataset_names[i], dataset_names[j])].append(jacc)
    results = {}
    for pair, values in pairwise_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        results[pair] = {'n': min_size, 'mean': mean_val, 'std': std_val}
    matrix = pd.DataFrame(index=dataset_names, columns=dataset_names, dtype=float)
    for i in range(len(dataset_names)):
        for j in range(len(dataset_names)):
            key = tuple(sorted((dataset_names[i], dataset_names[j])))
            matrix.iloc[i, j] = results[key]['mean']
    return results, matrix


def sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculates the similarity between two sequences using difflib's SequenceMatcher.

    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        float: A similarity ratio between 0 (no similarity) and 1 (identical).
    """
    return SequenceMatcher(None, seq1, seq2).ratio()


def compute_pairwise_weighted_jaccard(data_dict: dict, iterations: int = 20, random_state: int = 42):
    """
    Computes the weighted pairwise Jaccard index between multiple sequence datasets.

    For each dataset, unique sequences are extracted and the minimal set size is determined.
    For a specified number of iterations, subsamples of size equal to the minimal set size are drawn without replacement.
    For each pair, the weighted Jaccard index is computed as:
        weighted J(A, B) = soft_intersection / soft_union,
    where:
        soft_intersection = sum over a in A of (max_{b in B} similarity(a, b))
        soft_union = |A| + |B| - soft_intersection

    Parameters:
        data_dict (dict): Dictionary where keys are dataset names and values are lists or Pandas Series of sequences.
        iterations (int): Number of subsampling iterations (default is 20).
        random_state (int): Seed for reproducibility (default is 42).

    Returns:
        tuple: A tuple (results, matrix) where:
               - results (dict): Dictionary with keys as tuples (dataset1, dataset2) and values as a dict with:
                                 'n' (minimal set size), 'mean' (mean weighted Jaccard index), and 'std' (standard deviation).
               - matrix (pd.DataFrame): A symmetric DataFrame with dataset names as rows and columns containing the mean weighted Jaccard index.
    """
    unique_dict = {}
    for dataset, seqs in data_dict.items():
        if isinstance(seqs, pd.Series):
            seq_list = seqs.dropna().unique().tolist()
        else:
            seq_list = list(set(seqs))
        unique_dict[dataset] = seq_list

    min_size = min(len(seq_list) for seq_list in unique_dict.values())
    dataset_names = sorted(list(unique_dict.keys()))
    pairwise_results = {}
    for i in range(len(dataset_names)):
        for j in range(i, len(dataset_names)):
            pairwise_results[(dataset_names[i], dataset_names[j])] = []

    rng = np.random.RandomState(random_state)
    for _ in range(iterations):
        subsamples = {}
        for dataset, seq_list in unique_dict.items():
            subsample = rng.choice(seq_list, size=min_size, replace=False)
            subsamples[dataset] = set(subsample)
        for i in range(len(dataset_names)):
            for j in range(i, len(dataset_names)):
                setA = subsamples[dataset_names[i]]
                setB = subsamples[dataset_names[j]]
                soft_intersection = 0.0
                for a in setA:
                    best_sim = max(sequence_similarity(a, b) for b in setB)
                    soft_intersection += best_sim
                union_score = len(setA) + len(setB) - soft_intersection
                jacc = soft_intersection / union_score if union_score > 0 else 0
                pairwise_results[(dataset_names[i], dataset_names[j])].append(jacc)
    results = {}
    for pair, values in pairwise_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        results[pair] = {'n': min_size, 'mean': mean_val, 'std': std_val}
    matrix = pd.DataFrame(index=dataset_names, columns=dataset_names, dtype=float)
    for i in range(len(dataset_names)):
        for j in range(len(dataset_names)):
            key = tuple(sorted((dataset_names[i], dataset_names[j])))
            matrix.iloc[i, j] = results[key]['mean']
    return results, matrix

#=====================
# PLOTTING FUNCTIONS
#=====================

def create_heatmap_with_annotations(
        df: pd.DataFrame,
        colorscale: str = 'reds',
        margin_color='black',
        margin_width=1,
        color_bar_title="Colorbar"
) -> go.Figure:
    """
    Creates a heatmap with custom annotations for columns and rows.

    This function generates a Plotly heatmap from the values in a DataFrame. The DataFrame's columns
    and index are used as the labels for the x-axis and y-axis respectively. Additionally, the heatmap
    includes a horizontal colorbar with a custom title, and grid lines are drawn around each cell.

    Parameters:
        df (pd.DataFrame): DataFrame containing the numeric values for the heatmap. The DataFrame's columns
                           and index are used as labels.
        colorscale (str): Colorscale for the heatmap (default: 'reds').
        margin_color (str): Color used for the borders around cells and the colorbar outline.
        margin_width (int): Width of the borders around cells and the colorbar outline.
        color_bar_title (str): Title for the colorbar.

    Returns:
        go.Figure: A Plotly Figure object containing the heatmap with annotations.
    """
    n_rows, n_cols = df.shape

    # Calculate cell center coordinates: x = 0.5 .. (n_cols - 0.5), y = 0.5 .. (n_rows - 0.5)
    xvals = np.arange(n_cols) + 0.5
    yvals = np.arange(n_rows) + 0.5

    # Create the heatmap (using numeric axes) with a horizontal colorbar.
    fig = go.Figure(data=go.Heatmap(
        x=xvals,
        y=yvals,
        z=df.values,
        colorscale=colorscale,
        colorbar=dict(
            orientation='h',           # horizontal colorbar
            x=0,                       # left-aligned
            y=-0.1,                    # positioned closer to the heatmap
            xanchor='left',
            yanchor='top',
            thickness=10,              # thinner colorbar
            len=1,                     # spans the entire width of the heatmap domain
            outlinecolor=margin_color, # color of the colorbar border
            outlinewidth=margin_width, # width of the colorbar border
            title=dict(
                text=color_bar_title,  # text for the colorbar title
                side='bottom'          # display the title below the colorbar
            )
        )
    ))

    # Configure x-axis: place column labels at the top.
    fig.update_xaxes(
        side='top',
        tickmode='array',
        tickvals=xvals,
        ticktext=[str(col) for col in df.columns],
        range=[0, n_cols],  # allow drawing lines up to n_cols
        showgrid=False,
        zeroline=False,
        showline=False,
        scaleanchor='y',    # maintain a 1:1 aspect ratio
        scaleratio=1
    )

    # Configure y-axis: place row labels on the left and reverse the order so that the first row is at the top.
    fig.update_yaxes(
        tickmode='array',
        tickvals=yvals,
        ticktext=[str(idx) for idx in df.index],
        range=[0, n_rows],
        showgrid=False,
        zeroline=False,
        showline=False,
        autorange='reversed'  # ensures the top row corresponds to the first index
    )

    # Draw thin grid lines around each cell.
    # Horizontal lines.
    for i in range(n_rows + 1):
        fig.add_shape(
            type='line',
            x0=0, y0=i,
            x1=n_cols, y1=i,
            line=dict(color=margin_color, width=margin_width),
            xref='x', yref='y',
            layer='above'  # draw the line above the heatmap
        )
    # Vertical lines.
    for j in range(n_cols + 1):
        fig.add_shape(
            type='line',
            x0=j, y0=0,
            x1=j, y1=n_rows,
            line=dict(color=margin_color, width=margin_width),
            xref='x', yref='y',
            layer='above'
        )

    # Set compact layout dimensions.
    cell_size_px = 70  # pixels per cell
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
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=bottom_margin)
    )

    return fig

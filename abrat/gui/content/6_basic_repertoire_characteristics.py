# gui/content/6_basic_repertoire_characteristics.py

import os, glob, math
import time

import streamlit as st
import pandas as pd
import io, zipfile
import datetime


import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import altair as alt

import numpy as np
import peptides

from abrat.gui.gui_shared import folder_selectbox
from abrat.core.repertoire_characteristics import (compute_gravy_scores,
                                                   get_sample_information,
                                                   compute_diversity_indices_for_subsampling,
                                                   compute_histogram,
                                                   )

# FUNCTIONS: TODO: move to core, when set up

## Standard Columns

sample_col = 'SAMPLE_INFORMATION'
hc_col = 'HEAVY_CHAIN'
lc_col = 'LIGHT_CHAIN'

v_gene_subcol = 'V_GENE'
d_gene_subcol = 'D_GENE'
j_gene_subcol = 'J_GENE'

top_v_gene_subcol = 'TOP_V'
top_d_gene_subcol = 'TOP_D'
top_j_gene_subcol = 'TOP_J'

cdr3_aa_subcol = 'CDR3_AA'
cdr3_aa_len_subcol = 'CDR3_AA_LENGTH'
cdr3_nt_subcol = 'CDR3_NT'

v_ident_subcol = 'V_IDENTITY'

hc_cluster_subcol = 'HC_SUBCLUSTER'
lc_cluster_subcol = 'LC_SUBCLUSTER'
representative_bcr_col = 'CLUSTER_REPRESENTATIVE'
hclc_cluster_subcol = 'HC-LC_CLUSTER'
cluster_size_subcol = 'CLUSTER_SIZE'
clonal_subcol = 'IS_CLONAL'
clone_subcol = 'CLONE'
clone_color_subcol = 'CLONE_COLOR'

marker_list = ["circle-open", "square-open", "diamond-open", "cross", "x",
               "triangle-up-open", "triangle-down-open", "triangle-left-open", "triangle-right-open"]

def collapse_clone_data(df, target_column, clone_col=None, mode='unique_values'):
    """
    Collapses a DataFrame by reducing entries in a specified target column, optionally grouping by clone.

    If a clone column is provided and exists in the DataFrame, the function groups the DataFrame by that column.
    Within each group, it either computes the mean (if mode is 'mean') or retains only the unique values
    (if mode is 'unique_values') of the target column. If no clone column is provided, the function simply drops
    duplicate values in the target column.

    Parameters:
        df (pd.DataFrame): DataFrame containing gene segment data.
        target_column (tuple): A tuple representing the target column (e.g., ('HEAVY_CHAIN', 'TOP_V')).
        clone_col (tuple, optional): A tuple representing the clone identifier column (e.g., ('SAMPLE_INFORMATION', 'hclc_cluster_subcol')).
                                     If provided, collapse is performed within each clone. Default is None.
        mode (str, optional): Collapse mode. 'unique_values' (default) retains only unique values,
                              'mean' computes the mean value for the group.

    Returns:
        pd.DataFrame: A reduced DataFrame containing the collapsed values from the target column.
    """
    if clone_col is not None and clone_col in df.columns:
        if mode == 'mean':
            # Group by the clone identifier and compute the mean for the target column.
            collapsed = df.groupby(clone_col)[[target_column]].apply(lambda x: x.mean())
        else:
            # Group by the clone identifier and drop duplicate values for the target column.
            collapsed = df.groupby(clone_col)[[target_column]].apply(lambda x: x.drop_duplicates())

        # Reset index and remove the clone column from the result.
        collapsed_df = collapsed.reset_index().drop(columns=[clone_col])
    else:
        # If no clone column is provided, simply drop duplicate values in the target column.
        collapsed_df = df.drop_duplicates(subset=[target_column])

    return collapsed_df

def get_repertoire_statistics(df):
    """
    Computes basic repertoire statistics from the given DataFrame.

    This function calculates various statistics from the BCR repertoire data which can later be
    displayed on the dashboard. The statistics include clonality metrics (such as counts of undefined,
    non-clonal, and clonal sequences; number of clones; mean and median clone sizes; clone sizes; and
    clone color mapping), as well as V(D)J gene segment, isotype, and V gene identity statistics for
    heavy and light chains. Additionally, it computes CDR3-related statistics including length distribution,
    hydrophobicity (using both Eisenberg and Kyte-Doolittle scales), and net charge.

    Parameters:
        df (pd.DataFrame): A pandas DataFrame containing the repertoire data with MultiIndex columns.
                           Expected keys include:
                           - Sample information columns: sample_col, clone_subcol, clone_color_subcol,
                             cluster_size_subcol, and hclc_cluster_subcol.
                           - Heavy chain columns (hc_col) and light chain columns (lc_col) with various sub-columns
                             such as TOP_V, TOP_D, TOP_J, TOP_ISOTYPE (or PCR_ISOTYPE), V_IDENT (v_ident_subcol),
                             CDR3_AA (cdr3_aa_subcol), and CDR3_AA_LENGTH (cdr3_aa_len_subcol).

    Returns:
        dict: A dictionary with the following keys:
              - 'clonality': A dictionary with clonality statistics, including:
                  • 'undefined': Count of undefined clone entries.
                  • 'non_clonal': Count of non-clonal sequences.
                  • 'clonal': Number of clonal sequences.
                  • 'number_of_clones': Number of unique clones.
                  • 'mean_clone_size': Mean size of the clones.
                  • 'median_clone_size': Median size of the clones.
                  • 'clone_sizes': DataFrame of clone sizes.
                  • 'clone_colors': Mapping of clone identifiers to their assigned colors.
              - For each chain type ('heavy_chain', 'light_chain', 'kappa_light_chain', 'lambda_light_chain'):
                  For both 'not_collapsed' and 'collapsed' modes:
                    • Gene segment statistics for each gene in the chain (e.g., TOP_V, TOP_D, TOP_J, TOP_ISOTYPE)
                      including counts and percentages (with and without missing values).
                    • V gene identity statistics including histograms, mean, standard deviation, geometric mean,
                      and geometric standard deviation.
              - Additionally, for each chain type and for both 'all' and 'unique' modes:
                  • CDR3-related statistics including the CDR3 amino acid sequences, length distribution, hydrophobicity
                    (for both Kyte-Doolittle and Eisenberg scales), and net charge.

    Note:
        This function depends on external variables (e.g., sample_col, clone_subcol, clone_color_subcol, cluster_size_subcol,
        hc_col, lc_col, v_ident_subcol, cdr3_aa_subcol, cdr3_aa_len_subcol, hclc_cluster_subcol) and helper functions such as
        collapse_clone_data, compute_histogram, and compute_gravy_scores.
    """
    stats = {}

    # ----------
    # Clonality
    # ----------
    clone_df = df[df[(sample_col, clonal_subcol)].astype('boolean')]
    color_dict = (
        df[[(sample_col, clone_subcol), (sample_col, clone_color_subcol)]]
        .drop_duplicates(subset=[(sample_col, clone_color_subcol)])
        .set_index((sample_col, clone_subcol))[(sample_col, clone_color_subcol)]
        .to_dict()
    )

    clone_order = df[(sample_col, clone_subcol)].fillna('Undefined').drop_duplicates().to_list()
    clone_sizes = df[(sample_col, clone_subcol)].value_counts().reindex(clone_order).reset_index()

    stats['clonality'] = {
        'undefined': sum(df[(sample_col, clone_subcol)] == 'Undefined') + sum(df[(sample_col, clone_subcol)].isnull()),
        'non_clonal': sum(df[(sample_col, clone_subcol)] == 'Non-clonal'),
        'clonal': len(clone_df),
        'number_of_clones': len(clone_df[(sample_col, clone_subcol)].unique()),
        'mean_clone_size': clone_df.drop_duplicates((sample_col, clone_subcol))[
            (sample_col, cluster_size_subcol)].mean(),
        'median_clone_size': clone_df.drop_duplicates((sample_col, clone_subcol))[
            (sample_col, cluster_size_subcol)].median(),
        'clone_sizes': clone_sizes,
        'clone_colors': color_dict
    }

    # ------------------------------------------
    # V(D)J Gene Segments, Isotype, and V Gene Identities
    # ------------------------------------------
    chain_genes = {
        'heavy_chain': ['TOP_V', 'TOP_D', 'TOP_J', 'TOP_ISOTYPE'],
        'light_chain': ['TOP_V', 'TOP_J', 'PCR_ISOTYPE'],
        'kappa_light_chain': ['TOP_V', 'TOP_J', 'TOP_ISOTYPE'],
        'lambda_light_chain': ['TOP_V', 'TOP_J', 'TOP_ISOTYPE']
    }

    # Use composite cluster as the clone identifier.
    clone_identifier = (sample_col, hclc_cluster_subcol)

    for mode in ['not_collapsed', 'collapsed']:
        for chain, gene_list in chain_genes.items():
            # Initialize nested dictionary for the chain and mode.
            stats.setdefault(chain, {})[mode] = {}

            # Determine the appropriate DataFrame subset based on chain type.
            if chain == 'heavy_chain':
                chain_key = hc_col
                working_df = df.copy()
            elif chain in ['light_chain', 'kappa_light_chain', 'lambda_light_chain']:
                chain_key = lc_col
                if chain == 'kappa_light_chain':
                    working_df = df[df[(chain_key, 'PCR_ISOTYPE')] != "LC"]
                elif chain == 'lambda_light_chain':
                    working_df = df[df[(chain_key, 'PCR_ISOTYPE')] != "KC"]
                else:
                    working_df = df.copy()
            else:
                chain_key = chain.upper()
                working_df = df.copy()

            for gene in gene_list:
                v_gene_col = (chain_key, gene)
                if mode == 'collapsed':
                    working_df_mode = collapse_clone_data(working_df, target_column=v_gene_col, clone_col=clone_identifier).copy()
                else:
                    working_df_mode = working_df

                # Calculate counts and percentages including NaNs.
                with_nans_count = working_df_mode[v_gene_col].fillna('N.D.').value_counts()
                with_nans_percent = 100 * with_nans_count / with_nans_count.sum()

                # Calculate counts and percentages excluding NaNs.
                no_nans_count = working_df_mode[v_gene_col].replace('N.D.', np.nan).value_counts()
                no_nans_percent = 100 * no_nans_count / no_nans_count.sum()

                stats[chain][mode][gene.lower()] = {
                    'with_nans': {'counts': with_nans_count, 'percent': with_nans_percent},
                    'no_nans': {'counts': no_nans_count, 'percent': no_nans_percent}
                }

            # V gene identity.
            v_ident_col = (chain_key, v_ident_subcol)
            if mode == 'collapsed':
                working_df_vident = collapse_clone_data(working_df, target_column=v_ident_col, clone_col=clone_identifier).copy()
            else:
                working_df_vident = working_df

            # Calculate histograms for V gene identity; NaNs are replaced by 0 for one histogram.
            v_with_nans_count = compute_histogram(working_df_vident[v_ident_col].fillna(0), bin_size=1, min_val=-0, max_val=101)
            v_with_nans_percent = 100 * v_with_nans_count / v_with_nans_count.sum()

            v_no_nans_count = compute_histogram(working_df_vident[v_ident_col], bin_size=1, min_val=0, max_val=101)
            v_no_nans_percent = 100 * v_no_nans_count / v_no_nans_count.sum()

            stats[chain][mode]['v_identity'] = {
                'with_nans': {'counts': v_with_nans_count, 'percent': v_with_nans_percent},
                'no_nans': {
                    'identities': working_df_vident[v_ident_col].dropna(),
                    'counts': v_no_nans_count,
                    'percent': v_no_nans_percent,
                    'mean': working_df_vident[v_ident_col].dropna().mean(),
                    'std': working_df_vident[v_ident_col].dropna().std(),
                    'geomean': np.exp(np.log(working_df_vident[v_ident_col].dropna()).mean()),
                    'geostd': np.exp(np.log(working_df_vident[v_ident_col].dropna()).std())
                }
            }

    # -------------------
    # CDR3 and Mutations
    # -------------------
    for chain in chain_genes:
        for mode in ['all', 'unique']:
            # Initialize nested dictionary for the chain and mode.
            stats.setdefault(chain, {})[mode] = {}
            if chain == 'heavy_chain':
                chain_key = hc_col
                working_df = df.copy()
            elif chain in ['light_chain', 'kappa_light_chain', 'lambda_light_chain']:
                chain_key = lc_col
                if chain == 'kappa_light_chain':
                    working_df = df[df[(chain_key, 'PCR_ISOTYPE')] != "LC"]
                elif chain == 'lambda_light_chain':
                    working_df = df[df[(chain_key, 'PCR_ISOTYPE')] != "KC"]
                else:
                    working_df = df.copy()
            else:
                chain_key = chain.upper()
                working_df = df.copy()

            cdr3_aa_col = (chain_key, cdr3_aa_subcol)
            cdr3_aa_len_col = (chain_key, cdr3_aa_len_subcol)

            if mode == 'unique':
                working_df_mode = working_df.drop_duplicates(cdr3_aa_col)
            else:
                working_df_mode = working_df

            cdr3s = working_df_mode[cdr3_aa_col].dropna()
            cdr3_length_counts = working_df_mode[cdr3_aa_len_col].dropna().value_counts()
            cdr3_length_percent = 100 * cdr3_length_counts / cdr3_length_counts.sum()
            cdr3_gravy_eisenberg = compute_gravy_scores(cdr3s, 'eisenberg')
            cdr3_gravy_kyte = compute_gravy_scores(cdr3s, 'kyte-doolittle')
            cdr3_net_charges = cdr3s.apply(lambda seq: peptides.Peptide(seq).charge(pH=7.4))

            stats[chain][mode]['CDR3_AA'] = {
                'cdr3s': cdr3s,
                'cdr3_length_distribution': {'counts': cdr3_length_counts, 'percent': cdr3_length_percent},
                'cdr3_hydrophobicity': {'kyte-doolittle': cdr3_gravy_kyte, 'eisenberg': cdr3_gravy_eisenberg},
                'cdr3_net_charges': cdr3_net_charges
            }

    return stats

def get_gene_segment_stats(data_dict, chain, gene_segment, collapsed='not_collapsed', nans='no_nans', output='percent'):
    """
    Generates a DataFrame with gene segment statistics from the given data dictionary.

    The data_dict is expected to be a dictionary where each key corresponds to a dataset (e.g., a filename)
    and contains an 'alias' and a nested 'statistics' dictionary. This function extracts the specified gene segment
    statistic for the given chain from each dataset, using the provided parameters to determine whether to use
    collapsed data and whether to include missing values (NaNs). The output statistic can be chosen (e.g., 'percent' or 'counts').

    Parameters:
        data_dict (dict): Dictionary of datasets with their respective statistics.
        chain (str): The chain type to extract statistics for (e.g., 'heavy_chain', 'light_chain').
        gene_segment (str): The gene segment statistic to extract (e.g., 'top_v').
        collapsed (str, optional): Either 'not_collapsed' or 'collapsed', indicating whether to use collapsed data.
                                   Default is 'not_collapsed'.
        nans (str, optional): Specifies whether to use values including NaNs ('with_nans') or excluding them ('no_nans').
                              Default is 'no_nans'.
        output (str, optional): The output statistic to retrieve, such as 'percent' or 'counts'. Default is 'percent'.

    Returns:
        pd.DataFrame: A DataFrame containing the extracted gene segment statistics, with dataset aliases as columns.
    """
    gene_segment_data = {}
    for dataset in data_dict.keys():
        dataset_alias = data_dict[dataset]['alias']
        gene_segment_data[dataset_alias] = data_dict[dataset]['statistics'][chain][collapsed][gene_segment][nans][output]
    gene_segment_df = pd.DataFrame(data=gene_segment_data).fillna(0).sort_index()
    return gene_segment_df

def update_aliases():
    """
    Callback function that is executed immediately when the data editor is modified.

    This function reads the edited aliases from the 'edited_aliases_df' DataFrame in the session state,
    creates a dictionary mapping file names to aliases, and updates the corresponding alias in the
    'clustered_bcrs' section of 'bcr_results' in the session state.
    """
    alias_dict = st.session_state.edited_aliases_df.set_index('File Name')['Alias'].to_dict()
    for file_name, alias in alias_dict.items():
        st.session_state.bcr_results['clustered_bcrs'][file_name]['alias'] = alias

### PLOT FUNCTIONS
def plot_grouped_bar_chart(df, color_dict, x_axis_title="X-Axis", y_axis_title="Y-Axis"):
    """
    Plots a grouped bar chart using Plotly.

    The DataFrame's index is used as the x-values, and each column is plotted as a separate dataset,
    with the column names used as legend labels. The chart is set to have an approximately 3:1 width-to-height ratio,
    and the legend is horizontally centered below the chart.

    Parameters:
        df (pd.DataFrame): DataFrame where the index represents x-values and the columns contain y-values.
        color_dict (dict): Dictionary mapping each column name to a base color code.
        x_axis_title (str): Title for the x-axis.
        y_axis_title (str): Title for the y-axis.

    Returns:
        go.Figure: A Plotly Figure object representing the grouped bar chart.
    """
    traces = []
    # Create a bar trace for each column in the DataFrame.
    for col in df.columns:
        base_color = color_dict.get(col, "#000000")
        fill_color = hex_to_rgba(base_color, alpha=0.6)
        traces.append(go.Bar(
            x=df.index,
            y=df[col],
            name=str(col),
            marker=dict(
                color=fill_color,
                line=dict(color=base_color, width=0.5)
            )
        ))

    # Create the figure with the generated traces.
    fig = go.Figure(data=traces)

    # Update the layout: use group mode, set figure size, and position the legend below the chart.
    fig.update_layout(
        barmode='group',
        width=900,  # Approximately 3:1 aspect ratio.
        height=300,
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-1.1  # Place the legend below the plot.
        ),
        margin=dict(l=50, r=50, t=30, b=100),
        xaxis_title=x_axis_title,
        yaxis_title=y_axis_title
    )

    # Limit the number of y-axis ticks.
    fig.update_yaxes(nticks=6)

    return fig


def hex_to_rgba(hex_color, alpha=0.5):
    """
    Converts a hex color string to an rgba string.

    Parameters:
        hex_color (str): The color in hex format (e.g., "#FF5733").
        alpha (float): The alpha (opacity) value (default is 0.5).

    Returns:
        str: The color represented in rgba format (e.g., "rgba(255, 87, 51, 0.5)").
    """
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return f'rgba({r}, {g}, {b}, {alpha})'


def plot_interactive_line_plot(df, color_dict, x_axis_title="X-Axis", y_axis_title="Y-Axis",
                               fill_plot=False,
                               marker_list=None,
                               x_axis_minor_ticks=5,
                               x_axis_range=None,
                               title=''):
    """
    Creates an interactive line plot using Plotly.

    Each column in the DataFrame is plotted as a separate line. The colors for each line are assigned
    using the provided color_dict. Optionally, if fill_plot is True, the area under each line is filled with 50%
    opacity. Marker symbols can be specified via marker_list; if there are fewer markers than columns, the list
    is cycled through.

    Parameters:
        df (pd.DataFrame): DataFrame where the index represents x-values and columns represent y-values.
        color_dict (dict): Dictionary mapping column names to color codes (hex format).
        x_axis_title (str): Title for the x-axis.
        y_axis_title (str): Title for the y-axis.
        fill_plot (bool): If True, fills the area under the lines with a semi-transparent color (default: False).
        marker_list (list): List of marker symbols (e.g., ['circle', 'square', 'diamond']). Default is None.
        x_axis_minor_ticks (int): Number of minor ticks to display on the x-axis.
        x_axis_range (list or tuple): Range for the x-axis (e.g., [min, max]). Default is None.
        title (str): The chart title.

    Returns:
        go.Figure: A Plotly Figure object representing the interactive line plot.
    """
    traces = []

    # Create a trace for each column in the DataFrame.
    for i, col in enumerate(df.columns):
        color = color_dict.get(col, "#000000")
        line_color = hex_to_rgba(color, 0.8)
        # Determine the marker symbol if a marker list is provided.
        marker_symbol = marker_list[i % len(marker_list)] if marker_list else None

        # Set the mode to include markers if a marker symbol is provided.
        mode = "lines+markers" if marker_symbol is not None else "lines"

        # Create the scatter trace.
        trace = go.Scatter(
            x=df.index,
            y=df[col],
            name=str(col),
            mode=mode,
            line=dict(color=color, width=0.8),
            marker=dict(symbol=marker_symbol, size=8) if marker_symbol is not None else None,
            fill="tozeroy" if fill_plot else None,
            fillcolor=hex_to_rgba(color, 0.2) if fill_plot else None
        )
        traces.append(trace)

    # Create the figure with all traces.
    fig = go.Figure(data=traces)

    # Update the layout: axis titles, legend position, figure size, and margins.
    fig.update_layout(
        title=title,
        width=900,  # Approximately 3:1 aspect ratio.
        height=300,
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-0.3  # Legend placed below the plot.
        ),
        margin=dict(l=50, r=50, t=30, b=100),
        xaxis_title=x_axis_title,
        yaxis_title=y_axis_title
    )

    # Update y-axis: set a fixed number of major ticks and adjust standoff.
    fig.update_yaxes(nticks=6, title_standoff=13)

    # Update x-axis: set minor ticks, grid properties, and optionally the range.
    fig.update_xaxes(
        title_standoff=8,
        tickmode='auto',
        showgrid=True,
        linewidth=1,
        linecolor='grey',
        ticks="outside",  # Draw ticks outside the axis.
        ticklen=5,        # Tick length.
        tickwidth=1,      # Tick width.
        tickcolor="grey",
        minor=dict(
            showgrid=True,
            tickmode='auto',
            nticks=x_axis_minor_ticks
        )
    )

    if x_axis_range:
        fig.update_layout(xaxis_range=x_axis_range)

    return fig

# move to gui_shared
def plot_interactive_donut_chart(group_count_df, color_map, group='Group', inner_radius=50):
    """
    Creates an interactive donut chart with an associated legend using Altair.

    The function takes a DataFrame (group_count_df) containing group labels and counts, renames its first two columns
    to "Group" and "Count", and then generates a donut chart where the donut segments represent the counts per group.
    An interactive legend is created that allows multiple selections (via Shift+click) to filter the donut chart.
    Additional layers display the aggregated total in the center of the donut and an outer ring indicating clonal data.

    Parameters:
        group_count_df (pd.DataFrame): DataFrame with at least two columns, where the first column contains group labels
                                       and the second column contains numeric counts.
        color_map (dict): Dictionary mapping each group label to a corresponding hex color code.
        group (str): Title for the group in the legend (default: "Group").
        inner_radius (int): The inner radius (in pixels) of the donut chart (default: 50).

    Returns:
        alt.Chart: An Altair Chart object representing the combined interactive donut chart and legend.
    """
    data = group_count_df.copy()
    data.rename(columns={
        data.columns[0]: 'Group',
        data.columns[1]: 'Count'
    }, inplace=True)
    col1, col2 = data.columns[0:2]
    data['order'] = range(len(data))

    if "Undefined" not in color_map:
        color_map['Undefined'] = '#FFFFFF'

    data[col1] = data[col1].fillna('Undefined')

    if all(data[col2].isnull()):
        st.warning('No valid clone data.')

    # Preserve the order of groups as they appear in the DataFrame.
    groups = data[col1].to_list()
    data[col1] = pd.Categorical(data[col1], categories=groups, ordered=True)

    # Define an interactive multi-selection bound to the legend.
    selection = alt.selection_point(
        fields=[col1],
        bind='legend'
    )

    legend_columns = math.ceil(len(data) / 8)

    # Legend Chart (separate): This chart is used solely to display the legend.
    legend_chart = (
        alt.Chart(data)
        .mark_point(
            filled=True,
            shape='square',
            strokeWidth=0,
            size=100
        )
        .encode(
            color=alt.Color(
                f"{col1}:N",
                scale=alt.Scale(
                    domain=groups,
                    range=[color_map[g] for g in groups]
                ),
                legend=alt.Legend(
                    title=f"{group} (Shift+click for multiple selections)",
                    titleLimit=1000,
                    orient='right',
                    direction='horizontal',
                    columns=legend_columns,
                    symbolType='square',
                    symbolLimit=9999,
                    offset=-150,
                    values=groups
                )
            )
        )
        # Trick: No visible plot is displayed; only the legend is shown.
        .transform_filter("false")
        .add_params(selection)
        .properties(width=150, height=200)
    )

    # Donut Chart (without its own legend)
    donut = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(total=f"sum({col2})")  # Aggregate total for percentage calculation.
        .transform_calculate(percentage=f"datum.{col2} / datum.total")
        .mark_arc(innerRadius=inner_radius, outerRadius=75)
        .encode(
            theta=alt.Theta(field=col2, type="quantitative"),
            color=alt.Color(
                f"{col1}:N",
                scale=alt.Scale(
                    domain=groups,
                    range=[color_map[g] for g in groups]
                ),
                legend=None  # Disable the built-in legend.
            ),
            order=alt.Order('order:Q'),
            tooltip=[
                alt.Tooltip(f"{col1}:N", title=group),
                alt.Tooltip(f"{col2}:Q", title="Count"),
                alt.Tooltip("percentage:Q", format=".1%", title="Percent")
            ]
        )
        .properties(width=200, height=200)
    )

    # Central text displaying the aggregated total count.
    selected_text = alt.Chart(data).transform_filter(selection) \
        .transform_aggregate(total=f"sum({col2})") \
        .mark_text(size=24, align='left', dx=3, color='black') \
        .encode(
            y=alt.Y(f"{col1}:N"),
            x=alt.X(f"total:Q", title="Sum"),
            text=alt.Text("total:Q", format='d')
        )

    # Outer background arc for the total (used as a background for the outer clonal layer).
    outer_total = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(total_all='sum(Count)')
        # Only one aggregated value is present, filling the entire circle.
        .mark_arc(innerRadius=83, outerRadius=91, color="#f3f3f3")
        .encode(
            theta=alt.Theta("total_all:Q", stack=None)
        )
    )

    # Outer arc representing "clonal" counts (everything except Undefined and Non-clonal).
    outer_clonal = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(total_all='sum(Count)')
        .transform_calculate(
            clonal_expr="(datum.Group != 'Undefined' && datum.Group != 'Non-clonal') ? datum.Count : 0")
        .transform_joinaggregate(clonal='sum(clonal_expr)')
        .transform_calculate(
            percentage="datum.total_all > 0 ? datum.clonal / datum.total_all : 0")
        .mark_arc(innerRadius=83, outerRadius=91, color="gray")
        .encode(
            theta=alt.Theta("clonal:Q", stack=None),
            tooltip=[
                alt.Tooltip("clonal:Q", title="Clonal Count"),
                alt.Tooltip("percentage:Q", format=".1%", title="Percent (Clonal)")
            ]
        )
    )

    # Combine the outer arcs into one layer.
    donut_outer = alt.layer(outer_total, outer_clonal)

    # Combine the outer donut layer, the inner donut chart, and the central text.
    donut_layer = alt.layer(donut_outer, donut, selected_text)

    # Horizontally concatenate the donut chart and the legend chart.
    combined = alt.hconcat(
        donut_layer,
        legend_chart,
        spacing=1
    ).resolve_legend(
        color='independent'
    )

    return combined


def plot_interactive_box_plot(
        df: pd.DataFrame,
        plot_type: str = "box",
        color_dict: dict = None,
        vertical_gap: float = 1,
        figure_height: int = None,
        xaxis_title: str = "GRAVY score",
        xaxis_resolution: float = 0.5,
        yaxis_resolution: float = 1,
        yaxis_title: str = "Dataset",
        data_col: str = "GRAVY",
        group_col: str = "Group",
        orientation: str = 'h'
    ) -> go.Figure:
    """
    Creates an interactive Plotly box plot or violin plot from the given DataFrame.

    Each group (as defined by the 'group_col') is assigned a unique y position, and for each group,
    the distribution of the values in the specified 'data_col' is plotted as a box plot or violin plot.
    Optionally, individual data points are overlaid as scatter points. The plot supports both horizontal
    ('h') and vertical ('v') orientations.

    Parameters:
        df (pd.DataFrame): DataFrame containing at least the columns specified by 'group_col' and 'data_col'.
        plot_type (str): Type of plot to create; either "box" or "violin". Raises ValueError for other values.
        color_dict (dict): Dictionary mapping group names to color codes (hex strings). If not provided, a default color is used.
        vertical_gap (float): Numeric gap between groups on the y-axis. Default is 1.
        figure_height (int): Overall figure height in pixels. If None, the height is computed based on the number of groups.
        xaxis_title (str): Title for the x-axis.
        xaxis_resolution (float): Resolution for the x-axis ticks (major ticks are set at twice this value; default is 0.5).
        yaxis_resolution (float): Resolution for the y-axis minor ticks (default is 1).
        yaxis_title (str): Title for the y-axis.
        data_col (str): Column name in df to be used as data values (default is "GRAVY").
        group_col (str): Column name in df to be used as group labels (default is "Group").
        orientation (str): Plot orientation, 'h' for horizontal or 'v' for vertical. Default is 'h'.

    Returns:
        go.Figure: The interactive Plotly Figure object containing the box/violin plot.
    """
    if plot_type not in ["box", "violin"]:
        raise ValueError("plot_type must be either 'box' or 'violin'.")

    fig = go.Figure()

    # Sort groups and assign each a numeric y position.
    groups = sorted(df[group_col].unique())
    y_positions = {group: i * vertical_gap for i, group in enumerate(groups)}

    for group in groups:
        group_df = df[df[group_col] == group]
        color = color_dict.get(group, "#1f77b4") if color_dict else "#1f77b4"
        fill_color = hex_to_rgba(color, 0.3)
        y_val = y_positions[group]

        x = group_df[data_col]
        y = [y_val] * len(group_df)

        if orientation == 'v':
            # Swap x and y values for vertical orientation.
            new_y = x.copy()
            x = y.copy()
            y = new_y

        if plot_type == "box":
            fig.add_trace(go.Box(
                x=x,
                y=y,
                name=group,
                orientation=orientation,
                boxpoints=False,  # Custom scatter points added below.
                line=dict(width=0.8, color="black"),
                fillcolor=fill_color,
                showlegend=False,
                hoverlabel=dict(font=dict(color=color))
            ))
        elif plot_type == "violin":
            fig.add_trace(go.Violin(
                x=x,
                y=y,
                name=group,
                orientation=orientation,
                points=False,
                line=dict(width=0.8, color="black"),
                fillcolor=fill_color,
                box=dict(visible=True, fillcolor=fill_color, line=dict(width=1.3, color="black")),
                showlegend=False,
                hoverlabel=dict(font=dict(color=color))
            ))

        # Add individual data points as scatter traces.
        x_vals = group_df[data_col].values
        y_vals = [y_val] * len(group_df)

        if orientation == 'v':
            new_y_vals = x_vals.copy()
            x_vals = y_vals.copy()
            y_vals = new_y_vals

        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            name=group,
            mode="markers",
            marker=dict(
                symbol="circle",
                size=8,
                line=dict(width=0.5, color=color),
                color="rgba(0,0,0,0)",  # Transparent fill.
                opacity=0.6
            ),
            showlegend=False
        ))

    if orientation == 'v':
        # Configure y-axis for vertical orientation.
        fig.update_yaxes(
            type="linear",
            tickmode="linear",
            gridcolor="lightgrey",
            gridwidth=0.5,
            dtick=yaxis_resolution * 5,
            minor=dict(
                tickmode="linear",
                dtick=yaxis_resolution,
                showgrid=True,
                gridcolor="#F0F0F0"
            )
        )
        # Configure x-axis for vertical orientation.
        fig.update_xaxes(
            showline=True,
            linewidth=1,
            linecolor='grey',
            tickmode="array",
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            tickcolor="grey",
            tick0=0,
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys())
        )
        fig.update_layout(
            height=300,
            legend=dict(
                orientation="h",
                x=0.5,
                xanchor="center",
                y=-0.3
            ),
            margin=dict(l=50, r=50, t=30, b=100),
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title
        )
    else:
        # Configure y-axis for horizontal orientation.
        fig.update_yaxes(
            type="linear",
            tickmode="array",
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys()),
            range=[-0.5, (len(groups) - 1) * vertical_gap + 0.5],
            title_text=yaxis_title,
            autorange="reversed"  # Invert order.
        )
        # Configure x-axis for horizontal orientation.
        fig.update_xaxes(
            title_text=xaxis_title,
            showline=True,
            linewidth=1,
            linecolor='grey',
            tickmode="linear",
            showgrid=True,
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            tickcolor="grey",
            tick0=0,
            dtick=xaxis_resolution * 2,
            minor=dict(
                tickmode="linear",
                dtick=xaxis_resolution,
                showgrid=True,
                gridcolor="#F0F0F0"
            )
        )
        fig.update_layout(
            margin=dict(l=50, r=20, t=30, b=50),
            height=figure_height if figure_height is not None else 200 + 50 * len(groups)
        )

    return fig


def plot_interactive_bargraph_from_dict(results: dict,
                                        color_dict: dict,
                                        y_title: str = "Y axis",
                                        avg: str = 'Mean',
                                        std: str = 'Std'
                                        ) -> go.Figure:
    """
    Plots an interactive Plotly bar chart from a dictionary of results.

    The input dictionary should have keys as dataset names and values as dictionaries with keys
    'n', 'mean', and 'std'. Here, 'n' is the total sample count, 'mean' is the average value (e.g., diversity index),
    and 'std' is the standard deviation. The function assigns colors to the bars using the provided color_dict.

    Parameters:
        results (dict): Dictionary with dataset names as keys and dictionaries as values containing:
                        - 'n': Total sample count.
                        - 'mean': The average value.
                        - 'std': The standard deviation.
        color_dict (dict): Dictionary mapping dataset names to hex color codes.
        y_title (str): Title for the y-axis.
        avg (str): Key to use for the average value (default is 'Mean').
        std (str): Key to use for the standard deviation (default is 'Std').

    Returns:
        go.Figure: The resulting interactive Plotly bar chart.
    """
    # Extract labels, mean values, standard deviations, and sample counts.
    x_labels, y_values, std_values, customdata = [], [], [], []
    for key, val in results.items():
        x_labels.append(key)
        y_values.append(val[avg.lower()])
        std_values.append(val[std.lower()])
        customdata.append((val["n"], val[std.lower()]))

    fig = go.Figure(go.Bar(
        x=x_labels,
        y=y_values,
        error_y=dict(
            type="data",
            array=std_values,
            thickness=0.5,  # Error bar line thickness.
            width=15,       # Error bar cap width in pixels (~50% of the bar width).
            color="black",
            visible=True
        ),
        marker=dict(
            color=[hex_to_rgba(color_dict.get(key, "#1f77b4"), 0.3) for key in x_labels],
            line=dict(color="black", width=0.5)
        ),
        customdata=customdata,
        hovertemplate=(
            "<b>%{x}</b><br>"
            + avg + ": %{y:.2f}<br>"
            + std + ": %{customdata[1]:.2f}<br>"
            "n: %{customdata[0]}<extra></extra>"
        )
    ))

    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor='grey',
        tickmode="linear",
        showgrid=True,
        ticks="outside",  # Draws the ticks outside the axis.
        ticklen=5,        # Length of the ticks.
        tickwidth=1,      # Width of the ticks.
        tickcolor="grey"
    )

    fig.update_layout(
        xaxis_title="Dataset",
        yaxis_title=y_title,
        margin=dict(l=50, r=20, t=30, b=50)
    )

    return fig

# ==============================================
# Passing of global settings from session state
# ==============================================

page_name = os.path.splitext(os.path.basename(__file__))[0]
base_user_path = st.session_state.global_constants['base_user_path']

# Get default navigation start.
if page_name not in st.session_state.current_paths:
    st.session_state.current_paths[page_name] = base_user_path

# default data folder
default_data_folder = st.session_state.current_paths['latest_data_path']

# =======================
# Page content - Sidebar
# =======================

# Sidebar
st.sidebar.header("Settings", divider="rainbow")

# Datafolder selection
folder_selectbox(page_name, base_user_path, default_data_folder)
data_folder = st.session_state.current_paths[page_name]
output_folder = data_folder

# File selection
excel_files=[]
selected_files=[]

pre_filter = st.sidebar.toggle("Only show compatible 'clustered' files", value=True)
if os.path.isdir(data_folder):
    # pre-filter for standard name
    if pre_filter:
        # display only all-sequence files
        excel_files = glob.glob(os.path.join(data_folder, "*clustered*.xlsx"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*clustered*.xls")))
    else:
        # display all excel files
        excel_files = glob.glob(os.path.join(data_folder, "*.xls"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*.xlsx")))

    excel_files = [f for f in excel_files if not os.path.basename(f).startswith("~$")]

    if excel_files:
        selected_files = st.sidebar.multiselect(
            "Select Excel files for analysis",
            options=excel_files,
            default=[f for f in excel_files if "clustered" in f],
            format_func=os.path.basename
        )
    else:
        st.sidebar.warning("No Excel files found in the folder..")
else:
    st.warning("Data folder not found..")

## DATA LOADING ##

if st.sidebar.button("Load Datasets", disabled=False if selected_files else True):
    # reset loaded data:
    st.session_state.bcr_results['clustered_bcrs'] = {}

    header_cols = [0, 1]
    index_cols = [0]
    try:
        for i, file in enumerate(selected_files):
            file_name = os.path.basename(file)
            bcr_df = pd.read_excel(file, header=header_cols, index_col=index_cols)
            st.session_state.bcr_results['clustered_bcrs'][file_name] = {'dataframe': bcr_df}

        st.sidebar.success("Excel file(s) successfully loaded!")
    except Exception as e:
        st.sidebar.error(f"Error while loading excel file(s): {e}")

st.sidebar.divider()

# =========================
# Page content - Main page
# =========================

# Instructions
st.title("Basic Repertoire Characteristics")
#st.write(f"**session state 'current_path':** `{st.session_state.current_paths[page_name]}`")
st.markdown("""
### Module Workflow:
- Use data in memory from _'Clonal Clustering'_ or select your ***clustered** files and press "Load Datasets" to load
the datasets
- Click through the tabs and select the information you want to investigate.
<style> div.stButton {text-align:center} </style>""",
            unsafe_allow_html=True)

tab1, tab2, tab3, tab4, tab5 = st.tabs(['Data',
                                        'Clonality',
                                        'V(D)J Gene Segments & Isotypes',
                                        'CDR3 Properties',
                                        'Somatic Mutations']
                                       )

# =================
# Data Information
# =================
with tab1:
    st.subheader("Data:")

    if not st.session_state.bcr_results['clustered_bcrs']:
        st.warning("No data loaded yet.")

        st.markdown("""
            - Set the path to your filtered ***b-cell-receptors** files
            - Select the files you want to investigate
            - Press "Load" to load the data
            """)
    else:
        # check if datasets already contain statistics and colors (in first set), if not calculate and save in session state
        file_names = [file for file in st.session_state.bcr_results['clustered_bcrs'].keys()]

        if not 'color' in st.session_state.bcr_results['clustered_bcrs'][file_names[0]]:
            n = len(file_names)
            if n <= 10:
                cmap = plt.get_cmap('tab10')
                colors = [mcolors.rgb2hex(cmap(x)) for x in np.linspace(0, n / 10, n)]
            else:
                cmap = plt.get_cmap('GnBu')
                colors = [mcolors.rgb2hex(cmap(x)) for x in np.linspace(0, 1, n)]

            for i, file_name in enumerate(file_names):
                bcr_df = st.session_state.bcr_results['clustered_bcrs'][file_name]['dataframe']
                rep_stats = get_repertoire_statistics(bcr_df)
                sample_infos = get_sample_information(bcr_df)
                st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics'] = rep_stats
                st.session_state.bcr_results['clustered_bcrs'][file_name]['sample_infos'] = sample_infos
                st.session_state.bcr_results['clustered_bcrs'][file_name]['color'] = colors[i]
                st.session_state.bcr_results['clustered_bcrs'][file_name].setdefault('alias', f"Dataset-{i + 1}")

            # Reset results session state if calculated for the first time
            st.session_state.page_states[page_name]['results'] = {}

            # go through chains and determine diversity and overlap measures
            for chain in ['heavy_chain', 'light_chain', 'kappa_light_chain', 'lambda_light_chain']:
                st.session_state.page_states[page_name]['results'][chain] = {}
                for mode in ['all', 'unique']:
                    st.session_state.page_states[page_name]['results'][chain][mode] = {}
                    # prepare dictionary with cdr3s (unique or not unique)
                    series_dict = {file_name: st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics'] \
                        [chain][mode]['CDR3_AA']['cdr3s'] for file_name in file_names}

                    # Determin shannon and inverse simpons index and save for each dataset
                    for div_index in ["shannon", "inverse_simpson"]:
                        results = compute_diversity_indices_for_subsampling(series_dict =  series_dict,
                                                                             index_type=div_index,
                                                                             iterations = 20,
                                                                             random_state = 42,
                                                                             seq_column = "CDR3_AA")

                        for file_name in file_names:
                            st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics'] \
                                [chain][mode]['CDR3_AA'][div_index] = results[file_name]

        # Show datasets and allow renaming of datasets
        colors = [st.session_state.bcr_results['clustered_bcrs'][file]['color'] for
                  file in st.session_state.bcr_results['clustered_bcrs'].keys()]

        # Show datasets and allow renaming of datasets
        # Baue den DataFrame jedes Mal neu aus den aktuellen Alias-Werten
        aliases = []
        for i, file_name in enumerate(file_names):
            alias = st.session_state.bcr_results['clustered_bcrs'][file_name].get('alias', f"Dataset-{i + 1}")
            aliases.append(alias)
        aliases_df = pd.DataFrame({
            "Alias": aliases,
            "File Name": file_names
        })
        aliases_df.index = range(1, len(aliases_df) + 1)

        column_config = {
            "Alias": st.column_config.TextColumn("Alias", disabled=False),
            "File Name": st.column_config.TextColumn("File Name", disabled=True)
        }

        st.caption("Loaded Datasets (rename 'alias' as required)")
        edited_aliases_df = st.data_editor(
            aliases_df,
            key="aliases_editor",
            use_container_width=True,
            column_config=column_config
        )

        # Update the session state with the changes from the Data Editor
        alias_dict = edited_aliases_df.set_index('File Name')['Alias'].to_dict()
        for file_name, alias in alias_dict.items():
            st.session_state.bcr_results['clustered_bcrs'][file_name]['alias'] = alias

        ### COLOR SECTION ###
        # Bestimme die Anzahl der Zeilen (bis zu 6 Elemente pro Zeile)
        num_elements = len(file_names)
        max_cols = 6
        num_rows = math.ceil(num_elements / max_cols)
        st.caption('Change dataset colors')
        color_dict = {}
        for row in range(num_rows):
            cols = st.columns(max_cols)
            for i, filename in enumerate(file_names[row * max_cols:(row + 1) * max_cols]):
                idx = row * max_cols + i
                alias = st.session_state.bcr_results['clustered_bcrs'][filename]['alias']
                with cols[i]:
                    current_color = st.session_state.bcr_results['clustered_bcrs'][filename]['color']
                    new_color = st.color_picker(label=alias, value=current_color, key=f"color_{filename}")
                    st.session_state.bcr_results['clustered_bcrs'][filename]['color'] = new_color
                color_dict[alias] = new_color


        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            alias = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            color = st.session_state.bcr_results['clustered_bcrs'][f]['color']

        ### DATA SUMMARY ###
        # Sample Information
        sample_info_df = pd.DataFrame(
            data=[st.session_state.bcr_results['clustered_bcrs'][file_name]['sample_infos'] for
                  file_name in file_names],
            index=[st.session_state.bcr_results['clustered_bcrs'][file_name]['alias'] for
                   file_name in file_names])

        st.caption('Sample Information')
        st.dataframe(sample_info_df, use_container_width=True)


# ==========
# CLONALITY
# ==========
with (tab2):
    st.subheader("Clonality:")
    if not st.session_state.bcr_results['clustered_bcrs']:
        st.warning("No data loaded yet.")

    else:
        for file_name in [file for file in st.session_state.bcr_results['clustered_bcrs'].keys()]:
            count_df = st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics']['clonality']['clone_sizes']
            clone_colors = st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics']['clonality']['clone_colors']
            hex_color = st.session_state.bcr_results['clustered_bcrs'][file_name]['color']
            current_alias = st.session_state.bcr_results['clustered_bcrs'][file_name]['alias']
            st.markdown(f'<b><p style="color: {hex_color};">{current_alias}:</p></b>', unsafe_allow_html=True)
            st.altair_chart(plot_interactive_donut_chart(count_df, clone_colors, group='Clone', inner_radius=40),
                                                        use_container_width=True)

            #st.write(st.session_state.bcr_results['clustered_bcrs'][file_name]['statistics']['clonality'])

# ==============
# GENE SEGMENTS
# ==============

def display_gene_segment_panel(bcr_chain, gene_segment, color_dict,
                               x_label='Gene segment', y_label='Abundance', title="Gene segments"):
    """
    Displays a gene segment panel with a grouped bar chart for BCR gene segment statistics.

    This function writes a title to the Streamlit app, retrieves gene segment statistics from the clustered BCR results
    stored in the session state, and generates a grouped bar chart using the provided color mapping. For light chains,
    if the requested gene segment is 'top_isotype', it is replaced with 'pcr_isotype'. The function returns both the
    DataFrame containing gene segment statistics and the corresponding Plotly figure.

    Parameters:
        bcr_chain (str): The BCR chain type (e.g., 'heavy_chain' or 'light_chain').
        gene_segment (str): The gene segment to display (e.g., 'top_v', 'top_isotype').
        color_dict (dict): Dictionary mapping group names to color codes (hex strings).
        x_label (str): Label for the x-axis of the chart.
        y_label (str): Label for the y-axis of the chart.
        title (str): Title of the panel.

    Returns:
        tuple: A tuple (gene_segment_df, gene_segment_figure) where gene_segment_df is a DataFrame containing
               gene segment statistics and gene_segment_figure is the corresponding grouped bar chart as a Plotly Figure.
    """
    st.write("**" + title + "**")
    # Get data for the figure.
    if bcr_chain == 'light_chain' and gene_segment == 'top_isotype':
        gene_segment = 'pcr_isotype'
    gene_segment_df = get_gene_segment_stats(
        st.session_state.bcr_results['clustered_bcrs'],
        bcr_chain,
        gene_segment,
        collapsed=data_reduction,
        nans=gene_segment_nans,
        output=gene_segment_y_format
    )

    # Create the figure.
    gene_segment_figure = plot_grouped_bar_chart(
        gene_segment_df,
        color_dict,
        x_axis_title=x_label,
        y_axis_title=y_label
    )
    # Return both the data and the figure.
    return gene_segment_df, gene_segment_figure

with tab3:
    # re-initialize session state for tab 3 to hold only displayed data
    if "displayed_gene_segment_data" in st.session_state.page_states[page_name]:
        del st.session_state.page_states[page_name]['displayed_gene_segment_data']
    st.session_state.page_states[page_name]['displayed_gene_segment_data'] = {}

    st.subheader("V(D)J Gene Segments & Constant Region Isotypes:")
    if not st.session_state.bcr_results['clustered_bcrs']:
        st.warning("No data loaded yet.")

    else:
        ### DATA REDUCTION ###
        setting_cols = st.columns(3, vertical_alignment="bottom")
        with setting_cols[0]:
            if st.toggle('Collapse clonal sequences', value = False):
                data_reduction = 'collapsed'
            else:
                data_reduction = 'not_collapsed'
        with setting_cols[1]:
            if st.toggle('Include missing data', value = False):
                gene_segment_nans="with_nans"
            else:
                gene_segment_nans="no_nans"
        with setting_cols[2]:
            if st.toggle('Show absolute counts', value=False):
                gene_segment_y_format = "counts"
                gene_segment_y_label = "Counts"
            else:
                gene_segment_y_format = "percent"
                gene_segment_y_label = "Frequency (%)"

        gene_segment_settings = {'V': {'column': 'top_v',
                                       'x_label': 'V gene segment',
                                       'title': 'V Genes'
                                       },
                                 'D': {'column': 'top_d',
                                       'x_label': 'D gene segment',
                                       'title': 'D Genes'
                                       },
                                 'J': {'column': 'top_j',
                                       'x_label': 'J gene segment',
                                       'title': 'J Genes'
                                       },
                                 'C': {'column': 'top_isotype',
                                       'x_label': 'Isotype',
                                       'title': 'Constant Region'
                                       }

                                 }

        ### HEAVY CHAIN SECTION ###
        st.markdown("""#### Heavy Chain""")
        fixed_order = ["V", "D", "J", "C"]
        selected_hc_gene_segments = st.segmented_control('Show:',
                                                         options=['V', 'D', 'J','C'],
                                                         default=['V'],
                                                         selection_mode='multi')
        selected_hc_gene_segments_ordered = [seg for seg in fixed_order if seg in selected_hc_gene_segments]

        if selected_hc_gene_segments_ordered:
            # Determine column widths based on selection...
            if len(selected_hc_gene_segments_ordered) == 1:
                widths = [1]
            elif len(selected_hc_gene_segments_ordered) == 2:
                if "V" in selected_hc_gene_segments_ordered:
                    widths = [2 if seg == "V" else 1 for seg in selected_hc_gene_segments_ordered]
                else:
                    widths = [1, 1]
            elif len(selected_hc_gene_segments_ordered) == 3:
                widths = [2, 1, 1]
            elif len(selected_hc_gene_segments_ordered) == 4:
                widths = [2, 1, 1, 1]
            else:
                widths = []

            if widths:
                cols = st.columns(widths)
                for i, hc_gene_segment in enumerate(selected_hc_gene_segments_ordered):
                    with cols[i]:
                        heavy_data, heavy_figure = display_gene_segment_panel('heavy_chain',
                                                                              gene_segment_settings[hc_gene_segment]['column'],
                                                                              color_dict,
                                                                              gene_segment_settings[hc_gene_segment]['x_label'],
                                                                              gene_segment_y_label,
                                                                              gene_segment_settings[hc_gene_segment]['title']
                                                                              )

                        st.plotly_chart(heavy_figure, key='HC'+hc_gene_segment)
                        # save active element in session state for download
                        st.session_state.page_states[page_name][
                            'displayed_gene_segment_data']["HeavyChain_"+hc_gene_segment] = {
                            'data':heavy_data,
                            'figure':heavy_figure}

        else:
            st.write("Nothing selected")

        ### LIGHT CHAIN SECTION ###
        st.markdown("""#### Light Chain""")
        lc_setting_columns = st.columns([1,4], vertical_alignment="bottom")
        with lc_setting_columns[0]:
            fixed_order = ["V", "J", "C"]
            selected_lc_gene_segments = st.segmented_control('Show:',
                                                             options=['V', 'J', 'C'],
                                                             default=['V', 'C'],
                                                             selection_mode='multi')
            selected_lc_gene_segments_ordered = [seg for seg in fixed_order if seg in selected_lc_gene_segments]

        with lc_setting_columns[1]:
            split_light_chains = st.checkbox("Split Kappa and Lambda statistics", value=False)
        if selected_lc_gene_segments_ordered:
            # Determine column widths based on selection...
            if len(selected_lc_gene_segments_ordered) == 1:
                widths = [1]
            elif len(selected_lc_gene_segments_ordered) == 2:
                if 'V' in selected_lc_gene_segments_ordered:
                    widths = [2,1]
                else:
                    widths = [1,1]
            else:
                widths = [2,1,1]

            if widths:
                if split_light_chains:
                    kappa_cols = st.columns(widths)
                    for i, lc_gene_segment in enumerate(selected_lc_gene_segments_ordered):
                        with kappa_cols[i]:
                            kappa_data, kappa_figure = display_gene_segment_panel('kappa_light_chain',
                                                                                  gene_segment_settings[lc_gene_segment]['column'],
                                                                                  color_dict,
                                                                                  gene_segment_settings[lc_gene_segment]['x_label'],
                                                                                  gene_segment_y_label,
                                                       "Kappa " + gene_segment_settings[lc_gene_segment]['title']
                                                                                  )
                            st.plotly_chart(kappa_figure, key='KAPPA' + lc_gene_segment)
                            st.session_state.page_states[page_name][
                                'displayed_gene_segment_data']["KappaChain_" + lc_gene_segment] = {
                                'data': kappa_data,
                                'figure': kappa_figure}


                    lambda_cols = st.columns(widths)
                    for i, lc_gene_segment in enumerate(selected_lc_gene_segments_ordered):
                        with lambda_cols[i]:
                            lambda_data, lambda_figure = display_gene_segment_panel('lambda_light_chain',
                                                                                    gene_segment_settings[lc_gene_segment]['column'],
                                                                                    color_dict,
                                                                                    gene_segment_settings[lc_gene_segment]['x_label'],
                                                                                    gene_segment_y_label,
                                                       "Lambda " + gene_segment_settings[lc_gene_segment]['title']
                                                                                    )
                            st.plotly_chart(lambda_figure, key='LAMBDA' + lc_gene_segment)
                            st.session_state.page_states[page_name][
                                'displayed_gene_segment_data']["LambdaChain_" + lc_gene_segment] = {
                                'data': lambda_data,
                                'figure': lambda_figure}

                else:
                    cols = st.columns(widths)
                    for i, lc_gene_segment in enumerate(selected_lc_gene_segments_ordered):
                        with cols[i]:
                            light_data, light_figure = display_gene_segment_panel('light_chain',
                                                                                  gene_segment_settings[lc_gene_segment]['column'],
                                                                                  color_dict,
                                                                                  gene_segment_settings[lc_gene_segment]['x_label'],
                                                                                  gene_segment_y_label,
                                                                                  gene_segment_settings[lc_gene_segment]['title']
                                                                                  )

                            st.plotly_chart(light_figure, key='LC' + lc_gene_segment)
                            st.session_state.page_states[page_name][
                                'displayed_gene_segment_data']["LightChain_" + lc_gene_segment] = {
                                'data': light_data,
                                'figure': light_figure}

        else:
            st.write("Nothing selected")

        # Download button for displayed data
    if "displayed_gene_segment_data" in st.session_state.page_states[page_name] and st.session_state.page_states[page_name][
                            'displayed_gene_segment_data']:

        # Erstelle ein in-memory ZIP-Archiv
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zip_file:
            for key, content in st.session_state.page_states[page_name][
                            'displayed_gene_segment_data'].items():
                displayed_data = content["data"]
                fig_to_save = content["figure"]

                # DataFrame als CSV
                csv_data = displayed_data.to_csv(index=True).encode("utf-8")
                zip_file.writestr(f"{key}_data.csv", csv_data)

                # DataFrame als Excel
                excel_buffer = io.BytesIO()
                with pd.ExcelWriter(excel_buffer, engine="xlsxwriter") as writer:
                    displayed_data.to_excel(writer, index=True, sheet_name=key)
                zip_file.writestr(f"{key}_data.xlsx", excel_buffer.getvalue())

                # Plotly-Figur als PDF exportieren (Kaleido erforderlich)
                pdf_buffer = io.BytesIO()
                fig_to_save.write_image(pdf_buffer, format="pdf")
                zip_file.writestr(f"{key}_figure.pdf", pdf_buffer.getvalue())

        # Wichtig: Den Buffer zurück auf den Anfang setzen
        zip_buffer.seek(0)

        # Download data
        st.download_button(
            label="Download Displayed Data",
            data=zip_buffer,
            file_name=f"{datetime.datetime.now().strftime("%Y%m%d-%H%M%S")}_displayed_data.zip",
            mime="application/zip"
        )
# ===============
# CDR3 Properties
# ===============

with tab4:
    def cdr3_length_statistics(bcr_chain,
                               x_label="CDR3 length (aa)",
                               y_label="Frequency (%)",
                               y_format='percent',
                               cdr3_selection='all',
                               cdr3_marker_list=None,
                               fill_curves=False,
                               title=''):
        """
        Computes CDR3 length statistics from the BCR results and creates an interactive line plot.

        This function extracts CDR3 length distribution data for a specified BCR chain from the
        clustered BCR results stored in the session state. For each file in the results, it retrieves
        the length distribution (formatted according to the provided y_format), assigns the dataset alias
        as the series name, and concatenates all series into a single DataFrame. The DataFrame is then reindexed
        to cover the full range of observed CDR3 lengths (with a one-unit padding on both ends) and missing values
        are filled with zero.

        Next, a color dictionary is generated based on the alias and color information in the session state.
        Minor tick settings are specified for the x-axis (with a default for heavy_chain). Finally, the function
        creates an interactive line plot using the plot_interactive_line_plot function, with options for filling
        the area under the curve and specifying marker symbols.

        Parameters:
            bcr_chain (str): The BCR chain type (e.g., 'heavy_chain', 'light_chain') for which to compute statistics.
            x_label (str): Label for the x-axis (default: "CDR3 length (aa)").
            y_label (str): Label for the y-axis (default: "Frequency (%)").
            y_format (str): Format for the y-axis values, e.g., 'percent' (default: 'percent').
            cdr3_selection (str): Specifies which CDR3 data to use; 'all' or 'unique' (default: 'all').
            cdr3_marker_list (list): List of marker symbols to be used in the plot (default: None).
            fill_curves (bool): If True, fills the area under the line curves (default: False).
            title (str): Title of the plot (default is an empty string).

        Returns:
            tuple: A tuple (combined_df, cdr3_length_figure) where:
                - combined_df (pd.DataFrame): DataFrame with CDR3 length statistics.
                - cdr3_length_figure (go.Figure): Interactive Plotly figure of the CDR3 length distribution.
        """
        dfs = []
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            # Extract the CDR3 length distribution for the specified chain and selection.
            df_series = \
            st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA'] \
                ['cdr3_length_distribution'][y_format]
            df_series.name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            dfs.append(df_series)

        combined_df = pd.concat(dfs, axis=1)
        combined_df = combined_df.reindex(
            range(int(combined_df.index.min() - 1), int(combined_df.index.max() + 2))
        ).fillna(0)

        # Generate a color dictionary from the session state's clustered BCR results.
        color_dict = {
            st.session_state.bcr_results['clustered_bcrs'][f]['alias']:
                st.session_state.bcr_results['clustered_bcrs'][f]['color']
            for f in st.session_state.bcr_results['clustered_bcrs'].keys()
        }
        x_minor_ticks = {'heavy_chain': 5}

        # Create the interactive line plot.
        cdr3_length_figure = plot_interactive_line_plot(
            combined_df,
            color_dict,
            x_axis_title=x_label,
            y_axis_title=y_label,
            fill_plot=fill_curves,
            marker_list=cdr3_marker_list,
            x_axis_minor_ticks=x_minor_ticks.get(bcr_chain, 2),
            title=title
        )

        return combined_df, cdr3_length_figure


    def cdr3_hydro_statistics(bcr_chain,
                              x_label="CDR3 Hydrophobicity",
                              y_label="Frequency (%)",
                              cdr3_selection='all',
                              cdr3_marker_list=None,
                              scale='kyte-doolittle',
                              title=''):
        """
        Computes CDR3 hydrophobicity statistics for a given BCR chain and creates an interactive box plot.

        The function extracts hydrophobicity data (e.g., GRAVY scores) from the clustered BCR results stored in
        the session state for the specified BCR chain and CDR3 selection. For each file in the results, it retrieves
        the hydrophobicity series corresponding to the given scale ('kyte-doolittle' or 'eisenberg'), associates the
        dataset alias and its assigned color, and compiles these series into a combined DataFrame. Then, it creates an
        interactive box plot using the plot_interactive_box_plot function.

        Parameters:
            bcr_chain (str): The BCR chain type (e.g., 'heavy_chain' or 'light_chain').
            x_label (str): Label for the x-axis (default: "CDR3 Hydrophobicity").
            y_label (str): Label for the y-axis (default: "Frequency (%)").
            cdr3_selection (str): Specifies whether to use 'all' CDR3 sequences or a 'unique' subset (default: 'all').
            cdr3_marker_list (list): List of marker symbols for the plot (default is None).
            scale (str): Hydrophobicity scale to use ('kyte-doolittle' or 'eisenberg'; default is 'kyte-doolittle').
            title (str): Title for the panel (default is an empty string).

        Returns:
            tuple: A tuple (combined_df, cdr3_hydro_figure) where:
                - combined_df (pd.DataFrame): DataFrame with CDR3 hydrophobicity statistics.
                - cdr3_hydro_figure (go.Figure): Interactive Plotly figure (box plot) displaying the GRAVY scores.
        """
        df_list = []
        color_dict = {}
        # Define x-axis resolution based on the selected scale.
        x_res = {'kyte-doolittle': 0.5, 'eisenberg': 0.2}

        # Iterate over each file in the clustered BCR results.
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            # Extract the hydrophobicity series for the given chain and CDR3 selection.
            series = \
            st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA'] \
                ['cdr3_hydrophobicity'][scale]
            alias = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            color = st.session_state.bcr_results['clustered_bcrs'][f]['color']

            # Create a temporary DataFrame with the hydrophobicity values and group label.
            temp_df = pd.DataFrame({"GRAVY": series})
            temp_df["Group"] = alias
            df_list.append(temp_df)
            color_dict[alias] = color

        # Combine all temporary DataFrames into one.
        combined_df = pd.concat(df_list, ignore_index=True)

        # Create the interactive box plot.
        cdr3_hydro_figure = plot_interactive_box_plot(
            combined_df,
            plot_type="box",
            color_dict=color_dict,
            x_axis_title="GRAVY score",
            xaxis_resolution=x_res[scale],
            yaxis_title="Dataset",
            data_col="GRAVY",
            group_col="Group"
        )

        return combined_df, cdr3_hydro_figure


    def cdr3_charge_statistics(bcr_chain,
                               x_label="Net charge at pH7.4",
                               y_label="Frequency (%)",
                               cdr3_selection='all',
                               cdr3_marker_list=None,
                               title=''):
        """
        Computes CDR3 net charge statistics for a given BCR chain and creates an interactive box plot.

        This function extracts the CDR3 net charge data from the clustered BCR results stored in the session state
        for the specified BCR chain and CDR3 selection. For each file in the results, it retrieves the net charge series,
        associates the dataset alias and its assigned color, and compiles these series into a single DataFrame.
        Then, an interactive box plot is created using the plot_interactive_box_plot function.

        Parameters:
            bcr_chain (str): The BCR chain type (e.g., 'heavy_chain' or 'light_chain').
            x_label (str): Label for the x-axis (default: "Net charge at pH7.4").
            y_label (str): Label for the y-axis (default: "Frequency (%)").
            cdr3_selection (str): Specifies whether to use 'all' CDR3 sequences or a 'unique' subset (default: 'all').
            cdr3_marker_list (list): List of marker symbols for the plot (default: None).
            title (str): Title for the panel (default is an empty string).

        Returns:
            tuple: A tuple (combined_df, cdr3_charge_figure) where:
                - combined_df (pd.DataFrame): DataFrame with CDR3 net charge statistics.
                - cdr3_charge_figure (go.Figure): Interactive Plotly figure (box plot) displaying the net charge distribution.
        """
        df_list = []
        color_dict = {}

        # Iterate over each file in the clustered BCR results.
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            # Extract the CDR3 net charge series for the specified chain and selection.
            series = \
            st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA'] \
                ['cdr3_net_charges']
            alias = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            color = st.session_state.bcr_results['clustered_bcrs'][f]['color']

            # Create a temporary DataFrame for the current dataset.
            temp_df = pd.DataFrame({"Charge": series})
            temp_df["Group"] = alias
            df_list.append(temp_df)
            color_dict[alias] = color

        # Combine all temporary DataFrames into a single DataFrame.
        combined_df = pd.concat(df_list, ignore_index=True)

        # Create the interactive box plot for CDR3 net charge.
        cdr3_charge_figure = plot_interactive_box_plot(
            combined_df,
            plot_type="box",
            color_dict=color_dict,
            x_axis_title=x_label,
            xaxis_resolution=0.5,
            yaxis_title="Dataset",
            data_col="Charge",
            group_col="Group"
        )

        return combined_df, cdr3_charge_figure


    def cdr3_diversity_indices(bcr_chain,
                               color_dict,
                               index_type='shannon',
                               y_title='Y axis',
                               cdr3_selection='all',
                               resampling=20):
        """
        Computes CDR3 diversity indices for the specified BCR chain and generates an interactive bar chart.

        This function extracts diversity index values from the clustered BCR results stored in the session state.
        For each file in the 'clustered_bcrs', it retrieves the diversity index (e.g., Shannon or Inverse Simpson)
        associated with the given BCR chain and CDR3 selection. It then creates an interactive bar chart using the
        plot_interactive_bargraph_from_dict function, where each bar corresponds to a dataset's diversity index.

        Parameters:
            bcr_chain (str): The BCR chain type (e.g., 'heavy_chain' or 'light_chain').
            color_dict (dict): A dictionary mapping dataset aliases to color codes (hex strings).
            index_type (str): The type of diversity index to use ('shannon' or 'inverse_simpson'). Default is 'shannon'.
            y_title (str): Title for the y-axis of the bar chart.
            cdr3_selection (str): Specifies whether to use all CDR3 sequences ('all') or a unique subset ('unique'). Default is 'all'.
            resampling (int): Number of resampling iterations (currently not used in the computation). Default is 20.

        Returns:
            tuple: A tuple (diversity_results, div_figure) where:
                - diversity_results (dict): Dictionary with dataset aliases as keys and their diversity index values.
                - div_figure (go.Figure): Interactive Plotly bar chart of the diversity indices.
        """
        diversity_results = {}
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            diversity_results[name] = \
            st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA'][
                index_type]

        div_figure = plot_interactive_bargraph_from_dict(diversity_results,
                                                         color_dict,
                                                         y_title=y_title)
        return diversity_results, div_figure

    st.subheader("CDR3 Properties:")
    if not st.session_state.bcr_results['clustered_bcrs']:
        st.warning("No data loaded yet.")

    else:
        ### DATA REDUCTION ###
        setting_cols = st.columns(3, vertical_alignment="bottom")
        with setting_cols[0]:

            if st.toggle('Unique CDR3s only', value=False):
                cdr3_mode = 'unique'
            else:
                cdr3_mode = 'all'

            #if st.toggle('Show marker symbols', value=False):
            #    cdr3_marker_list = ["circle-open", "square-open", "diamond-open", "cross", "x",
            #"triangle-up-open", "triangle-down-open", "triangle-left-open", "triangle-right-open"]
            #else:
            #    cdr3_marker_list = None
        with setting_cols[1]:
            cdr3_fill_curves = st.toggle('Colorize curve areas', value=False)
        with setting_cols[2]:
            if st.toggle('Show absolute counts', key='cdr3_ylabel', value=False):
                cdr3_y_format = "counts"
                cdr3_segment_y_label = "Counts"
            else:
                cdr3_y_format = "percent"
                cdr3_segment_y_label = "Frequency (%)"

        fixed_cdr3_order = ['heavy_chain', 'kappa_light_chain', 'lambda_light_chain']
        chain_dict = {'heavy_chain' : 'Heavy Chain',
                      'kappa_light_chain' : 'Kappa Chain',
                      'lambda_light_chain': 'Lambda Chain'}
        selected_chain = st.segmented_control('Show:',
                                                         options=fixed_cdr3_order,
                                                         default=fixed_cdr3_order[0],
                                                         key='cdr3_chain_to_display',
                                                         format_func= lambda x:chain_dict[x],
                                                         selection_mode='multi')
        selected_chain_ordered = [seg for seg in fixed_cdr3_order if seg in selected_chain]

        if selected_chain_ordered:
            # Determine column widths based on selection...
            if len(selected_chain_ordered) == 1:
                widths = [1]
            elif len(selected_chain_ordered) == 2:
                if 'Heavy chain' in selected_chain_ordered:
                    widths = [2 if seg == 'Heavy chain' else 1 for seg in selected_chain_ordered]
                else:
                    widths = [1, 1]
            elif len(selected_chain_ordered) == 3:
                widths = [2, 1, 1]
            else:
                widths = []

            if widths:
                with st.container():
                    st.markdown('''#### CDR3 length''')
                    cols = st.columns(widths)
                    for i, cdr3_chain in enumerate(selected_chain_ordered):
                        with cols[i]:
                            st.write(f"**{chain_dict[cdr3_chain]}**")
                            cdr3_df, cdr3_fig = cdr3_length_statistics(cdr3_chain,
                                                                       title=chain_dict[cdr3_chain],
                                                               y_label=cdr3_segment_y_label,
                                                               y_format=cdr3_y_format,
                                                               cdr3_selection=cdr3_mode,
                                                               fill_curves=cdr3_fill_curves
                                                               )
                            st.plotly_chart(cdr3_fig)
                with st.container():
                    st.markdown('''#### CDR3 Hydrophobicity''')

                    hydro_scale = st.radio("Hydrophobicity scale",
                             ['kyte-doolittle', 'eisenberg'],
                             format_func=lambda x: {'kyte-doolittle':'Kyte-Doolittle',
                                                    'eisenberg':'Eisenberg'}.get(x),
                                           horizontal=True
                             )
                    cols = st.columns(widths)
                    for i, cdr3_chain in enumerate(selected_chain_ordered):
                        with cols[i]:
                            st.write(f"**{chain_dict[cdr3_chain]}**")
                            hydro_df, hydro_fig = cdr3_hydro_statistics(cdr3_chain,
                                                                       title=chain_dict[cdr3_chain],
                                                                       y_label=cdr3_segment_y_label,
                                                                       cdr3_selection=cdr3_mode,
                                                                       scale=hydro_scale
                                                                       )
                            st.plotly_chart(hydro_fig)

                with st.container():
                    st.markdown('''#### CDR3 Charge''')
                    cols = st.columns(widths)
                    for i, cdr3_chain in enumerate(selected_chain_ordered):
                        with cols[i]:
                            st.write(f"**{chain_dict[cdr3_chain]}**")
                            charge_df, charge_fig = cdr3_charge_statistics(cdr3_chain,
                                                                           cdr3_selection=cdr3_mode)
                            st.plotly_chart(charge_fig)

                with st.container():
                    st.markdown('''#### CDR3 Diversity''')
                    index_dict = {'shannon': 'Shannon Entropy',
                     'inverse_simpson': 'Inverse Simpson Index'}
                    diversity_index = st.radio("Diversity Index",
                                           ['shannon', 'inverse_simpson'],
                                           format_func=lambda x: index_dict.get(x),
                                           horizontal=True
                                           )
                    cols = st.columns(widths)
                    for i, cdr3_chain in enumerate(selected_chain_ordered):
                        with cols[i]:
                            st.write(f"**{chain_dict[cdr3_chain]}**")
                            diversity_data, diversity_figure = cdr3_diversity_indices(cdr3_chain,color_dict,
                                                   index_type=diversity_index,
                                                   cdr3_selection=cdr3_mode,
                                                   resampling=20,
                                                    y_title=index_dict.get(diversity_index)
                                                   )
                            st.plotly_chart(diversity_figure)


# ==================
# SOMATIC MUTATIONS
# ==================
def v_ident_statistics(bcr_chain,
                       color_dict,
                       x_label="V gene germline identity (%)",
                       y_label="Frequency (%)",
                       y_format='percent',
                       collapsed='non_collapsed',
                       curve_fill=False,
                       nans='no_nans',
                       clone_marker_list=None,
                       avg='mean',
                       std='std',
                       title=''):
    """
    Computes V gene identity statistics for a given BCR chain and generates interactive plots.

    This function extracts V gene identity distribution data from the clustered BCR results stored in the
    session state for the specified BCR chain and data reduction mode (collapsed or non-collapsed). It constructs
    two interactive plots:
      1. A line plot displaying the V gene identity distribution across datasets.
      2. A violin plot summarizing the V gene identity values (referred to as 'identities') for each dataset.
    The x-axis range for the line plot is dynamically determined based on the minimum non-zero index of the distribution.

    Parameters:
        bcr_chain (str): The BCR chain type (e.g., 'heavy_chain' or 'light_chain').
        color_dict (dict): Dictionary mapping dataset aliases to color codes (hex strings).
        x_label (str): Label for the x-axis of the distribution plot (default: "V gene germline identity (%)").
        y_label (str): Label for the y-axis of the distribution plot (default: "Frequency (%)").
        y_format (str): Format key for the histogram data (default: 'percent').
        collapsed (str): Data reduction mode, either 'non_collapsed' or 'collapsed' (default: 'non_collapsed').
        curve_fill (bool): If True, fills the area under the distribution curves (default: False).
        nans (str): Specifies whether to use values including NaNs ('with_nans') or excluding them ('no_nans').
        clone_marker_list (list): List of marker symbols for the interactive line plot (default is None).
        avg (str): Key for the average value to be used in the statistics (default: 'mean').
        std (str): Key for the standard deviation (default: 'std').
        title (str): Title for the panel (default is an empty string).

    Returns:
        tuple: A tuple (combined_df, v_ident_distribution, avg_figure) where:
            - combined_df (pd.DataFrame): DataFrame containing the V gene identity distribution for each dataset.
            - v_ident_distribution (go.Figure): Interactive line plot showing the V gene identity distribution.
            - avg_figure (go.Figure): Interactive violin plot summarizing the V gene identity values per dataset.
    """
    dist_dfs = []
    flat_dfs = []
    for f in st.session_state.bcr_results['clustered_bcrs'].keys():
        # Extract the distribution of V gene identities.
        df_series = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][collapsed]['v_identity'][nans][y_format]
        name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
        df_series.name = name
        dist_dfs.append(df_series)

        # Extract the raw identity values.
        series = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][collapsed]['v_identity'][nans]['identities']
        temp_df = pd.DataFrame({"Identities": series})
        temp_df["Group"] = name
        flat_dfs.append(temp_df)

    flat_df = pd.concat(flat_dfs, ignore_index=True)
    combined_df = pd.concat(dist_dfs, axis=1)

    # Determine the x-axis range based on the minimal index where data exists.
    max_index = combined_df.index[(combined_df != 0).any(axis=1)].min()
    if max_index > 5:
        max_index -= 5
    else:
        max_index = 0
    x_range = [max_index, 100]

    # Create an interactive line plot for the V gene identity distribution.
    v_ident_distribution = plot_interactive_line_plot(
        combined_df,
        color_dict,
        x_axis_title=x_label,
        y_axis_title=y_label,
        fill_plot=curve_fill,
        marker_list=clone_marker_list,
        x_axis_minor_ticks=5,
        x_axis_range=x_range,
        title=''
    )

    # Create an interactive violin plot summarizing the V gene identity values.
    avg_figure = plot_interactive_box_plot(
        flat_df,
        plot_type="violin",
        color_dict=color_dict,
        figure_height=None,
        xaxis_title="Datasets",
        xaxis_resolution=0.5,
        yaxis_title="%",
        data_col="Identities",
        group_col="Group",
        orientation='v'
    )

    return combined_df, v_ident_distribution, avg_figure

with tab5:
    st.subheader("Somatic Mutations:")
    if not st.session_state.bcr_results['clustered_bcrs']:
        st.warning("No data loaded yet.")

    else:
        ### DATA REDUCTION ###
        setting_cols = st.columns(3, vertical_alignment="bottom")
        with setting_cols[0]:
            if st.toggle('Collapse clonal sequences', value=False, key='v_ident_toggle_collapse'):
                data_reduction = 'collapsed'
            else:
                data_reduction = 'not_collapsed'
        with setting_cols[1]:
            v_ident_fill_curves = st.toggle('Colorize curve areas', value=False, key='v_ident_curve_fill')
        with setting_cols[2]:
            if st.toggle('Show absolute counts', value=False, key='v_ident_toggle_counts'):
                v_ident_y_format = "counts"
                v_ident_y_label = "Counts"
            else:
                v_ident_y_format = "percent"
                v_ident_y_label = "Frequency (%)"

        fixed_vident_order = ['heavy_chain', 'kappa_light_chain', 'lambda_light_chain']
        chain_dict = {'heavy_chain': 'Heavy Chain',
                      'kappa_light_chain': 'Kappa Chain',
                      'lambda_light_chain': 'Lambda Chain'}
        selected_chain = st.segmented_control('Show:',
                                              options=fixed_vident_order,
                                              default=fixed_vident_order[0],
                                              key='vident_chain_to_display',
                                              format_func=lambda x: chain_dict[x],
                                              selection_mode='multi')

        selected_chain_ordered = [seg for seg in fixed_cdr3_order if seg in selected_chain]

        for i, v_ident_chain in enumerate(selected_chain_ordered):
            st.write(f"**{chain_dict[v_ident_chain]}**")
            mut_dist_data, mut_dist_fig, mut_avg_fig = v_ident_statistics(v_ident_chain,
                                                                      color_dict,
                                                                      x_label="V gene germline identity (%)",
                                                                      y_label=v_ident_y_label,
                                                                      y_format=v_ident_y_format,
                                                                      collapsed=data_reduction,
                                                                      curve_fill=v_ident_fill_curves,
                                                                      clone_marker_list=None)
            cols = st.columns([2,1], vertical_alignment="top")
            with cols[0]:
                st.plotly_chart(mut_dist_fig, use_container_width=True)

            with cols[1]:
                st.plotly_chart(mut_avg_fig, use_container_width=True)
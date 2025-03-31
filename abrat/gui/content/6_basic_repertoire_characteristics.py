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
                                                   create_heatmap_with_annotations)

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
    Reduces the DataFrame by returning only the unique TOP gene segment values per clone.

    If a clone column is provided, groups the DataFrame by that clone column and
    takes the unique values of the specified top_gene_col within each group.
    If no clone column is provided, it simply drops duplicate values in the top_gene_col.

    :param mode: 'unique_values' to reduce data in target_column to unique values or 'mean' to calculate a mean value
    :param df: pandas DataFrame containing gene segment data.
    :param target_column: Tuple representing the target column (e.g., ('HEAVY_CHAIN', 'TOP_V')).
    :param clone_col: Tuple representing the clone identifier column (e.g., (SAMPLE_INFORMATION, hclc_cluster_subcol)).
                      If provided, collapse is done within each clone.
    :return: A reduced DataFrame containing only the unique TOP gene segments.
    """
    if clone_col is not None and clone_col in df.columns:
        if mode == 'mean':
            # Group by the clone identifier and take unique gene segment values per group
            collapsed = df.groupby(clone_col)[[target_column]].apply(lambda x: x.mean())
        else: # for everything else make 'unique_values'
            # Group by the clone identifier and take unique gene segment values per group
            collapsed = df.groupby(clone_col)[[target_column]].apply(lambda x: x.drop_duplicates())

        # collapsed is a Series with a MultiIndex (clone, row index) – you may want to reset the index
        collapsed_df = collapsed.reset_index().drop(columns=[clone_col])
    else:
        # No clone column provided, simply drop duplicates in the top_gene column
        collapsed_df = df.drop_duplicates(subset=[target_column])

    return collapsed_df

def get_repertoire_statistics(df):
    """
    Computes basic repertoire statistics from the given DataFrame.

    Returns a dictionary with various statistics that can later be displayed
    on the dashboard.

    :param df: pandas DataFrame containing the repertoire data.
    :return: Dictionary with calculated statistics.
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
        'undefined':sum(df[(sample_col, clone_subcol)]=='Undefined')+sum(df[(sample_col, clone_subcol)].isnull()),
        'non_clonal':sum(df[(sample_col, clone_subcol)]=='Non-clonal'),
        'clonal':len(clone_df),
        'number_of_clones':len(clone_df[(sample_col, clone_subcol)].unique()),
        'mean_clone_size':clone_df.drop_duplicates((sample_col, clone_subcol))[
            (sample_col, cluster_size_subcol)].mean(),
        'median_clone_size':clone_df.drop_duplicates((sample_col, clone_subcol))[
            (sample_col, cluster_size_subcol)].median(),
        'clone_sizes':clone_sizes,
        'clone_colors':color_dict
    }

    # ------------------------------------------
    # V(D)J gene segments, Isotype and V gene identities
    # ------------------------------------------
    chain_genes = {
        'heavy_chain': ['TOP_V', 'TOP_D', 'TOP_J', 'TOP_ISOTYPE'],
        'light_chain': ['TOP_V', 'TOP_J', 'PCR_ISOTYPE'],
        'kappa_light_chain': ['TOP_V', 'TOP_J', 'TOP_ISOTYPE'],
        'lambda_light_chain': ['TOP_V', 'TOP_J', 'TOP_ISOTYPE']
    }

    # Composite cluster as clone identifier
    clone_identifier = (sample_col, hclc_cluster_subcol)

    for mode in ['not_collapsed', 'collapsed']:
        for chain, gene_list in chain_genes.items():
            # Initialize nested dictionary for the chain and mode
            stats.setdefault(chain, {})[mode] = {}

            # Set the chain key and determine working DataFrame based on chain type
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
                    working_df = df.copy()  # for chain 'light_chain'
            else:
                chain_key = chain.upper()
                working_df = df.copy()

            for gene in gene_list:
                v_gene_col = (chain_key, gene)
                if mode == 'collapsed':
                    # If a representative column is not available, collapse based on unique gene segment values
                    working_df_mode = collapse_clone_data(working_df, target_column=v_gene_col, clone_col=clone_identifier).copy()
                else:
                    working_df_mode = working_df

                # Calculate counts and percentages including NaNs
                with_nans_count = working_df_mode[v_gene_col].fillna('N.D.').value_counts()
                with_nans_percent = 100 * with_nans_count / with_nans_count.sum()

                # Calculate counts and percentages excluding NaNs
                no_nans_count = working_df_mode[v_gene_col].replace('N.D.',np.nan).value_counts() # value_counts() does not autommatically include NaNs
                no_nans_percent = 100 * no_nans_count / no_nans_count.sum()

                stats[chain][mode][gene.lower()] = {
                    'with_nans': {'counts': with_nans_count, 'percent': with_nans_percent},
                    'no_nans': {'counts': no_nans_count, 'percent': no_nans_percent}
                }

            # V gene identity
            v_ident_col = (chain_key, v_ident_subcol)
            if mode == 'collapsed':
                working_df_vident = collapse_clone_data(working_df, target_column=v_ident_col,clone_col=clone_identifier).copy()
            else:
                working_df_vident = working_df

            # V ident => nans are replaced by 0
            v_with_nans_count = compute_histogram(working_df_vident[v_ident_col].fillna(0), bin_size=1, min_val=-0,
                                                  max_val=101)
            v_with_nans_percent = 100 * v_with_nans_count / v_with_nans_count.sum()

            v_no_nans_count = compute_histogram(working_df_vident[v_ident_col], bin_size=1, min_val=0, max_val=101)
            v_no_nans_percent = 100 * v_no_nans_count / v_no_nans_count.sum()

            stats[chain][mode]['v_identity'] = {
                'with_nans': {'counts': v_with_nans_count,
                              'percent': v_with_nans_percent,
                              },
                'no_nans': {'identities': working_df_vident[v_ident_col].dropna(),
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
            # Initialize nested dictionary for the chain and mode
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
                    working_df = df.copy()  # for chain 'light_chain'
            else:
                chain_key = chain.upper()
                working_df = df.copy()

            cdr3_aa_col =  (chain_key, cdr3_aa_subcol)
            cdr3_aa_len_col = (chain_key, cdr3_aa_len_subcol)

            if mode == 'unique':
                working_df_mode = working_df.drop_duplicates(cdr3_aa_col)
            else:
                working_df_mode = working_df

            cdr3s = working_df_mode[cdr3_aa_col].dropna()
            cdr3_length_counts = working_df_mode[cdr3_aa_len_col].dropna().value_counts()
            cdr3_length_percent = 100 * cdr3_length_counts / cdr3_length_counts.sum()
            cdr3_gravy_eisenberg = compute_gravy_scores(cdr3s, 'eisenberg')
            cdr3_gravy_kyte= compute_gravy_scores(cdr3s, 'kyte-doolittle')

            cdr3_net_charges = cdr3s.apply(lambda seq: peptides.Peptide(seq).charge(pH=7.4))

            stats[chain][mode]['CDR3_AA'] = {
                'cdr3s': cdr3s,
                'cdr3_length_distribution': {'counts': cdr3_length_counts, 'percent': cdr3_length_percent},
                'cdr3_hydrophobicity': {'kyte-doolittle':cdr3_gravy_kyte,
                                        'eisenberg': cdr3_gravy_eisenberg
                                        },
                'cdr3_net_charges':cdr3_net_charges
                }

    return stats

def get_gene_segment_stats(data_dict, chain, gene_segment, collapsed = 'not_collapsed', nans='no_nans', output='percent'):
    """Generates a dataframe with statistics taken from data_dict.
    Data_dict should be the dictionary of filenames that carry the respective stats.
    """
    gene_segment_data = {}
    for dataset in data_dict.keys():
        dataset_alias = data_dict[dataset]['alias']
        gene_segment_data[dataset_alias] = data_dict[dataset]['statistics'][chain][collapsed][gene_segment][nans][output]

    gene_segment_df = pd.DataFrame(data=gene_segment_data).fillna(0).sort_index()

    return gene_segment_df

def update_aliases():
    # Dieser Callback wird direkt ausgeführt, wenn sich der Data Editor ändert.
    alias_dict = st.session_state.edited_aliases_df.set_index('File Name')['Alias'].to_dict()
    for file_name, alias in alias_dict.items():
        st.session_state.bcr_results['clustered_bcrs'][file_name]['alias'] = alias

### PLOT FUNCTIONS
def plot_grouped_bar_chart(df, color_dict, x_axis_title="X-Axis", y_axis_title="Y-Axis"):
    """
    Plots a grouped bar chart using Plotly.

    The DataFrame's index is used as x-values and each column is plotted as a separate dataset,
    using the column names as legend labels. The figure will be approximately three times as wide as high.
    The legend is placed horizontally centered below the graph.

    :param df: pandas DataFrame where index are x values and columns are y values.
    :return: A Plotly Figure object.
    """
    traces = []
    # Create a bar trace for each column in the DataFrame
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

    # Create figure with the traces
    fig = go.Figure(data=traces)

    # Update layout: use 'group' mode, set size and legend position (legend horizontal below the chart)
    fig.update_layout(
        barmode='group',
        width=900,  # Adjust width for a 3:1 aspect ratio
        height=300,
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-1.1  # Legende unterhalb des Plots
        ),
        margin=dict(l=50, r=50, t=30, b=100),
        xaxis_title=x_axis_title,
        yaxis_title=y_axis_title
        # extra margin for legend if needed
    )

    fig.update_yaxes(nticks=6)

    return fig


def hex_to_rgba(hex_color, alpha=0.5):
    """
    Converts a hex color string to an rgba string.

    Parameters:
        hex_color (str): The color in hex format (e.g. "#FF5733").
        alpha (float): The alpha (opacity) value.

    Returns:
        str: The color in rgba format.
    """
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return f'rgba({r}, {g}, {b}, {alpha})'


def plot_interactive_line_plot(df, color_dict, x_axis_title="X-Axis", y_axis_title="Y-Axis",
                               fill_plot=False,
                               marker_list=None,
                               x_axis_minor_ticks = 5,
                               x_axis_range=None,
                               title=''
                               ):
    """
    Erstellt einen interaktiven Linienplot mit Plotly.

    Jede Spalte des DataFrames wird als eigener Datensatz (Linie) dargestellt. Die Farben werden
    über color_dict zugewiesen. Optional kann per fill_plot=True die Fläche unter der Linie mit 50%
    Opazität gefüllt werden. Über marker_list können Marker-Symbole spezifiziert werden, wobei bei
    zu wenigen Symbolen die Liste zyklisch wiederverwendet wird.

    :param df: pandas DataFrame, wobei der Index die x-Werte und die Spalten die y-Werte darstellen.
    :param color_dict: Dictionary, das Spaltennamen auf Farbwerte (Hex-Code) abbildet.
    :param x_axis_title: Titel der x-Achse.
    :param y_axis_title: Titel der y-Achse.
    :param fill_plot: Bool, ob die Fläche unter der Linie gefüllt werden soll (default: False).
    :param marker_list: Liste von Marker-Symbolen (z.B. ['circle', 'square', 'diamond']). Default ist None.
    :return: Ein Plotly Figure Objekt.
    """
    traces = []

    for i, col in enumerate(df.columns):
        color = color_dict.get(col, "#000000")
        line_color = hex_to_rgba(color, 0.8)
        # Bestimme den Marker, falls eine Liste übergeben wurde
        marker_symbol = marker_list[i % len(marker_list)] if marker_list else None

        # Bestimme den Modus: "lines" oder "lines+markers"
        mode = "lines+markers" if marker_symbol is not None else "lines"

        # Marker werden als nicht gefüllte Symbole mit Linienumriss definiert
        marker_props = dict(
            symbol=marker_symbol,
            color="rgba(0,0,0,0)",  # transparente Füllung
            line=dict(color=line_color, width=0.5)
        )

        # Erstelle das Scatter-Objekt
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

    # Erstelle die Figur
    fig = go.Figure(data=traces)

    # Aktualisiere das Layout: Achsentitel, Legendenposition, Größe etc.
    fig.update_layout(
        title=title,
        width=900,  # ca. 3:1 Verhältnis
        height=300,
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-0.3  # Legende unterhalb des Plots
        ),
        margin=dict(l=50, r=50, t=30, b=100),
        xaxis_title=x_axis_title,
        yaxis_title=y_axis_title
    )

    # Versuche, mehr Major-Ticks auf der y-Achse zu zeigen
    fig.update_yaxes(nticks=6,
                     title_standoff=13
                     )

    # Minor-Ticks auf der x-Achse
    fig.update_xaxes(
        title_standoff=8,
        tickmode='auto',
        showgrid=True,
        linewidth=1,
        linecolor='grey',
        ticks="outside",  # Zeichnet die Ticks außerhalb der Achse
        ticklen=5,  # Länge der Ticks
        tickwidth=1,  # Breite der Ticks
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
def plot_interactive_donut_chart(group_count_df,
                                 color_map,
                                 group='Group',
                                 inner_radius=50):
    data = group_count_df.copy()
    data.rename(columns={
        data.columns[0]: 'Group',
        data.columns[1]: 'Count'
    }, inplace=True)
    col1, col2 = data.columns[0:2]
    data['order'] = range(len(data))

    if not "Undefined" in color_map:
        color_map['Undefined'] = '#FFFFFF'

    data[col1] = data[col1].fillna('Undefined')

    if all(data[col2].isnull()):
        st.warning('No valid clone data.')

    # set the order of groups as in dataframe
    groups = data[col1].to_list()

    data[col1] = pd.Categorical(data[col1], categories=groups, ordered=True)

    # set selection
    selection = alt.selection_point(
        fields=[col1],
        bind = 'legend'
        #empty=True,
        #toggle=True,
        #clear=False
    )

    legend_columns = math.ceil(len(data) / 8)

    # 2) Legend-Chart (separat)
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
        # Trick: kein sichtbarer Plot, wir wollen nur die Legende
        .transform_filter("false")
        .add_params(selection)
        .properties(width=150, height=200)
    )

    # 3) Donut-Chart (ohne eigene Legende)
    donut = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(total=f"sum({col2})") # Summe für die prozentuale Berechnung
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
                legend=None  # <--- Keine integrierte Legende!
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

    # 4) Zentraler Text mit aktualisierter Summe
    selected_text = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_aggregate(total=f"sum({col2})")
        .mark_text(size=24, align="center", baseline="middle", color="black")
        .encode(
            text=alt.Text("total:Q", format="d")
        )
    )

    # Background for outer clonal layer
    outer_total = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(total_all='sum(Count)')
        # Hier haben wir nur 1 Zeile mit total_all
        .mark_arc(innerRadius=83, outerRadius=91, color="#f3f3f3")
        .encode(
            # Ein einzelner Wert => füllt den ganzen Kreis
            theta=alt.Theta("total_all:Q", stack=None)
        )
    )

    # Outer Arc = "Clonal": everything except Undefined & Non-clonal
    outer_clonal = (
        alt.Chart(data)
        .transform_filter(selection)
        .transform_joinaggregate(
            total_all='sum(Count)')
        .transform_calculate(
            clonal_expr="(datum.Group != 'Undefined' && datum.Group != 'Non-clonal') ? datum.Count : 0")
        .transform_joinaggregate(
            clonal='sum(clonal_expr)')
        .transform_calculate(
            percentage="datum.total_all > 0 ? datum.clonal / datum.total_all : 0")
        .mark_arc(innerRadius=83, outerRadius=91, color="gray")
        .encode(
            theta=alt.Theta("clonal:Q", stack=None),
            tooltip=[
                alt.Tooltip("clonal:Q", title="Clonal Count"),
#                alt.Tooltip("total_all:Q", title="Total Count"),
                alt.Tooltip("percentage:Q", format=".1%", title="Percent (Clonal)")
            ]
        )
    )

    # Outer donut layer
    donut_outer = alt.layer(
        outer_total,
        outer_clonal
    )

    donut_layer = alt.layer(donut_outer, donut, selected_text)

    # Combine chart and legend next to each other
    combined = alt.hconcat(
        donut_layer,
        legend_chart,
        spacing=1
    ).resolve_legend(
        color='independent'
    )

    # 6) In Streamlit anzeigen
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
        Creates an interactive Plotly boxplot or violin plot of GRAVY scores.

        The plot displays each group on a numeric y-axis so that the plots can be placed
        closer together. The GRAVY score is on the x-axis. For each group:
          - The summary (box or violin) is drawn with a black border (line width 0.8).
          - The fill color (with 0.5 alpha) is taken from the provided color_dict.
          - The median is highlighted by an inner box (line width 1.3).
          - Individual data points are overlaid as circles with line width 0.5, colored according 
            to the color dictionary, with an opacity of 0.5.

        Parameters:
            df (pd.DataFrame): DataFrame containing at least the columns "Group" and "GRAVY".
            plot_type (str): Either "box" or "violin", to select the plot type.
            color_dict (dict): Dictionary mapping group names to colors (hex strings).
            vertical_gap (float): The numeric gap between groups on the y-axis. Default is 0.3.
            figure_height (int): Overall figure height. If None, it's calculated based on the number of groups.

        Returns:
            go.Figure: The interactive Plotly figure.
        """
    if plot_type not in ["box", "violin"]:
        raise ValueError("plot_type must be either 'box' or 'violin'.")

    fig = go.Figure()

    # Sort groups and assign each a numeric y position
    groups = sorted(df[group_col].unique())
    y_positions = {group: i * vertical_gap for i, group in enumerate(groups)}

    for group in groups:
        group_df = df[df[group_col] == group]
        color = color_dict.get(group, "#1f77b4") if color_dict else "#1f77b4"
        fill_color = hex_to_rgba(color, 0.3)
        y_val = y_positions[group]

        x = group_df[data_col]
        y = [y_val] * len(group_df)

        if orientation=='v':
            # change x and y values
            new_y = x.copy()
            x = y.copy()
            y =  new_y

        if plot_type == "box":
            fig.add_trace(go.Box(
                x=x,
                y=y,
                name=group,
                orientation=orientation,
                boxpoints=False,  # We'll add custom scatter points below.
                line=dict(width=0.8, color="black"),
                fillcolor=fill_color,
                showlegend=False,
                hoverlabel=dict(
                    font=dict(color=color)
                    )
                )
            )

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
                hoverlabel=dict(
                    font=dict(color=color),
                )
            ))

        # Add individual data points as scatter traces.
        x_vals = group_df[data_col].values
        y_vals = [y_val] * len(group_df)

        if orientation=='v':
            # change x and y values
            new_y_vals = x_vals.copy()
            x_vals = y_vals.copy()
            y_vals =  new_y_vals

        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            name=group,
            mode="markers",
            marker=dict(
                symbol="circle",
                size=8,
                line=dict(width=0.5, color=color),
                color="rgba(0,0,0,0)",
                opacity=0.6
            ),
            showlegend=False
        ))

    if orientation == 'v':
        # Configure the y-axis
        fig.update_yaxes(
            type="linear",
            tickmode="linear",
            gridcolor="lightgrey",
            gridwidth= 0.5,
            dtick=yaxis_resolution * 5,  # Abstand der Major Ticks (z.B. alle 5 Einheiten)
            minor=dict(  # Einstellungen für Minor Ticks
                tickmode="linear",
                dtick=yaxis_resolution,
                showgrid=True,
                gridcolor="#F0F0F0"
            )
        )

        # Configure the x-axis
        fig.update_xaxes(
            showline=True,
            linewidth=1,
            linecolor='grey',
            tickmode="array",
            ticks="outside",  # Zeichnet die Ticks außerhalb der Achse
            ticklen=5,  # Länge der Ticks
            tickwidth=1,  # Breite der Ticks
            tickcolor="grey",
            tick0=0,  # Startwert für die Major Ticks
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys()),
        )

        fig.update_layout(
            height=300,
            legend=dict(
                orientation="h",
                x=0.5,
                xanchor="center",
                y=-0.3  # Legende unterhalb des Plots
            ),
            margin=dict(l=50, r=50, t=30, b=100),
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title
        )

    else:
        # Configure the y-axis
        fig.update_yaxes(
            type="linear",
            tickmode="array",
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys()),
            range=[-0.5, (len(groups) - 1) * vertical_gap + 0.5],
            title_text=yaxis_title,
            autorange="reversed" # to invert the order
        )

        # Configure the x-axis
        fig.update_xaxes(
            title_text=xaxis_title,
            showline=True,
            linewidth=1,
            linecolor='grey',
            tickmode="linear",
            showgrid=True,
            ticks="outside",  # Zeichnet die Ticks außerhalb der Achse
            ticklen=5,  # Länge der Ticks
            tickwidth=1,  # Breite der Ticks
            tickcolor="grey",
            tick0=0,  # Startwert für die Major Ticks
            dtick=xaxis_resolution*2,  # Abstand der Major Ticks (z.B. alle 5 Einheiten)
            minor=dict(  # Einstellungen für Minor Ticks
                tickmode="linear",
                dtick=xaxis_resolution,
                showgrid=True,
                gridcolor="#F0F0F0"
            )
        )

        # Configure layout
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

    Parameters:
        results (dict): Dictionary with keys as dataset names and values as dicts with keys
                        'n', 'mean' and 'std'. 'n' is the total sample count, 'mean' the average
                        diversity index and 'std' its standard deviation.
        color_dict (dict): Dictionary mapping dataset names to hex colors.

    Returns:
        go.Figure: The resulting interactive bar chart.
    """

    # Extrahiere Labels, Mittelwerte, Standardabweichungen und sample counts.
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
            thickness=0.5,  # Error bar line thickness
            width=15,  # Error bar cap width in pixels (~50% of the bar width)
            color="black",
            visible = True
        ),
        marker=dict(
            color=[hex_to_rgba(color_dict.get(key, "#1f77b4"), 0.3) for key in x_labels],
            line=dict(color="black", width=0.5)
        ),
        customdata=customdata,
        hovertemplate=(
            "<b>%{x}</b><br>"
            +avg+": %{y:.2f}<br>"
            +std+": %{customdata[1]:.2f}<br>"
            "n: %{customdata[0]}<extra></extra>")
    ))

    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor='grey',
        tickmode="linear",
        showgrid=True,
        ticks="outside",  # Zeichnet die Ticks außerhalb der Achse
        ticklen=5,  # Länge der Ticks
        tickwidth=1,  # Breite der Ticks
        tickcolor="grey",
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
    st.write("**"+title+"**")
    # get data for figure
    if bcr_chain == 'light_chain' and gene_segment == 'top_isotype':
        gene_segment = 'pcr_isotype'
    gene_segment_df = get_gene_segment_stats(st.session_state.bcr_results['clustered_bcrs'],
                                             bcr_chain,
                                             gene_segment,
                                             collapsed=data_reduction,
                                             nans=gene_segment_nans,
                                             output=gene_segment_y_format)

    # make figure
    gene_segment_figure = plot_grouped_bar_chart(gene_segment_df,
                                                                 color_dict,
                                                                 x_axis_title=x_label,
                                                                 y_axis_title=y_label)
    #return data and figure
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
                               fill_curves = False,
                               title=''):
        '''returns dataframe and figure'''
        # get data for figure

        dfs = []
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            df = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA']\
            ['cdr3_length_distribution'][y_format]

            df.name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            dfs.append(df)

        combined_df = pd.concat(dfs, axis=1)
        combined_df = combined_df.reindex(range(int(combined_df.index.min()-1),
                                                int(combined_df.index.max()+2))).fillna(0)

        # get colors
        color_dict = {st.session_state.bcr_results['clustered_bcrs'][f]['alias']:
                          st.session_state.bcr_results['clustered_bcrs'][f]['color']
                      for f in st.session_state.bcr_results['clustered_bcrs'].keys()
                      }
        x_minor_ticks = {'heavy_chain':5}

        # make figure
        cdr3_length_figure = plot_interactive_line_plot(combined_df,
                                                        color_dict,
                                                        x_axis_title=x_label,
                                                        y_axis_title=y_label,
                                                        fill_plot=fill_curves,
                                                        marker_list=cdr3_marker_list,
                                                        x_axis_minor_ticks=x_minor_ticks.get(bcr_chain, 2),
                                                        title=''
        )

        # return data and figure
        return combined_df, cdr3_length_figure

    def cdr3_hydro_statistics(bcr_chain,
                               x_label="CDR3 Hydrophobicity",
                               y_label="Frequency (%)",
                               cdr3_selection='all',
                               cdr3_marker_list=None,
                               scale='kyte-doolittle',
                               title=''):
        '''returns dataframe and figure'''
        # get data and color for figure

        df_list = []
        color_dict = {}
        x_res = {'kyte-doolittle':0.5,
                 'eisenberg':0.2}
        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            series = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA']\
            ['cdr3_hydrophobicity'][scale]
            alias = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            color = st.session_state.bcr_results['clustered_bcrs'][f]['color']

            temp_df = pd.DataFrame({"GRAVY": series})
            temp_df["Group"] = alias
            df_list.append(temp_df)
            color_dict[alias]=color

        combined_df = pd.concat(df_list, ignore_index=True)

        # make figure
        cdr3_hydro_figure = plot_interactive_box_plot(combined_df,
                                                      plot_type= "box",
                                                      color_dict = color_dict,
                                                      xaxis_title = "GRAVY score",
                                                      xaxis_resolution = x_res[scale],
                                                      yaxis_title="Dataset",
                                                      data_col="GRAVY",
                                                      group_col="Group")

        # return data and figure
        return combined_df, cdr3_hydro_figure

    def cdr3_charge_statistics(bcr_chain,
                               x_label="CDR3 Hydrophobicity",
                               y_label="Frequency (%)",
                               cdr3_selection='all',
                               cdr3_marker_list=None,
                               title=''):
        '''returns dataframe and figure'''
        # get data and color for figure

        df_list = []
        color_dict = {}

        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            series = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][cdr3_selection]['CDR3_AA']\
            ['cdr3_net_charges']
            alias = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            color = st.session_state.bcr_results['clustered_bcrs'][f]['color']

            temp_df = pd.DataFrame({"Charge": series})
            temp_df["Group"] = alias
            df_list.append(temp_df)
            color_dict[alias]=color

        combined_df = pd.concat(df_list, ignore_index=True)

        # make figure
        cdr3_charge_figure = plot_interactive_box_plot(combined_df,
                                                       plot_type= "box",
                                                       color_dict = color_dict,
                                                       xaxis_title = "Net charge at pH7.4",
                                                       xaxis_resolution = 0.5,
                                                       yaxis_title="Dataset",
                                                       data_col="Charge",
                                                       group_col="Group")

        # return data and figure
        return combined_df, cdr3_charge_figure

    def cdr3_diversity_indices(bcr_chain,
                               color_dict,
                               index_type = 'shannon',
                               y_title = 'Y axis',
                               cdr3_selection='all',
                               resampling=20
                               ):

        diversity_results = {}

        for f in st.session_state.bcr_results['clustered_bcrs'].keys():
            name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
            diversity_results[name] = st.session_state.bcr_results['clustered_bcrs'][f]['statistics']\
                [bcr_chain][cdr3_selection]['CDR3_AA'][index_type]

        div_figure = plot_interactive_bargraph_from_dict(diversity_results,
                                                         color_dict,
                                                         y_title = y_title
                                                         )

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
                       curve_fill = False,
                       nans='no_nans',
                       clone_marker_list=None,
                       avg='mean',
                       std='std',
                       title=''):
    '''returns dataframe and figure'''
    # get data for figure

    dist_dfs = []
    flat_dfs = []
    for f in st.session_state.bcr_results['clustered_bcrs'].keys():
        df = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][collapsed]['v_identity']\
        [nans][y_format]
        name = st.session_state.bcr_results['clustered_bcrs'][f]['alias']
        df.name = name
        dist_dfs.append(df)

        series = st.session_state.bcr_results['clustered_bcrs'][f]['statistics'][bcr_chain][collapsed]['v_identity']\
        [nans]['identities']
        temp_df = pd.DataFrame({"Identities": series})
        temp_df["Group"] = name
        flat_dfs.append(temp_df)


    flat_df = pd.concat(flat_dfs, ignore_index=True)

    combined_df = pd.concat(dist_dfs, axis=1)

    max_index = combined_df.index[(combined_df != 0).any(axis=1)].min()
    if max_index > 5:
        max_index -= 5
    else:
        max_index = 0
    x_range = [max_index, 100]
    # make figure
    v_ident_distribution = plot_interactive_line_plot(combined_df,
                                                color_dict,
                                                x_axis_title=x_label,
                                                y_axis_title=y_label,
                                                fill_plot = curve_fill,
                                                marker_list=clone_marker_list,
                                                x_axis_minor_ticks=5,
                                                x_axis_range=x_range,
                                                title=''
    )

    avg_figure = plot_interactive_box_plot(flat_df,
                                            plot_type= "violin",
                                            color_dict= color_dict,
                                            figure_height = None,
                                            xaxis_title = "Datasets",
                                            xaxis_resolution = 0.5,
                                            yaxis_title = "%",
                                            data_col = "Identities",
                                            group_col = "Group",
                                            orientation = 'v'
    )

    # return data and figure
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
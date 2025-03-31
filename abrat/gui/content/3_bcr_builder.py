# gui/content/3_bcr_builder.py

import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from pathlib import Path
import os
import glob

from abrat.gui.gui_shared import apply_theme_configurations, folder_selectbox, reset_page_initialization

from abrat.core.utils import date_stamp
from abrat.core.compile_bcr_sequences import compile_bcrs

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

#other pages
other_pages = [page for page in st.session_state.page_initialized.keys() if page != page_name]

hc_column = "HEAVY_CHAIN"
lc_column = "LIGHT_CHAIN"
kappa_column = "KAPPA_CHAIN"
lambda_column = "LAMBDA_CHAIN"
qc_column = "QCHECK_PASSED" #True/False
prod_column = "PRODUCTIVE" #Yes/No
hc_found = "HC_FOUND"
kc_found = "KC_FOUND"
lc_found = "LC_FOUND"
cohort_column = ('SAMPLE_INFORMATION','COHORT')
subject_column = ('SAMPLE_INFORMATION','SUBJECT')

# activate themes for graphical output
apply_theme_configurations()

# ========================================
# Default settings for Repertoire builder
# ========================================

pd.set_option('future.no_silent_downcasting', True)

default_collision_cutoff = 3
default_setting_2 = False

if not st.session_state.sequence_xlsx_verified:
    st.session_state.sequence_xlsx_to_combine = False

hc_identifier = 'HC'
kc_identifier = 'KC'
lc_identifier = 'LC'

# Colors for plots
color_map = {
        'Unpaired heavy chains': "#4d4d4d",        # dunkles Grau
        'Paired heavy/kappa chains': "#4682B4",      # steelblue
        'Paired heavy/lambda chains': "#B0C4DE",     # helles steelblue
        'Double positives (HKL)': "#FA8072",          # salmon
        'Two light chains (KL)': "#FFB19A",          # light salmon
        'Unpaired kappa chains': "#808080",           # mittleres Grau
        'Unpaired lambda chains': "#D3D3D3"           # helles Grau
    }

# =================
# Define functions
# =================

def bcr_df_quality_statistics(bcr_df, hc_col, kc_col, lc_col, qc_col, prod_col):
    """
    Extracts information from bcr_df and returns an Altair bar chart that shows
    the total number of quality passed and productive sequences for heavy and light chains.
    The legend is interactive: when categories (groups) are selected/deselected in the legend,
    both the bars and the aggregated total update dynamically.

    The tooltip displays the Group, Count, and Percentage (of the row's total).
    """

    # Summarize the data
    summary_data = {
        "Passed, productive": [
            sum((bcr_df[(hc_col, qc_col)] == True) & (bcr_df[(hc_col, prod_col)] == "Yes")),
            sum((bcr_df[(kc_col, qc_col)] == True) & (bcr_df[(kc_col, prod_col)] == "Yes")),
            sum((bcr_df[(lc_col, qc_col)] == True) & (bcr_df[(lc_col, prod_col)] == "Yes"))
        ],
        "Passed, non-productive": [
            sum((bcr_df[(hc_col, qc_col)] == True) & (bcr_df[(hc_col, prod_col)] == "No")),
            sum((bcr_df[(kc_col, qc_col)] == True) & (bcr_df[(kc_col, prod_col)] == "No")),
            sum((bcr_df[(lc_col, qc_col)] == True) & (bcr_df[(lc_col, prod_col)] == "No"))
        ],
        "Not passed, productive": [
            sum((bcr_df[(hc_col, qc_col)] == False) & (bcr_df[(hc_col, prod_col)] == "Yes")),
            sum((bcr_df[(kc_col, qc_col)] == False) & (bcr_df[(kc_col, prod_col)] == "Yes")),
            sum((bcr_df[(lc_col, qc_col)] == False) & (bcr_df[(lc_col, prod_col)] == "Yes"))
        ],
        "Not passed, non-productive": [
            sum((bcr_df[(hc_col, qc_col)] == False) & (bcr_df[(hc_col, prod_col)] == "No")),
            sum((bcr_df[(kc_col, qc_col)] == False) & (bcr_df[(kc_col, prod_col)] == "No")),
            sum((bcr_df[(lc_col, qc_col)] == False) & (bcr_df[(lc_col, prod_col)] == "No"))
        ],
    }

    # Create a DataFrame with rows for Heavy Chains and Light Chains
    df = pd.DataFrame(summary_data, index=['Heavy Chains', 'Kappa Chains', 'Lambda Chains'])
    # Convert to long format for Altair
    df_long = df.reset_index().melt(id_vars='index', var_name='Group', value_name='Count')

    # Define an interactive multi-selection bound to the legend.
    # With empty='all', if nothing is selected, all data is shown.
    selection = alt.selection_point(fields=['Group'], bind='legend', empty=True)

    # Define the color encoding
    color = alt.Color(
        'Group:N',
        scale=alt.Scale(
            domain=[
                "Passed, productive",
                "Passed, non-productive",
                "Not passed, productive",
                "Not passed, non-productive"
            ],
            range=[
                "#caebca",  # for "Passed, productive"
                "#ffa95b",  # for "Passed, non-productive"
                "#FFB19A",  # for "Not passed, productive"
                "#FA8072"  # for "Not passed, non-productive"
            ]
        ),
        legend=alt.Legend(
            title="Group (Shift+click for multiple selections)",
            titleLimit=1000,
            orient='bottom',
            direction='vertical',  # Arrange items vertically
            columns=1,
            titleAnchor='start',
            labelLimit=200
        )
    )

    # Create the bar chart.
    # First, filter the data according to the selection.
    # Then use transform_joinaggregate to compute the total Count per row (grouped by 'index'),
    # and transform_calculate to compute the percentage for each row.
    bars = alt.Chart(df_long).transform_filter(selection) \
        .transform_joinaggregate(total='sum(Count)', groupby=['index']) \
        .transform_calculate(percentage='datum.Count / datum.total') \
        .mark_bar(orient='horizontal') \
        .encode(
        y=alt.Y('index:N', axis=alt.Axis(title=None)),
        x=alt.X('Count:Q',
                stack='zero',
                title="Sum",
                scale=alt.Scale(domain=[0, df_long['Count'].max() * 1.1])
                ),
        color=color,
        tooltip=[
            alt.Tooltip('Group:N', title="Group"),
            alt.Tooltip('Count:Q', title="Count"),
            alt.Tooltip('percentage:Q', format=".1%", title="Percent")
        ]
    ).add_params(selection)

    # Create the central text that shows the aggregated sum per row.
    # We aggregate the total Count per row (Heavy and Light) from the filtered data.
    sum_text = alt.Chart(df_long).transform_filter(selection) \
        .transform_aggregate(total='sum(Count)', groupby=['index']) \
        .mark_text(align='left', dx=3, color='black') \
        .encode(
        y=alt.Y('index:N'),
        x=alt.X('total:Q', title="Sum"),
        text=alt.Text('total:Q', format='d')
    )

    # Combine the layers: bars plus aggregated total text.
    final_chart = (bars + sum_text).properties(
        width=200,
        height=300
    )

    return final_chart

def reset_verification():
    """
    Reset session state for files.
    """
    st.session_state.sequence_xlsx_verified = False
    st.session_state.sequence_xlsx_to_combine = False
    # st.session_state.bcr_results = False

def merge_sequence_excel_files(paths_to_files):
    """
        Check input files for duplicates
        :return:
        """
    if len(paths_to_files) > 1:
    # check if sequence names are identical in different lists
        return pd.concat([pd.read_excel(path) for path in paths_to_files], ignore_index=True)
    else:
        return pd.read_excel(paths_to_files[0])
    pass

def check_input_files(input_folder):
    """
    Check input files before loading
    :return:
    """
    try:
        files = [p for p in Path(input_folder).rglob("*"+st.session_state.annotations_file_name)]
        if not files:
            raise FileNotFoundError("No '"+st.session_state.annotations_file_name+"'-files found in "+input_folder)
        st.session_state.sequence_xlsx_verified = True
        #check_files_for_duplicates(input_folder)
        st.success("Files successfully evaluated!")
        return files
    except FileNotFoundError as e:
        st.error(f"Error: {e}")
        return False
    except ValueError as e:
        st.error(f"Formatting error: {e}")
        return False

def get_productive_chain_pairing_counts(df):
    """Takes a bcr dataframe with
    hc_column = "HEAVY_CHAIN"
    kappa_column = "KAPPA_CHAIN"
    lambda_column = "LAMBDA_CHAIN"
    qc_column = "QCHECK_PASSED" #True/False
    prod_column = "PRODUCTIVE" #Yes/No

    and filters for productive and qc-passed sequences to generate a statistics dictionary.

    Returns  prod_pairing_stats and all_passed_sequences
    """
    passed_hcs = df[(hc_column,prod_column)].replace({'Yes': True, 'No': False, 'N/A':np.nan}).fillna(False).astype(bool)&\
                 df[(hc_column, qc_column)]
    passed_kcs = df[(kappa_column, prod_column)].replace({'Yes': True, 'No': False, 'N/A':np.nan}).fillna(False).astype(bool)&\
                 df[(kappa_column, qc_column)]
    passed_lcs = df[(lambda_column, prod_column)].replace({'Yes': True, 'No': False, 'N/A':np.nan}).fillna(False).astype(bool)&\
                 df[(lambda_column, qc_column)]
    prod_pairing_stats = {
        'Unpaired heavy chains': sum(passed_hcs&~passed_kcs&~passed_lcs),
        'Unpaired kappa chains': sum(~passed_hcs&passed_kcs&~passed_lcs),
        'Unpaired lambda chains': sum(~passed_hcs&~passed_kcs&passed_lcs),
        'Two light chains (KL)': sum(~passed_hcs&passed_kcs&passed_lcs),
        'Paired heavy/kappa chains': sum(passed_hcs&passed_kcs&~passed_lcs),
        'Paired heavy/lambda chains': sum(passed_hcs&~passed_kcs&passed_lcs),
        'Double positives (HKL)': sum(passed_hcs&passed_kcs&passed_lcs)
    }

    total_sum = sum(prod_pairing_stats.values())

    return prod_pairing_stats, total_sum

def plot_productive_chain_pairing_donut(prod_pairing_stats, color_map):
    """
    Erzeugt einen interaktiven Donut-Chart mit Legende. Beim Hover über Segmente
    werden sowohl der Count als auch der prozentuale Anteil angezeigt.
    Wenn nichts selektiert ist, zeigt der Chart die Gesamtdaten.
    """
    # DataFrame erstellen
    data = pd.DataFrame({
        'Group': list(prod_pairing_stats.keys()),
        'Count': list(prod_pairing_stats.values())
    })

    # Interaktive Mehrfachselektion, die an die Legende gebunden ist.
    # Mit empty='all' werden bei keiner Selektion alle Daten angezeigt.
    selection = alt.selection_point(fields=['Group'], bind='legend', empty=True)

    # Donut-Chart:
    # Zuerst filtern wir die Daten mit der Selektion.
    # Dann aggregieren wir (mit joinaggregate) den Gesamtwert (total) der Count-Spalte
    # und berechnen anschließend den Prozentwert für jede Zeile.
    donut = alt.Chart(data).transform_filter(selection) \
        .transform_joinaggregate(total='sum(Count)') \
        .transform_calculate(percentage='datum.Count / datum.total') \
        .mark_arc(innerRadius=50, outerRadius=80) \
        .encode(
        theta=alt.Theta(field='Count', type='quantitative'),
        color=alt.Color(
            'Group:N',
            scale=alt.Scale(
                domain=list(prod_pairing_stats.keys()),
                range=[color_map[key] for key in prod_pairing_stats.keys()]
            ),
            legend=alt.Legend(
                title="Group (Shift+click for multiple selections)",
                titleLimit=1000,
                orient='bottom',
                direction='vertical',  # Items vertikal anordnen
                columns=2,
                symbolType='square'
            )
        ),
        tooltip=[
            alt.Tooltip('Group:N', title="Group"),
            alt.Tooltip('Count:Q', title="Count"),
            alt.Tooltip('percentage:Q', format=".1%", title="Percent")
        ]
    ).add_params(selection)

    # Zentraler Text: Aggregiere die Summe der selektierten Daten
    selected_text = alt.Chart(data).transform_filter(selection) \
        .transform_aggregate(total='sum(Count)') \
        .mark_text(
        size=24,
        align='center',
        baseline='middle',
        color='black'
    ).encode(
        text=alt.Text('total:Q', format='d')
    )

    # Layering: Donut-Chart plus zentraler Text
    chart = alt.layer(
        donut,
        selected_text
    ).properties(
        width=300,
        height=300
    )

    return chart

# =======================
# Page content - Sidebar
# =======================

# Sidebar
st.sidebar.header("Settings", divider="rainbow")

# Get Project Name and Folder
project_name = st.sidebar.text_input("Project name", st.session_state.project_name, placeholder="Enter a project name")
st.session_state.project_name = project_name

folder_selectbox(page_name, base_user_path, default_data_folder, reset_verification)
data_folder = st.session_state.current_paths[page_name]

output_folder = data_folder

# File selection
excel_files=[]
selected_files=[]

pre_filter = st.sidebar.toggle("Only show campatible 'all-sequences' files", value=True)
if os.path.isdir(data_folder):
    # pre-filter for standard name
    if pre_filter:
        # display only all-sequence files
        excel_files = glob.glob(os.path.join(data_folder, "*all-sequences.xlsx"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*all-sequences.xls")))
    else:
        # display all excel files
        excel_files = glob.glob(os.path.join(data_folder, "*.xls"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*.xlsx")))

    excel_files = [f for f in excel_files if not os.path.basename(f).startswith("~$")]

    if excel_files:
        selected_files = st.sidebar.multiselect(
            "Select Excel files for analysis:",
            options=excel_files,
            default=[f for f in excel_files if "all-sequences.xls" in f],
            format_func=os.path.basename
        )
    else:
        st.sidebar.warning("No Excel files found in the folder..")
else:
    st.sidebar.warning("Data folder not found..")

st.sidebar.divider()
# Set Cutoffs
with st.sidebar.expander("Advanced settings", icon="⚙️"):
    collision_cutoff = st.slider("Collision cut off", 0, 30, value = default_collision_cutoff, key="collision_cut_off")
    setting_2 = st.toggle("Setting_2", value=default_setting_2, key="setting_2")

# Summary of settings as dict
settings= {
    "Collision cutoff": collision_cutoff,
    "Setting 2": setting_2
    }

# =========================
# Page content - Main page
# =========================

# Instructions
st.title("Compile B Cell Receptors from Annotated Sequences")

st.markdown("""
### Module Workflow:
- Set the path to your ***all-sequences.xlsx** files 
- Select the files you want to include _(by default, all *all-sequences*-files are pre-selected)_  
- Press "Compile" to start combination of heavy and light chains 

_Note: Compiled BCR data will be saved to the selected data folder and automatically kept in memory for downstream modules_
<style> div.stButton {text-align:center}</style>""", unsafe_allow_html=True)

st.subheader("Data:")

if selected_files:
    st.session_state.sequence_xlsx_verified = True
    st.session_state.sequence_xlsx_to_combine = selected_files
    st.info("**Selected file(s):**  \n" + "  \n".join(
        [os.path.basename(f) for f in st.session_state.sequence_xlsx_to_combine]))
else:
    st.warning("No files selected")

if st.button("Compile", disabled=not st.session_state.sequence_xlsx_verified):
    with st.spinner("Combining best reads for heavy and light chains..."):
        try:
            output_path = os.path.join(str(output_folder), str(date_stamp(project_name) +"_"+st.session_state.bcr_file_name))

            # Check if output folder is available
            os.makedirs(output_folder, exist_ok=True)

            # make combined df
            combined_df = merge_sequence_excel_files(st.session_state.sequence_xlsx_to_combine)
            bcr_df = compile_bcrs(combined_df, hc_identifier, kc_identifier, lc_identifier, collision_cutoff)

            with pd.ExcelWriter(output_path) as writer:
                # all B cell infos
                bcr_df.to_excel(writer, sheet_name='Combined_Chains')
                # make something different to demultiplex again?!
                bad_midi_df = bcr_df[bcr_df[('SAMPLE_INFORMATION', 'WARNING')].str.contains('MIDI')]
                bad_midi_df.to_excel(writer, sheet_name='Bad_Midis')
                bad_fwr4_df = bcr_df[bcr_df[('SAMPLE_INFORMATION', 'WARNING')].str.contains('FWR4')]
                bad_fwr4_df.to_excel(writer, sheet_name='Bad_FWR4s')

                unproductives_df = bcr_df[
                    ((bcr_df[('HEAVY_CHAIN', 'PRODUCTIVE')] == "No") & (bcr_df[('HEAVY_CHAIN', 'QCHECK_PASSED')])) |
                    ((bcr_df[('LIGHT_CHAIN', 'PRODUCTIVE')] == "No") & (bcr_df[('LIGHT_CHAIN', 'QCHECK_PASSED')]))
                ]
                unproductives_df.to_excel(writer, sheet_name='Final_Unproductives')

                # productive and complete only
                complete_only = bcr_df[(bcr_df[('HEAVY_CHAIN', 'PRODUCTIVE')] == "Yes") &
                                       (bcr_df[('HEAVY_CHAIN', 'QCHECK_PASSED')]) &
                                       (bcr_df[('LIGHT_CHAIN', 'PRODUCTIVE')] == "Yes") &
                                       (bcr_df[('LIGHT_CHAIN', 'QCHECK_PASSED')])
                                       ]

                complete_only = complete_only.loc[:,
                                pd.IndexSlice[['SAMPLE_INFORMATION', 'HEAVY_CHAIN', 'LIGHT_CHAIN'], :]]
                complete_only.to_excel(writer, sheet_name='Complete_Abs')
            bcr_success = True

        except Exception as e:
            bcr_success = False
            st.error(f"Error while building repertoire: {e}")
    if bcr_success:
        # make statistics from rep_df of quality check
        qc_df = bcr_df_quality_statistics(bcr_df, hc_column, kappa_column, lambda_column, qc_column, prod_column)

        # make statistics for paired productive and qc_passed chains
        chain_pairs_statistics, all_passed = get_productive_chain_pairing_counts(bcr_df)

        # save qc_df for plots
        st.session_state.bcr_results = {'bcr_df':bcr_df,
                                        'bcr_file_name':os.path.basename(output_path),
                                        'qc_df':qc_df,
                                        'chain_pairs_statistics':chain_pairs_statistics,
                                        'all_passed':all_passed
                                        }
        st.success('DONE.')

    # set latest data folder
    st.session_state.current_paths['latest_data_path'] = output_folder

    # reset all other page initializations to preload latest data folder
    reset_page_initialization(other_pages)

# display the most recent summary data (data actively saved)
if st.session_state.bcr_results and "qc_df" in st.session_state.bcr_results:
    st.subheader("Results:")
    st.write(st.session_state.bcr_results['bcr_file_name'])
    left_result, right_result = st.columns(2)
    combined_stat_df = st.session_state.bcr_results['qc_df']
    with left_result:
        st.caption("Quality of selected heavy and light chains:")
        st.altair_chart(st.session_state.bcr_results['qc_df'], use_container_width=True)
    with right_result:
        st.caption("Numbers of paired/unpaired productive chains:")
        # display number of productive heavy and light chain pairs
        donut_plot = plot_productive_chain_pairing_donut(st.session_state.bcr_results['chain_pairs_statistics'],
                                                         color_map)
        st.altair_chart(donut_plot, use_container_width=True)
        #st.table(st.session_state.bcr_results['chain_pairs_statistics'])
    # display the latest combined bcr_df
    with st.expander("🔍 Show combined heavy and light chains"):
        st.dataframe(st.session_state.bcr_results['bcr_df'])

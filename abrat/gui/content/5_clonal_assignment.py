import os, glob
import streamlit as st
import pandas as pd
from pathlib import Path
import concurrent.futures
import datetime

from abrat.gui.gui_shared import folder_selectbox, reset_page_initialization
from abrat.core.utils import timestamp_filename, dict_to_str, write_clustering_log, plot_pie_chart
from abrat.core.clonal_assignment import (process_bcr_list,
                                          save_excel_with_row_colors,
                                          plot_clonality_donut_chart,
                                          clone_subcol,
                                          clone_color_subcol)

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

# other pages
other_pages = [page for page in st.session_state.page_initialized.keys() if page != page_name]

# ===========================================
# Default settings for basic clonal analysis
# ===========================================

default_number_of_iterations = {
"Iterative CDR3 similarity":20,
"Matrix CDR3 similarity":10,
"Hierarchical CDR3 Clustering":1
}

default_hc_minimum_identity_aa_fixed = 75
default_lc_minimum_identity_aa_fixed = 50

# allowed length difference
default_cdrh3_aa_length_difference_threshold = 1
default_cdrl3_aa_length_difference_threshold = 1

default_collision_cutoff = 3

info_text = []
log_file_message = []

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
cdr3_nt_subcol = 'CDR3_NT'

settings_to_columns = {
        'V': top_v_gene_subcol,
        'D': top_d_gene_subcol,
        'J': top_j_gene_subcol
    }

clone_color_dict = {'palette': 'GnBu_d',
                    'non_clonal': '#D3D3D3',  # light gray
                    'undefined': '#FFFFFF'    # white
                   }

# =================
# Define functions
# =================

def reset_verification():
    """
    Resets the session state for file verification.

    This function is intended to reset any session state flags related to file verification.
    Currently, the relevant session state updates are commented out.

    Returns:
        None
    """
    # st.session_state.sequence_xlsx_verified = False
    # st.session_state.sequence_xlsx_to_combine = False
    pass


def check_input_files(input_folder, search_string):
    """
    Checks for input files in the specified folder that match a given search string.

    This function recursively searches the specified folder for files whose names contain the provided
    search string. If matching files are found, the function sets the session state flag 'bcr_xlsx_verified'
    to True, displays a success message, and returns a list of Path objects representing the found files.
    If no matching files are found or if an error occurs, it displays an error message and returns False.

    Parameters:
        input_folder (str or Path): The folder in which to search for input files.
        search_string (str): The string to search for in the file names.

    Returns:
        list or bool: A list of Path objects for the matching files if found; otherwise, False.
    """
    try:
        files = [p for p in Path(input_folder).rglob("*" + search_string)]
        if not files:
            raise FileNotFoundError("No '" + search_string + "' files found in " + str(input_folder))
        st.session_state.bcr_xlsx_verified = True
        # Optionally, check for duplicates here.
        st.success("Files found!")
        return files
    except FileNotFoundError as e:
        st.error(f"Error: {e}")
        return False
    except Exception as e:
        st.error(f"Exception: {e}")
        return False

# =======================
# Page content - Sidebar
# =======================

# Sidebar
st.sidebar.header("Settings", divider="rainbow")

# Datafolder selection
folder_selectbox(page_name, base_user_path, default_data_folder)
data_folder = st.session_state.current_paths[page_name]
output_folder = os.path.join(data_folder, 'clustered_repertoires')

# File selection
excel_files=[]
selected_files=[]

pre_filter = st.sidebar.toggle("Only show compatible 'b-cell-receptors' files", value=True)
if os.path.isdir(data_folder):
    # pre-filter for standard name
    if pre_filter:
        # display only all-sequence files
        excel_files = glob.glob(os.path.join(data_folder, "*b-cell-receptors*.xlsx"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*b-cell-receptors*.xls")))
    else:
        # display all excel files
        excel_files = glob.glob(os.path.join(data_folder, "*.xls"))
        excel_files.extend(glob.glob(os.path.join(data_folder, "*.xlsx")))

    excel_files = [f for f in excel_files if not os.path.basename(f).startswith("~$")]

    if excel_files:
        selected_files = st.sidebar.multiselect(
            "Select Excel files for analysis",
            options=excel_files,
            default=[f for f in excel_files if "b-cell-receptors" in f],
            format_func=os.path.basename
        )
    else:
        st.sidebar.warning("No Excel files found in the folder..")
else:
    st.warning("Data folder not found..")

if st.sidebar.button("Load BCR excel-file(s)", disabled=False if selected_files else True):
    # reset loaded data:
    st.session_state.bcr_results['filtered_bcrs'] = {}

    header_cols = [0, 1]
    index_cols = [0]
    try:
        for file in selected_files:
            file_name = os.path.basename(file)
            st.session_state.bcr_results['filtered_bcrs'][file_name] = pd.read_excel(file, header=header_cols, index_col=index_cols)
        st.sidebar.success("Excel file(s) successfully loaded!")
    except Exception as e:
        st.sidebar.error(f"Error while loading excel file(s): {e}")

st.sidebar.divider()

# Set Cutoffs
st.sidebar.caption("Select chains for clustering")
setting_expander = st.sidebar.expander("Clustering settings", icon="⚙️")

hc_selected = st.sidebar.checkbox('Heavy Chain', value=True)

if hc_selected:
    # Selection of VDJ for pre-grouping
    # hc_setting_expander = st.sidebar.expander("Heavy chain settings", icon="⚙️")
    with st.sidebar.expander("Heavy chain settings", icon="⚙️"):
        hc_vdj_selection = st.segmented_control('Gene segments for clustering',
                                                    ['V','D','J'],
                                                    default=['V','J'],
                                                    selection_mode='multi',
                                                   help='Gene segments that need to match between sequences to be '
                                                        'considered as clonal')

        heavy_chain_settings = {'group_by': [settings_to_columns[gene_segment]
                                             for gene_segment in hc_vdj_selection]}

        # selection of clustering algorithm heavy chains
        hc_algorithm = st.radio("Clustering algorithm",
                 key="hc_cluster",
                 options=["Iterative CDR3 similarity",
                          "Matrix CDR3 similarity",
                          "Hierarchical CDR3 Clustering"],
                 help="Clustering algorithm for CDR3s")

        heavy_chain_settings['algorithm']=hc_algorithm
        heavy_chain_settings['params'] = {}

        if hc_algorithm: #== "Iterative CDR3 similarity" or hc_algorithm == "Matrix CDR3 similarity":
            if st.toggle('Restrict CDR3 length difference for clustering',
                         value=True,
                         key='hc_length_difference',
                         help='Sequences will only be considered clonal, if the difference in their CDR3 lengths is below '
                              'the slider-specified number of amino acids'):
                heavy_chain_settings['params']['length_threshold'] = st.slider("Max. difference", 0, 20,
                                                                               value=default_cdrh3_aa_length_difference_threshold,
                                                                               key="cdrh3_length_difference_threshold",
                                                                               help='Define the maximum amino acid length difference '
                                                                                    'for CDR3s to be considered as clonal')
            else:
                heavy_chain_settings['params']['length_threshold'] = None
            heavy_chain_settings['params']['lev_threshold'] = (100 - st.slider("Minimum % CDRH3 identity for clustering", 0, 100,
                                                         value=default_hc_minimum_identity_aa_fixed,
                                                         step=1, key="fixed_cdrh3_aa_identity",
                                                         help="Similarity cutoff for heavy chain CDR3s to be considered as clonal")) / 100

            heavy_chain_settings['params']['iterations'] = st.slider("Number of iterations", 1, 100,
                                                                     value=default_number_of_iterations[hc_algorithm],
                                                                     key="number_of_hc_iterations",
                                                                     help="Number of repeated clonal assignments to "
                                                                          "reduce seeding effects")

else:
    heavy_chain_settings = {}

lc_selected = st.sidebar.checkbox('Light Chain', value=False)

if lc_selected:
    with (st.sidebar.expander("Light chain settings", icon="⚙️")):
        lc_vj_selection = st.segmented_control('Gene segments for clustering',
                                                   ['V', 'J'],
                                                   default=['V'],
                                                   selection_mode='multi',
                                                   help="Gene segments that need to match between sequences to be "
                                                        "considered as clonal"
                                               )
        light_chain_settings = {
            'group_by': [settings_to_columns[gene_segment]
                         for gene_segment in lc_vj_selection] if lc_vj_selection else []
        }

        # selection of clustering algorithm for light chains
        lc_algorithm = st.radio("Clustering algorithm",
                 key="lc_cluster",
                 options=["Iterative CDR3 similarity",
                          "Matrix CDR3 similarity",
                          "Hierarchical CDR3 Clustering"],
                 help="Clustering algorithm for CDR3s")
        light_chain_settings['lc_algorithm'] = lc_algorithm
        light_chain_settings['params'] = {}

        if lc_algorithm: # == "Iterative CDR3 similarity" or lc_algorithm == "Matrix CDR3 similarity":
            if st.toggle('Restrict CDR3 length difference for clustering',
                         value=True,
                         key="lc_length_difference",
                         help='Sequences will only be considered clonal, if the difference in their CDR3 lengths is below '
                              'the slider-specified number of amino acids'):
                light_chain_settings['params']['length_threshold'] = st.slider("Max. difference", 0, 10,
                                      value = default_cdrl3_aa_length_difference_threshold,
                                      key="cdrl3_length_difference_threshold",
                                      help='Define the maximum amino acid length difference '
                                      'for CDR3s to be considered as clonal')
            else:
                light_chain_settings['params']['length_threshold'] = None

            light_chain_settings['params']['lev_threshold'] = (100 - st.slider(
                "Minimum % CDRL3 identity for clustering", 1, 100,
                value=default_lc_minimum_identity_aa_fixed,
                step=1, key="fixed_cdrl3_aa_identity",
                help="Similarity cutoff for light chain CDR3s to be considered as clonal")) / 100

            light_chain_settings['params']['iterations'] = st.slider("Number of iterations", 0, 100,
                                      value = default_number_of_iterations[lc_algorithm],
                                      key="number_of_lc_iterations",
                                      help="Number of repeated clonal assignments to "
                                           "reduce seeding effects"
                                      )

else:
    light_chain_settings = {}

if not hc_selected and not lc_selected:
    st.sidebar.warning('At least one chain must be selected.')

# Summary of settings as dict
settings= {'heavy_chain': heavy_chain_settings,
           'light_chain': light_chain_settings}

# =========================
# Page content - Main page
# =========================

# Instructions
st.title("Clonal Clustering")

st.markdown("""
### Module Workflow:
- Use data in memory from _'Filter & Split B-Cell Receptors'_ or select your filtered ***b-cell-receptors** files
- Check your settings.
- Press "Cluster" to start clonal analysis.
<style> div.stButton {text-align:center} </style>""",
            unsafe_allow_html=True)

st.subheader("Data:")

if not st.session_state.bcr_results['filtered_bcrs']:
    st.warning("No data loaded yet.")

    st.markdown("""
        - Set the path to your filtered ***b-cell-receptors** files
        - Select the files you want to filter
        - Press "Load" to load the data
        """)
else:
    st.info("**Selected file(s):**  \n" + "  \n".join(
        [os.path.basename(file) for file in st.session_state.bcr_results['filtered_bcrs'].keys()])
            )

    # Settings visualization
    st.subheader("Settings:")
    with st.expander('Clustering settings', expanded=True):
        settings_lines = dict_to_str(settings, indent=0)
        settings_text = "\n".join(settings_lines)
        st.markdown("```\n" + settings_text + "\n```")#, use_container_width=True)

    # ===========
    # Clustering
    # ===========

    if st.button("Perform clonal analysis", disabled=False if st.session_state.bcr_results['filtered_bcrs'] else True):
        # Check if output folder is available
        os.makedirs(output_folder, exist_ok=True)

        bcr_dict = st.session_state.bcr_results['filtered_bcrs']
        total = len(bcr_dict)
        results = {}
        errors = {}  # Hier speichern wir Fehlermeldungen pro Datei

        # Widgets für Info-Text und Fortschrittsbalken
        progress_text = st.empty()
        progress_bar = st.progress(0)
        progress_text.text(f"0/{total} files processed, {len(results)}/{total} successfully evaluated.")

        # Erstelle ein Mapping von Future zu Dateinamen
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
            future_to_file = {
                executor.submit(process_bcr_list, df, settings): filename
                for filename, df in bcr_dict.items()
            }

            while future_to_file:
                # Wait at least 0.1s.
                done, _ = concurrent.futures.wait(
                    future_to_file, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in done:
                    filename = future_to_file.pop(future)
                    try:
                        result_df = future.result()
                        results[timestamp_filename(filename, 'clustered')] = result_df
                    except Exception as e:
                        st.error(f"Error processing file {filename}: {e}")
                        errors[filename] = str(e)

                processed = total - len(future_to_file)
                progress_text.text(f"{processed}/{total} files processed, {len(results)}/{total} successfully evaluated.")
                progress_bar.progress(processed / total)

        if errors:
            st.error("Errors during clustering!")

        if results:
            st.session_state.bcr_results['clustered_bcrs'] = {}
            for filename, df in results.items():
                save_excel_with_row_colors(df, os.path.join(output_folder, os.path.basename(filename)))
                #reset and save the latest results in session state
                st.session_state.bcr_results['clustered_bcrs'][os.path.basename(filename)] = {'dataframe':df}
                # automatically plot pie chart and export
                pie_fig = plot_clonality_donut_chart(df, (sample_col, clone_subcol),
                                                     (sample_col, clone_color_subcol))
                # export
                pie_chart_path = os.path.join(output_folder, os.path.basename(filename).split('.')[0]+".pdf")
                pie_fig.savefig(pie_chart_path, bbox_inches='tight')
                # save in session state
                st.session_state.bcr_results['clustered_bcrs'][os.path.basename(filename)]['clonality_pie'] = pie_fig

            timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
            log_file_name = timestamp + "_Clustering-log-file.txt"
            path_to_logfile = os.path.join(output_folder, log_file_name)
            write_clustering_log(path_to_logfile,
                                 settings=settings,
                                 input_files=st.session_state.bcr_results['filtered_bcrs'].keys(),
                                 output_files=st.session_state.bcr_results['clustered_bcrs'].keys()
                                 )

            st.balloons()
            st.success("Clustering completed!")
    # display the latest results
    if st.session_state.bcr_results['clustered_bcrs']:
        st.subheader("Results:")
        for filename, data in st.session_state.bcr_results['clustered_bcrs'].items():
            df = st.session_state.bcr_results['clustered_bcrs'][filename]['dataframe']
            clonality_pie = st.session_state.bcr_results['clustered_bcrs'][filename]['clonality_pie']
            st.caption(filename)
            st.pyplot(clonality_pie)
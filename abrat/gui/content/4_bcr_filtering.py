# gui/content/4_bcr_filtering.py

import streamlit as st
import pandas as pd
import datetime
import os
import glob
from pathlib import Path

from abrat.gui.gui_shared import folder_selectbox, reset_page_initialization

# ==============================================
# Passing of global settings from session state
# ==============================================

page_name = os.path.splitext(os.path.basename(__file__))[0]
base_user_path = st.session_state.global_constants['base_user_path']

# Get default navigation start.
if page_name not in st.session_state.current_paths:
    st.session_state.current_paths[page_name] = base_user_path

# default folders
default_data_folder = st.session_state.current_paths['latest_data_path']
default_output_folder = "filtered_repertoires"

# other pages
other_pages = [page for page in st.session_state.page_initialized.keys() if page != page_name]

# ========================================
# Default settings for Repertoire filter
# ========================================

cohort_col = ('SAMPLE_INFORMATION','COHORT')
subject_col = ('SAMPLE_INFORMATION','SUBJECT')
prod_hc_col = ('HEAVY_CHAIN','PRODUCTIVE')
qc_passed_hc_col = ('HEAVY_CHAIN','QCHECK_PASSED')
prod_lc_col = ('LIGHT_CHAIN','PRODUCTIVE')
qc_passed_lc_col = ('LIGHT_CHAIN','QCHECK_PASSED')

standard_filter_columns = [prod_hc_col,
                           qc_passed_hc_col,
                           prod_lc_col,
                           qc_passed_lc_col
                           ]
standard_split_columns = [cohort_col, subject_col]

default_prod_hc = ["Yes"]
default_qc_passed_hc = [True]
default_prod_lc = ["Yes"]
default_qc_passed_lc = [True]

default_settings = {
    prod_hc_col: default_prod_hc,
    qc_passed_hc_col: default_qc_passed_hc,
    prod_lc_col: default_prod_lc,
    qc_passed_lc_col: default_qc_passed_lc
}

# Initialize bcr_results and filtered_bcr if not done yet
if 'bcr_results' not in st.session_state:
    st.session_state['bcr_results'] = {}
if 'filtered_bcrs' not in st.session_state.bcr_results:
    st.session_state.bcr_results['filtered_bcrs'] = {}

info_text = []
log_file_message = []
default_filter_col = []
default_split_col = []
missing_columns = []
selected_filter_columns = []
selected_split_columns = []
active_filters = {}

delete_message = "  \n=> delete complete BCR if not matching  \n"

# =================
# Define functions
# =================
def delete_alternate_message(col):
    """
    Returns a message instructing to delete only the entries for the specified column level if they do not match.

    Parameters:
        col (list or tuple): A collection where the first element represents the column name.

    Returns:
        str: A formatted message.
    """
    return "  \n=> delete only " + str(col[0]) + " entries, if not matching  \n"


def reset_verification():
    """
    Resets the session state for file verification.

    Currently, the session state variables for file verification are commented out.

    Returns:
        None
    """
    # st.session_state.sequence_xlsx_verified = False
    # st.session_state.sequence_xlsx_to_combine = False
    # st.session_state.bcr_results = False


def get_default_delete_toggle(col):
    """
    Determines the default deletion toggle for a column based on its level 0 name.

    For "HEAVY_CHAIN", the default is to delete the entire row.
    For "LIGHT_CHAIN", "KAPPA_CHAIN", and "LAMBDA_CHAIN", the default is to delete only the specific value.
    For all other columns, the default is to delete the entire row.

    Parameters:
        col (list or tuple): A collection where the first element represents the column name.

    Returns:
        bool: True if the default is to delete the entire row; False if only the individual value should be deleted.
    """
    if col[0] == "HEAVY_CHAIN":
        return True
    elif col[0] in ["LIGHT_CHAIN", "KAPPA_CHAIN", "LAMBDA_CHAIN"]:
        return False
    else:
        return True


def get_group_columns(input_df, col):
    """
    Finds all columns in the DataFrame that share the same level 0 value as the provided column.

    Parameters:
        input_df (pd.DataFrame): The DataFrame containing columns with a MultiIndex.
        col (list or tuple): A column identifier where the first element is used for matching.

    Returns:
        list: A list of columns from the DataFrame with the same level 0 value as col[0].
    """
    return [c for c in input_df.columns if c[0] == col[0]]

# =======================
# Page content - Sidebar
# =======================

# Sidebar
st.sidebar.header("Settings", divider="rainbow")

folder_selectbox(page_name, base_user_path, default_data_folder)
data_folder = st.session_state.current_paths[page_name]

output_folder = os.path.join(data_folder, default_output_folder)

#st.write("Current folder:", data_folder)
excel_files=[]
selected_file=[]
if os.path.isdir(data_folder):
    excel_files = glob.glob(os.path.join(data_folder, "*.xlsx"))
    excel_files.extend(glob.glob(os.path.join(data_folder, "*.xls")))

    excel_files = [f for f in excel_files if not os.path.basename(f).startswith("~$")]

    if excel_files:
        selected_file = st.sidebar.selectbox(
            "Choose a file.",
            options=excel_files,
            format_func=lambda x: os.path.basename(x)
        )

    else:
        st.sidebar.info("No excel files found in folder..")
else:
    st.sidebar.error("Folder does not exist.")

if st.sidebar.button("Load BCR excel-file", disabled=False if selected_file else True):
    try:
        header_cols = [0, 1]
        index_cols = [0]
        df_loaded = pd.read_excel(selected_file, header=header_cols, index_col=index_cols)
        st.session_state['bcr_results']['bcr_df'] = df_loaded
        st.session_state['bcr_results']['bcr_file_name'] = os.path.basename(selected_file)
        st.sidebar.success("Excel file successfully loaded!")
    except Exception as e:
        st.sidebar.error(f"Error while loading excel file: {e}")

st.sidebar.divider()

# =========================
# Page content - Main page
# =========================

st.title("Filter compiled B Cell Receptors")
st.markdown("""
### Module Workflow:
- Use data in memory from _'Compile B-Cell Receptors'_ or load a ***b-cell-receptors.xlsx** file
- Set conditions for filtering and data separation
- Save your filtered and separated dataset(s)
  <style>
  div.stButton {text-align:center}
  </style>""", unsafe_allow_html=True)

st.subheader("Data:")

if not 'bcr_df' in st.session_state['bcr_results']:
    st.warning("No data loaded yet.")

    st.markdown("""
    - Set the path to your ***b-cell-receptors.xlsx** files
    - Select the file you want to filter
    - Press "Load" to load the data
    """)
else:
    st.info("**Selected file:**  \n" + st.session_state.bcr_results['bcr_file_name'])

    with st.expander("🔍 Inspect loaded dataset:", expanded=False):
        st.dataframe(st.session_state['bcr_results']['bcr_df'])

    orig_df = st.session_state['bcr_results']['bcr_df']
    filtered_df = orig_df.copy()
    st.subheader("Filtering and Data separation:")

    default = st.toggle("Use default settings", True)

    if default: # Default filtering
        for default_col in standard_filter_columns:
            if default_col in orig_df.columns:
                selected_filter_columns.append(default_col)
                default_delete = get_default_delete_toggle(default_col)
                # Hier wird auch der Hinweis in active_filters abgespeichert:
                active_filters[default_col] = f"{default_settings[default_col]}" + (
                    delete_message if default_delete else delete_alternate_message(default_col))
            else:
                missing_columns.append(default_col)
        for column in selected_filter_columns:
            selected_options = default_settings[column]
            default_delete = get_default_delete_toggle(column)
            if default_delete:
                filtered_df = filtered_df[filtered_df[column].isin(selected_options)]
            else:
                cond = filtered_df[column].isin(selected_options)
                group_cols = get_group_columns(filtered_df, column)
                filtered_df.loc[~cond, group_cols] = pd.NA

    else: # Custom filtering
        for default_col in standard_filter_columns:
            if default_col in orig_df.columns:
                default_filter_col.append(default_col)
            else:
                missing_columns.append(default_col)
        selected_filter_columns = st.multiselect(
            "**Choose the columns to be filtered:**",
            options=orig_df.columns,
            default=default_filter_col,
            format_func=lambda x: ", ".join([str(a) for a in x])
        )

        # ------------------------------
        # Erweiterte Filter-Section mit zusätzlichem Toggle
        # ------------------------------
        for column in selected_filter_columns:
            col_name = ", ".join(column)
            non_missing = orig_df[column].dropna()
            unique_vals = list(non_missing.unique())
            if orig_df[column].isnull().any():
                unique_vals.append("NaN")

            # Check if all non-NaNs are exclusively 0 or 1 (as a number)
            non_nan_values = [x for x in unique_vals if x != "NaN"]
            if non_nan_values and all(x in {0, 0.0, 1, 1.0} for x in non_nan_values):
                # Replace 1/1.0 by True and 0/0.0 by False
                unique_vals = [True if x in {1, 1.0} else False if x in {0, 0.0} else x for x in unique_vals]

            # Unterscheidung der Filter-Widgets je nach Anzahl der einzigartigen Werte:
            if len(unique_vals) <= 3:
                sorted_options = sorted(unique_vals, key=lambda x: str(x))
                pre_selection = default_settings[column] if column in standard_filter_columns else sorted_options

                col1, col2 = st.columns([3, 1])
                with col1:
                    selected_options = st.pills(
                        f"Filter for '{col_name}':",
                        options=sorted_options,
                        selection_mode="multi",
                        default=pre_selection,
                        key=f"pills_{column}"
                    )
                with col2:
                    default_delete = get_default_delete_toggle(column)
                    delete_complete = st.checkbox("Delete complete BCR", value=default_delete, key=f"delete_{column}")

                # Hier wird zusätzlich der Hinweis angehängt, wenn delete_complete aktiv ist:
                active_filters[column] = f"{selected_options}" +\
                                         (delete_message if delete_complete else delete_alternate_message(column))
                if selected_options:
                    # Erstellen der Filterbedingung
                    cond = filtered_df[column].isin([opt for opt in selected_options if opt != "NaN"])
                    if "NaN" in selected_options:
                        cond |= filtered_df[column].isnull()
                    if delete_complete:
                        filtered_df = filtered_df[cond]
                    else:
                        group_cols = get_group_columns(filtered_df, column)
                        filtered_df.loc[~cond, group_cols] = pd.NA
                else:
                    filtered_df = filtered_df.iloc[0:0]

            elif pd.api.types.is_numeric_dtype(orig_df[column]):
                min_val = float(orig_df[column].min())
                max_val = float(orig_df[column].max())
                col1, col2 = st.columns([3, 1])
                with col1:
                    range_values = st.slider(
                        f"Filter for '{col_name}' (range):",
                        min_value=min_val,
                        max_value=max_val,
                        value=(min_val, max_val),
                        key=f"slider_{column}"
                    )
                with col2:
                    default_delete = get_default_delete_toggle(column)
                    delete_complete = st.checkbox("Delete complete BCR", value=default_delete, key=f"delete_{column}")
                active_filters[column] = f"{range_values}" +\
                                         (delete_message if delete_complete else delete_alternate_message(column))
                if delete_complete:
                    filtered_df = filtered_df[
                        (filtered_df[column] >= range_values[0]) & (filtered_df[column] <= range_values[1])
                    ]
                else:
                    cond = (filtered_df[column] >= range_values[0]) & (filtered_df[column] <= range_values[1])
                    group_cols = get_group_columns(filtered_df, column)
                    filtered_df.loc[~cond, group_cols] = pd.NA

            else:
                unique_options = orig_df[column].unique().tolist()
                pre_selection = default_settings[column] if column in standard_filter_columns else unique_options

                col1, col2 = st.columns([3, 1])
                with col1:
                    selected_options = st.multiselect(
                        f"Filter for '{col_name}':",
                        options=unique_options,
                        default=pre_selection,
                        key=f"multiselect_{column}"
                    )
                with col2:
                    default_delete = get_default_delete_toggle(column)
                    delete_complete = st.checkbox("Delete complete BCR", value=default_delete, key=f"delete_{column}")
                if selected_options:
                    if "NaN" in selected_options:
                        cond = filtered_df[column].isin([opt for opt in selected_options if opt != "NaN"]) | \
                               filtered_df[column].isnull()
                    else:
                        cond = filtered_df[column].isin(selected_options)
                    if delete_complete:
                        filtered_df = filtered_df[cond]
                    else:
                        group_cols = get_group_columns(filtered_df, column)
                        filtered_df.loc[~cond, group_cols] = pd.NA
                else:
                    filtered_df = filtered_df.iloc[0:0]
                active_filters[column] = f"{selected_options}" +\
                                         (delete_message if delete_complete else delete_alternate_message(column))

    with st.expander("🔍 Inspect filtered dataset:", expanded=False):
        st.dataframe(filtered_df)

    # Split datasets by selected columns
    if default:
        for default_col in standard_split_columns:
            if default_col in orig_df.columns:
                selected_split_columns.append(default_col)
            else:
                selected_split_columns.append(default_col)
    else:
        for default_col in standard_split_columns:
            if default_col in orig_df.columns:
                default_split_col.append(default_col)
            else:
                default_split_col.append(default_col)
        selected_split_columns = st.multiselect(
            "**Choose the columns by which the data should be separated** (default: cohort and subject):",
            options=orig_df.columns,
            default=default_split_col,
            format_func=lambda x: ", ".join([str(a) for a in x])
        )

    # add active filters to info text
    info_text.append("**Active filter:**")
    for k, v in active_filters.items():
        info_text.append(str(k) + ": " + str(v))

    # initialize split_df as dictionary
    split_dfs = {"noSplit": filtered_df}
    if selected_split_columns:
        split_dfs = filtered_df.groupby(selected_split_columns)
        unique_groups = list(split_dfs.groups.keys())
        # add splitting infos to text, if dataframe is splitted
        info_text.append("**Subgroup(s) to be exported:** " + str(len(unique_groups)) + "  \n" +
                         "  \n".join([str(x) for x in unique_groups]))

    # check for any missing columns and give warning
    if missing_columns:
        st.warning("**Column(s) not found:**  \n" + "  \n".join([str(x) for x in missing_columns]), icon="⚠️")

    # print info text
    if info_text:
        st.info("  \n".join(info_text))

    # Export data
    if st.button("Save filtered Excel file(s)"):
        # Check if output folder is available
        os.makedirs(output_folder, exist_ok=True)

        # Reset filtered dataframes in memory
        st.session_state.bcr_results['filtered_bcrs'] = {}

        timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
        if 'bcr_file_name' in st.session_state['bcr_results']:
            original_name = Path(st.session_state['bcr_results']['bcr_file_name']).stem
        else:
            original_name = "NewData"

        # generate log file message
        log_file_name = timestamp + "_FILTER-LOG_" + original_name + ".txt"
        log_file_message.append("TIMESTAMP: " + timestamp + "\n")
        log_file_message.append("\n".join(info_text).replace('**',''))

        # generate and save splitted dataframes
        exported_files = []
        for group, df in split_dfs:
            file_name = f"{original_name}_{timestamp}_filtered_{'_'.join(list(group))}.xlsx"
            file_path = os.path.join(output_folder, file_name)

            try:
                with pd.ExcelWriter(file_path) as writer:
                    df.to_excel(writer, sheet_name='Filtered Data')
                    st.session_state.bcr_results['filtered_bcrs'][file_path] = df
                exported_files.append(file_name)
            except Exception as e:
                st.error(f"Error while saving: {e}")
        if exported_files:
            st.success(f"Filtered Excel file(s) saved to folder {output_folder}:  \n" + "  \n".join(exported_files))
            # add exported files to log file
            log_file_message.append("\n**Exported files:**")
            log_file_message.append("\n".join(exported_files))

            # save output folder
            st.session_state.current_paths['latest_data_path'] = output_folder
        # save log file message
        with open(os.path.join(output_folder, log_file_name), 'w') as outfile:
            outfile.write("\n".join(log_file_message))

        # reset other pages to pass latest data folder
        reset_page_initialization(other_pages)
# abrat-app.py
import os
import streamlit as st
from importlib.resources import files

def main():
    # ================
    # Set Page Config
    # ================

    st.set_page_config(
        layout='centered',
        page_title='AbRAT: The Antibody Repertoire Analysis Toolkit',
        menu_items={
            'About': '''**AbRAT v1.0.0**        
                The Antibody Repertoire Analysis Toolkit.
                © 2025 Christoph Kreer
                '''
        }
    )

    # =========
    # Set logo
    # =========

    logo_path = files('abrat.gui.assets') / 'abrat_logo.png'
    st.logo(str(logo_path), size="large")

    # ==========================
    # Initialize session states
    # ==========================

    ### Standard names for output data ###

    if 'annotations_file_name' not in st.session_state:
        st.session_state.annotations_file_name = "all-sequences.xlsx"

    if 'bcr_file_name' not in st.session_state:
        st.session_state.bcr_file_name = "b-cell-receptors.xlsx"

    ### Checks for file verification to unlock analysis buttons ###

    if 'ab1files_verified' not in st.session_state:
        st.session_state.ab1files_verified = False  # checks, if ab1-files for ab1_analyzer have been verified

    if 'sequence_xlsx_verified' not in st.session_state:
        st.session_state.sequence_xlsx_verified = False  # checks, if all-sequence.xlsx for bcr_builder are available

    if 'ab1file_info' not in st.session_state:
        st.session_state.ab1file_info = {}  # Contains all infos from ab1 files

    if 'project_name' not in st.session_state:
        st.session_state.project_name = "NewProject"  # To set the project name globally and pass it between scripts

    if 'ab1_results_path' not in st.session_state:
        st.session_state.ab1_results_path = False

    if 'annotations_info' not in st.session_state:
        st.session_state.annotations_info = {}

    if 'sequence_xlsx_to_combine' not in st.session_state:
        st.session_state.sequence_xlsx_to_combine = False

    if 'bcr_xlsx_verified' not in st.session_state:
        st.session_state.bcr_xlsx_verified = False

    if 'bcr_xlsx_to_analyze' not in st.session_state:
        st.session_state.bcr_xlsx_to_analyze = False

    # saves the latest bcr_df
    if 'bcr_results' not in st.session_state:
        st.session_state.bcr_results = {} # will contain the latest bcrs_builder results
    if 'filtered_bcrs' not in st.session_state.bcr_results:
        st.session_state.bcr_results['filtered_bcrs'] = {} # will contain the latest filtered and splitted bcrs
    if 'clustered_bcrs' not in st.session_state.bcr_results:
        st.session_state.bcr_results['clustered_bcrs'] = {} # will contain the latest clustered bcrs

    ### Some global constants

    # Define the allowed base directory and standard output path
    base_path = os.path.abspath(os.path.join("data", "userdata"))
    data_output_path = os.path.join(base_path, "output")

    # Define pages folder and get all pages names
    pages_folder = os.path.join("abrat", "gui", "content")
    pages = [os.path.splitext(os.path.basename(file))[0] for file in os.listdir(pages_folder) \
             if file.endswith(".py") and file != "__init__.py"]

    if 'global_constants' not in st.session_state:
        st.session_state.global_constants = {
            'base_user_path': base_path,
            'pages_folder': pages_folder
        }

    ### Path/folder handling to pass paths across pages
    if 'current_paths' not in st.session_state:
        st.session_state.current_paths = {
        "latest_data_path": data_output_path
        }

    ### Page initialization
    if 'page_initialized' not in st.session_state:
        st.session_state.page_initialized = {
            page: False for page in pages}

    ### Page specific states
    if 'page_states' not in st.session_state:
        st.session_state.page_states = {page: {} for page in pages}

    # ================
    # Initialize page
    # ================

    # clock_loader_40
    # genetics
    # lists
    # hub
    # monitoring
    # analytics
    # grouped_bar_chart
    # search_insights
    # reorder


    home = st.Page("content/1_home.py",
                   title="AbRAT",
                   icon=":material/home:",
                   default=True
                   )
    qc_annotation = st.Page("content/2_ab1_analyzer.py",
                            title="Quality Control & Annotation",
                            icon=":material/troubleshoot:"
                            )
    build_repertoire = st.Page("content/3_bcr_builder.py",
                               title="Compile B-Cell Receptors",
                               icon=":material/genetics:"
                               )
    filter_repertoire = st.Page("content/4_bcr_filtering.py",
                               title="Filter & Split B-Cell Receptors",
                               icon=":material/filter_alt:"
                               )
    clonal_assignment = st.Page("content/5_clonal_assignment.py",
                                title="Clonal Clustering",
                                icon=":material/hub:")
    basic_rep_stats = st.Page("content/6_basic_repertoire_characteristics.py",
                                title="Basic Repertoire Characteristics",
                                icon=":material/grouped_bar_chart:")
    documentation = st.Page("content/7_quickguide.py",
                             title="AbRAT - Quick Guide",
                             icon=":material/menu_book:")

    # Generate navigation menu
    pg = st.navigation(
            {
                "Home": [home],
                "Data Preparation": [qc_annotation, build_repertoire, filter_repertoire],
                "Clonal Assignment": [clonal_assignment],
                "Exploratory & Comparative Analysis": [basic_rep_stats],
                "Documentation":[documentation],
            }
    )

    pg.run()

    st.sidebar.markdown("""_AbRAT v1.0.0 © 2025 Christoph Kreer_""")

if __name__ == "__main__":
    main()
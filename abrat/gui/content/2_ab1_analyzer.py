# gui/content/2_ab1_analyzer.py

import os
import datetime
import streamlit as st
import subprocess
import pandas as pd
from pathlib import Path

from abrat.gui.gui_shared import folder_selectbox, reset_page_initialization

from abrat.core.utils import date_stamp, plot_qc_donut_chart, plot_qc_statistics

# ==============================================
# Passing of global settings from session state
# ==============================================

page_name = os.path.splitext(os.path.basename(__file__))[0]
base_user_path = st.session_state.global_constants['base_user_path']

# Get default navigation start.
if page_name not in st.session_state.current_paths:
    st.session_state.current_paths[page_name] = base_user_path

# default data folder
default_data_folder = os.path.join(base_user_path, "ab1files")
default_output_folder = os.path.join(base_user_path, "output")

other_pages = [page for page in st.session_state.page_initialized.keys() if page!=page_name]

# ========================
# Default settings for QC
# ========================

default_q_cut_off = 16
default_mean_q_cut_off = 28
default_min_length = 240
default_inner_n = 15
default_isotype_determination = True
default_airr_export = False

# =================
# Define functions
# =================

def verify_files(folder_path):
    """
    Checks input files in the specified folder and extracts information from their file names.

    The function searches for all .ab1 files within the folder (recursively) and verifies that each file name
    contains at least 13 segments (separated by underscores). It then compiles a set of unique values for each
    segment position and stores the summary in the session state.

    Parameters:
        folder_path (str or Path): The path to the folder containing the .ab1 files.

    Returns:
        None

    Side Effects:
        Updates st.session_state.ab1file_info with file information and sets st.session_state.ab1files_verified to True.
        Displays success or error messages via Streamlit.
    """
    try:
        files = [p.stem for p in Path(folder_path).rglob("*.ab1")]
        if not files:
            raise FileNotFoundError("No .ab1-files found in " + str(folder_path))

        num_files = len(files)
        splitted_files = []

        for f in files:
            segments = f.split('_')
            if len(segments) < 13:
                raise ValueError(f"File '{f}' has {len(segments)} segments, expected: 13.")
            else:
                splitted_files.append(segments)

        transposed = zip(*splitted_files)
        all_positions = [set(column) for column in transposed]

        # Additional checks can be added here.
        st.session_state.ab1file_info = {
            "Number of files": num_files,
        }
        for i, pos_set in enumerate(all_positions, start=1):
            st.session_state.ab1file_info[f"Pos. {i}"] = ", ".join(sorted(pos_set))

        st.session_state.ab1files_verified = True
        st.success("Files successfully evaluated!")
    except FileNotFoundError as e:
        st.error(f"Error: {e}")
    except ValueError as e:
        st.error(f"Formatting error: {e}")


def reset_verification():
    """
    Resets the session state for file verification.

    Clears the verification status and file information stored in the session state.

    Returns:
        None
    """
    st.session_state.ab1files_verified = False
    st.session_state.ab1file_info = {}


def run_process(cmd):
    """
    Runs the analysis script as a subprocess and returns the process object.

    Parameters:
        cmd (list): The command to execute as a list of arguments.

    Returns:
        subprocess.Popen: The process object for the running subprocess.
    """
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return process


def display_output(process, output_container, error_container):
    """
    Displays output and error messages from a subprocess in real time.

    Reads lines from the subprocess's stdout and stderr, appending them to buffers stored in the session state,
    and updates the provided Streamlit text areas with the latest output.

    Parameters:
        process (subprocess.Popen): The running subprocess.
        output_container: Streamlit container for displaying standard output.
        error_container: Streamlit container for displaying error output.

    Returns:
        None
    """
    if 'output_buffer' not in st.session_state:
        st.session_state.output_buffer = ""
    if 'error_buffer' not in st.session_state:
        st.session_state.error_buffer = ""

    while True:
        # Read output and error from the subprocess.
        output = process.stdout.readline()
        error = process.stderr.readline()

        if output:
            st.session_state.output_buffer += output
            output_container.text_area("Progress:", value=st.session_state.output_buffer, height=100)

        if error:
            st.session_state.error_buffer += error
            error_container.text_area("Errors during processing (stderr):", value=st.session_state.error_buffer, height=100)

        # Break the loop if the process has terminated.
        if output == '' and error == '' and process.poll() is not None:
            break


def run_annotation(s, f):
    """
    Runs the quality check and annotation script based on the provided settings and file information.

    This function creates the output folder (if it doesn't exist), exports the settings to a CSV file, constructs
    the command for running the analysis script, executes it as a subprocess, and displays the output in real time.

    Parameters:
        s (dict): Dictionary of settings.
        f (dict): File information (not directly used here, but typically contains additional file parameters).

    Returns:
        int: The return code of the subprocess.
    """
    with st.spinner("Running..."):
        # Ensure the output folder exists.
        Path(output_folder).mkdir(parents=True, exist_ok=True)

        # Build the command for running the analysis script.
        cmd = [
            "python", "abrat/core/analyze_ab1_files.py",
            "--input", str(data_folder),
            "--output", str(output_folder),
            "--project_name", str(st.session_state.project_name),
            "--q_cut_off", str(q_cut_off),
            "--mean_q_cut_off", str(mean_q_cut_off),
            "--min_length", str(min_length),
            "--inner_n", str(inner_n),
            "--isotype_determination", str(isotype_determination),
            "--airr_export", str(airr_export)
        ]

        # Export settings with a timestamp.
        stamped_settings = {'TIME-STAMP': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
        for k, v in s.items():
            stamped_settings[k] = v
        settings_csv_file = Path(output_folder).joinpath(date_stamp("qc-settings.csv"))
        pd.DataFrame([stamped_settings]).to_csv(settings_csv_file, index=False)

        # Run the command.
        process = run_process(cmd)

        # Create Streamlit elements for output.
        output_text = st.empty()
        error_text = st.empty()

        # Display the output in real time.
        display_output(process, output_text, error_text)

        # Get the return code from the subprocess.
        return_code = process.poll()

        return return_code


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

output_folder = os.path.join(default_output_folder, st.sidebar.text_input(
    "Output folder",
    date_stamp(st.session_state.project_name),
    key="qc_output_path",
    placeholder="path/to/output/folder"
    ))

st.sidebar.divider()
# Set Cutoffs
with st.sidebar.expander("Advanced settings", icon="⚙️"):
    min_length = st.slider("Minimum sequence length",
                           50,
                           300,
                           default_min_length,
                           help="Shorter sequences will not pass the quality control, "
                                "_default:"+str(default_min_length)+"_")
    mean_q_cut_off = st.slider("Mean sequence quality score cutoff",
                               1,
                               66,
                               default_mean_q_cut_off,
                               help="Sequences with a lower mean Phred score will not pass the quality control, "
                                    "_default:"+str(default_mean_q_cut_off)+"_")
    q_cut_off = st.slider("Per base quality score cut off",
                          1,
                          66,
                          default_q_cut_off,
                          help="Base calls below this Phred score will be considered uncertain, "
                               "_default:"+str(default_q_cut_off)+"_")
    inner_n = st.slider("Maximum number of uncertain bases",
                        0,
                        100,
                        default_inner_n,
                        help="Sequences with more uncertain bases will **not pass** the quality control, "
                             "_default:"+str(default_inner_n)+"_")

    isotype_determination = st.toggle("Determine chain isotype", value = default_isotype_determination, key="isotype_selection")
    airr_export = st.toggle("Export IgBLAST results in airr format", value=default_airr_export,
                                      key="airr_export")
# Summary of settings as dict
settings= {
    "Project name": st.session_state.project_name,
    "Input folder": data_folder,
    "Output folder": output_folder,
    "Base call cutoff": q_cut_off,
    "Mean Phred cutoff": mean_q_cut_off,
    "Minimum length": min_length,
    "Inner N": inner_n,
    "Isotype determination": isotype_determination,
    "Export airr results": airr_export
    }

# =========================
# Page content - Main page
# =========================

# Header
st.title("Quality Control and Sequence Annotation")
# Instructions
st.markdown("""
### Module Workflow:
- Copy all ***.ab1 files** in the data folder _(see **Documentation** for *.ab1 file formatting guidelines)_
- Type in a project name, select your input folder and rename the output folder, if required
- Adjust quality filtering settings, if needed _(default settings should work in most scenarios)_
- Review your settings and pre-evaluate your input files
- Review the details of your input files and press "Start Analysis" to start quality checks and sequence annotation
""")

st.markdown("""
  <style>
  div.stButton {text-align:center}
  </style>""", unsafe_allow_html=True)

# Settings visualization
with st.expander("🔍 Show settings"):
    st.dataframe(settings, use_container_width=True)  # Alternative: st.write(settings) oder st.table(settings)

if st.button("Pre-evaluate input files", key='evaluation_button'):
    verify_files(data_folder)
    # reset displayed results
    st.session_state.ab1_results_path = False

if st.session_state.ab1files_verified:
    st.subheader("Data:")
    with st.expander("🔍 Show file information"):
        st.table(st.session_state.ab1file_info)
    # trigger analysis if files were verified
    if st.button("Start Analysis", key="annotation_button"):
        # reset output and error buffer
        st.session_state.output_buffer = ""
        st.session_state.error_buffer = ""

        # run annotation
        annotation_return_code = run_annotation(settings, st.session_state.ab1file_info)
        if annotation_return_code == 0:
            st.success("Analysis successful.")
            st.session_state.ab1_results_path = Path(output_folder) / date_stamp(project_name+\
                                                    "_all-sequences.xlsx")  # might not work if analysis runs at midnight
            # save the latest data path globally
            st.session_state.current_paths['latest_data_path'] = output_folder
            reset_page_initialization(other_pages)
        else:
            st.error(f"Analysis failed with error code {annotation_return_code}.")


if st.session_state.ab1_results_path:
    st.subheader("Latest Results:")
    excel_file_path = st.session_state.ab1_results_path
    st.write(os.path.basename(st.session_state.ab1_results_path))
    if excel_file_path.exists():
        try:
            results_df = pd.read_excel(excel_file_path, index_col=[0], header=[0])
            left_result, right_result = st.columns(2)
            qc_donut_chart = plot_qc_donut_chart(results_df)
            qc_statistics = plot_qc_statistics(results_df)

            #save files:
            donut_chart_path = Path(output_folder).joinpath(date_stamp("qc-statistics.pdf"))
            qc_donut_chart.savefig(donut_chart_path, bbox_inches='tight')
            path_to_summary_plot_file = Path(output_folder).joinpath(date_stamp('qc-fail-statistics.pdf'))
            qc_statistics.savefig(path_to_summary_plot_file, bbox_inches='tight')

            # display results
            with left_result:
                st.caption("Quality check summary:")
                st.pyplot(qc_donut_chart, use_container_width=True)
            with right_result:
                st.caption("Reasons for failed quality checks:")
                st.pyplot(qc_statistics, use_container_width=True)
            st.caption("Exported summary table:")
            st.dataframe(results_df, key="all_sequences_df", use_container_width=True)
        except Exception as e:
            st.error(f"Error while opening results file: {e}")
    else:
        st.error("Results file not found.")


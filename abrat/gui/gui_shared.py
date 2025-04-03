import os
import streamlit as st
import altair as alt
import plotly.io as pio
import matplotlib.pyplot as plt

abrat = '<span style="color: #4682B4; font-style: italic;">AbRAT</span>'

def apply_theme_configurations():
    """
    Applies global theme configurations for Matplotlib, Altair, and Plotly.

    This function sets default styling parameters for Matplotlib (e.g., font, line widths, tick widths),
    defines and enables a custom Altair theme, and creates a custom Plotly template based on the "plotly_white" template.
    The custom Plotly template sets the global font and axis properties, and is then set as the default template.

    Returns:
        None
    """
    # Matplotlib configuration
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.5

    # Altair: define and enable a custom theme
    def my_altair_theme():
        return {
            "config": {
                "font": "DejaVu Sans, sans-serif",
                "title": {
                    "font": "DejaVu Sans, sans-serif",
                    "fontSize": 12,
                    "fontColor": "black"
                },
                "axis": {
                    "labelFont": "DejaVu Sans, sans-serif",
                    "labelFontSize": 12,
                    "titleFont": "DejaVu Sans, sans-serif",
                    "titleFontSize": 12,
                    "domainColor": "black",
                    "domainWidth": 0.5,
                    "tickColor": "black",
                    "tickSize": 5
                },
                "legend": {
                    "labelFont": "DejaVu Sans, sans-serif",
                    "labelFontSize": 12,
                    "titleFont": "DejaVu Sans, sans-serif",
                    "titleFontSize": 12,
                    "fontColor": "black",
                    "symbolSize": 100
                }
            }
        }

    alt.themes.register("my_altair_theme", my_altair_theme)
    alt.themes.enable("my_altair_theme")

    # Plotly: create a custom template based on "plotly_white"
    pio.templates["my_plotly_template"] = pio.templates["plotly_white"]

    # Set global font properties for Plotly
    pio.templates["my_plotly_template"].layout.font.family = "DejaVu Sans, sans-serif"
    pio.templates["my_plotly_template"].layout.font.size = 12
    pio.templates["my_plotly_template"].layout.font.color = "gray"

    # Configure axis lines for Plotly
    pio.templates["my_plotly_template"].layout.xaxis.linewidth = 0.5
    pio.templates["my_plotly_template"].layout.yaxis.linewidth = 0.5

    # Set legend font properties for Plotly
    pio.templates["my_plotly_template"].layout.legend.font.family = "DejaVu Sans, sans-serif"
    pio.templates["my_plotly_template"].layout.legend.font.size = 12
    pio.templates["my_plotly_template"].layout.legend.font.color = "gray"

    pio.templates.default = "my_plotly_template"

### Logic for sidebar select folder

def get_subfolders(path):
    """
    Returns a sorted list of subfolder names from the given path.

    Parameters:
        path (str): The directory path to search for subfolders.

    Returns:
        list: A sorted list of subfolder names. If an error occurs, an error message is displayed and an empty list is returned.
    """
    try:
        folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        return sorted(folders)
    except Exception as e:
        st.error(f"Error while reading subfolders: {e}")
        return []


def update_path(page_name):
    """
    Callback to update the current path based on the user's selection.

    This function retrieves the current path from the session state for the given page, determines the selected folder,
    and updates the session state's current path accordingly. If the selected folder is "..", it navigates up one level.

    Parameters:
        page_name (str): The page name used as key to retrieve and update the current path in the session state.

    Returns:
        None
    """
    select_key = f"{page_name}_selected_folder"
    current_path = os.path.normpath(st.session_state.current_paths[page_name])
    selected_folder = st.session_state[select_key]

    if selected_folder == "..":
        # Navigate up one directory level.
        new_path = os.path.dirname(current_path)
    else:
        # Navigate to the chosen subfolder.
        new_path = os.path.join(current_path, selected_folder)

    # Update the current path in the session state.
    st.session_state.current_paths[page_name] = os.path.normpath(new_path)

def folder_selectbox_on_change_callback(page_name, reset_callback=None):
    """
    Wrapper function for handling folder selectbox changes.

    This callback is triggered when the folder selectbox value changes.
    It calls an optional reset_callback, then updates the current path for the given page.

    Parameters:
        page_name (str): The key for the page in the session state.
        reset_callback (callable, optional): A function to reset additional session state if needed.
    """
    if reset_callback is not None:
        reset_callback()
    update_path(page_name)


def folder_selectbox(page_name, base_user_path, default_folder=None, reset_callback=None, caption="Select folder"):
    """
    Generates a sidebar selectbox for folder selection and updates the session state accordingly.

    This function creates a selectbox in the Streamlit sidebar to allow the user to navigate folders.
    The current path is stored in st.session_state.current_paths under the given page name.
    On the first load, the current path is initialized to the default_folder (if provided) or to the base_user_path.
    The selectbox options include the ".." option for navigating up one directory (if not at the base path)
    along with all subfolders in the current directory.

    Parameters:
        page_name (str): The key for the page to store the current path.
        base_user_path (str or Path): The base directory path from which navigation starts.
        default_folder (str or Path, optional): The default folder for initial load. Defaults to None.
        reset_callback (callable, optional): A function to call when the folder selection changes. Defaults to None.
        caption (str): The label for the selectbox. Default is "Select folder".

    Returns:
        None
    """
    if page_name not in st.session_state.current_paths:
        st.session_state.current_paths[page_name] = base_user_path

    # If the page is loaded for the first time, initialize the current path.
    if not st.session_state.page_initialized[page_name]:
        st.session_state.current_paths[page_name] = os.path.normpath(default_folder) if default_folder else os.path.normpath(base_user_path)
        st.session_state.page_initialized[page_name] = True

    current_path = st.session_state.current_paths[page_name]
    # Define possible options: include ".." for navigating up if not at the base path.
    options = []
    if os.path.abspath(current_path) != os.path.abspath(base_user_path):
        options.append("..")
    # Add all subfolders to the options.
    options.extend(get_subfolders(current_path))

    # Display the selectbox in the sidebar.
    st.sidebar.selectbox(
        label=caption,
        options=options,
        index=None,
        key=f"{page_name}_selected_folder",
        placeholder=os.path.relpath(current_path, base_user_path),
        on_change=lambda: folder_selectbox_on_change_callback(page_name, reset_callback)
    )


def reset_page_initialization(page_names):
    """
    Resets the initialization flags for the specified pages.

    For each page in the provided list, this function sets its 'page_initialized' flag in the session state to False,
    indicating that the page should be reinitialized.

    Parameters:
        page_names (list): List of page keys to reset.

    Returns:
        None
    """
    for page_name in page_names:
        st.session_state.page_initialized[page_name] = False
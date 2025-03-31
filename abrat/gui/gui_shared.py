import os
import streamlit as st
import altair as alt
import plotly.io as pio
import matplotlib.pyplot as plt

def apply_theme_configurations():
    """Wendet die globalen Theme-Konfigurationen für Matplotlib, Altair und Plotly an."""
    # Matplotlib
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.5

    # Altair: definiere und aktiviere das Theme
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

    # Theme for plotly
    # Starte mit einem Basis-Template, hier "plotly_white"
    pio.templates["my_plotly_template"] = pio.templates["plotly_white"]

    # Globale Schriftdefinition
    pio.templates["my_plotly_template"].layout.font.family = "DejaVu Sans, sans-serif"
    pio.templates["my_plotly_template"].layout.font.size = 12
    pio.templates["my_plotly_template"].layout.font.color = "gray"

    # Achsenlinien ähnlich einstellen
    pio.templates["my_plotly_template"].layout.xaxis.linewidth = 0.5
    pio.templates["my_plotly_template"].layout.yaxis.linewidth = 0.5

    # Legende: explizit die gleichen Werte setzen
    pio.templates["my_plotly_template"].layout.legend.font.family = "DejaVu Sans, sans-serif"
    pio.templates["my_plotly_template"].layout.legend.font.size = 12
    pio.templates["my_plotly_template"].layout.legend.font.color = "gray"

    pio.templates.default = "my_plotly_template"

### Logic for sidebar select folder

def get_subfolders(path):
    """get sorted subfolders from a given path"""
    try:
        folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        return sorted(folders)
    except Exception as e:
        st.error(f"Error while reading subfolders: {e}")
        return []

def update_path(page_name):
    """Callback to change the new path according to the selection"""

    select_key = f"{page_name}_selected_folder"
    current_path = os.path.normpath(st.session_state.current_paths[page_name])
    selected_folder = st.session_state[select_key]

    if selected_folder == "..":
        # navigate up if ".." was selected
        new_path = os.path.dirname(current_path)
    else:
        # if subfolder was selected, navigate to chosen folder
        new_path = os.path.join(current_path, selected_folder)

    # update current_path
    st.session_state.current_paths[page_name] = os.path.normpath(new_path)

def folder_selectbox_on_change_callback(page_name, reset_callback=None):
    """Wrapper function for callback on change of the folder selectbox"""
    if reset_callback is not None:
        reset_callback()
    update_path(page_name)

def folder_selectbox(page_name, base_user_path, default_folder=None, reset_callback=None, caption="Select folder"):
    """Generate a sidebar selectbox to select a folder that is saved in the pages session state"""
    if page_name not in st.session_state.current_paths:
        st.session_state.current_paths[page_name] = base_user_path

    # check initialization and override with default folder, if page is loaded for the first time
    if not st.session_state.page_initialized[page_name]:
        st.session_state.current_paths[page_name] = os.path.normpath(default_folder) if default_folder \
            else os.path.normpath(base_user_path)
        st.session_state.page_initialized[page_name] = True

    current_path = st.session_state.current_paths[page_name]
    # define possible options
    options = []
    # as long as we are not in the base path, we will allow ".." to navigate up and add this as an option
    if os.path.abspath(current_path) != os.path.abspath(base_user_path):
        options.append("..")

    # add all subfolders to the options
    options.extend(get_subfolders(current_path))

    # display the selectbox and get selected folder
    st.sidebar.selectbox(
        label=caption,
        options=options,
        index=None,
        key=f"{page_name}_selected_folder",
        placeholder=os.path.relpath(current_path, base_user_path),
        on_change=lambda: folder_selectbox_on_change_callback(page_name, reset_callback)
    )

def reset_page_initialization(page_names):
    for page_name in page_names:
        st.session_state.page_initialized[page_name] = False
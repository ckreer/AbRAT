# gui/content/7_quickguide.py

import os
import streamlit as st
from importlib.resources import files
from pathlib import Path
from abrat.gui.gui_shared import abrat

def remove_sphinx_code_and_images(text):
    """
    Removes lines containing Sphinx-specific directives and image embeddings from the given text.

    This function removes:
      - Lines that start with a code fence for Sphinx directives (e.g., "```{panels}")
      - Lines that start with "::::" (used for panels and columns in Sphinx)
      - Lines that contain image embeddings (lines that include "![")

    Parameters:
        text (str): The input multiline string (e.g. Markdown content).

    Returns:
        str: The modified text with the Sphinx-specific lines removed.
    """
    filtered_lines = []
    for line in text.splitlines():
        if line.lstrip().startswith("```{"):
            continue
        if line.lstrip().startswith("::::"):
            continue
        if "![" in line and "](" in line:
            continue
        filtered_lines.append(line)
    return "\n".join(filtered_lines)

# ==============================================
# Global settings passed from session state
# ==============================================
page_name = os.path.splitext(os.path.basename(__file__))[0]

workflow_path = files('abrat.gui.assets') / 'workflow.png'
workflow_caption = ("The typical workflow of AbRAT consists of (1) Data Preparation, (2) Clonal Assignment, "
                    "and (3) Exploratory & Comparative Analysis.")
#workflow_md = Path("../../../docs/source/quickguide/quickguide_workflow.md").read_text(encoding="utf-8")

docs_path = Path("docs/source/quickguide")

# load chapters
workflow_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_workflow.md").read_text(encoding="utf-8")).replace(
    "### |AbRAT| Workflow", "").replace("|AbRAT|", abrat)
folder_structure_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_folder_structure.md").read_text(encoding="utf-8")).replace("|AbRAT|", abrat)
data_format_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_data_format.md").read_text(encoding="utf-8")).replace("|AbRAT|", abrat)
clonal_assignment_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_clonal_assignment.md").read_text(encoding="utf-8")).replace("|AbRAT|", abrat)
repertoire_char_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_repertoire_characteristics.md").read_text(encoding="utf-8")).replace("|AbRAT|", abrat)
citation_md = remove_sphinx_code_and_images(
    (docs_path / "quickguide_citation_references.md").read_text(encoding="utf-8")).replace("|AbRAT|", abrat)


st.title("Quick Guide")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Workflow",
    "Folder Structure",
    "Data Format",
    "Clonal Assignment",
    "Repertoire Characteristics",
    "Citation & References"
])

with tab1:
    st.markdown("### "+abrat+" Workflow", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 2], vertical_alignment="top")
    with col1:
        st.image(str(workflow_path), caption=workflow_caption)
    with col2:
        st.markdown(workflow_md, unsafe_allow_html=True)

with tab2:
    st.markdown(folder_structure_md, unsafe_allow_html=True)
with tab3:
    st.markdown(data_format_md, unsafe_allow_html=True)
with tab4:
    st.markdown(clonal_assignment_md, unsafe_allow_html=True)
with tab5:
    st.markdown(repertoire_char_md, unsafe_allow_html=True)
with tab6:
    st.markdown(citation_md, unsafe_allow_html=True)


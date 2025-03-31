# content/1_home.py

import streamlit as st

# Begrüßungstext und Workflow-Erklärung
st.title("Welcome to the Antibody Repertoire Analysis Toolkit - :blue[_AbRAT_]")
st.markdown("""
:blue[_AbRAT_] is a suite of tools designed to streamline B cell receptor (BCR) repertoire analysis. Built around the 
single-cell BCR analysis pipeline from [Gieselmann and Kreer et al. (2021)](https://doi.org/10.1038/s41596-021-00554-w), 
this toolkit simplifies sequence data processing and supports a comprehensive BCR analysis workflow.

**Key Features:**

- **Data Preparation:** Transform raw _*.ab1_ sequence data into a well-annotated table of quality-checked, paired BCR 
heavy and light chains.
- **Clonal Assignment:** Execute robust clonal clustering on BCR data.
- **Exploratory & Comparative Analysis:** Examine and compare essential BCR repertoire statistics through interactive visualizations.

Explore :blue[_AbRAT_] to boost your B cell antibody research and accelerate your discoveries!

___

For more information check out the related publication [Kreer et al. (2025)](https://doi.org/10.1038/s41596-021-00554-w)
 and the **documentation**.
""")

st.page_link("content/7_quickguide.py", label="_**Quick Guide**_ 📄")

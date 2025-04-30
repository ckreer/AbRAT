# AbRAT

[![Build Status](https://github.com/ckreer/AbRAT/actions/workflows/ci.yml/badge.svg)](https://github.com/ckreer/AbRAT/actions)
[![Docs](https://readthedocs.org/projects/abrat/badge/?version=latest)](https://abrat.readthedocs.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License](https://img.shields.io/github/license/ckreer/AbRAT)](LICENSE)

## 🌟 What is AbRAT?

AbRAT (Antibody Repertoire Analysis Toolkit) is a Python package and Streamlit app for analyzing single-cell BCR sequencing data.  
It covers data preparation, clonal assignment, and exploratory & comparative analysis.

## 🚀 Quick Start

1) Clone the repo
```bash
git clone https://github.com/ckreer/AbRAT.git
cd AbRAT
```

2) Create and activate your environment
```bash
conda env create -f environment.yml
conda activate abrat_env
```

3) Run the Streamlit GUI
```bash
streamlit run abrat/gui/abrat_app.py
```

Or with Docker:

```bash
docker pull ckreer/abrat:latest
docker run -d -p 8501:8501 ckreer/abrat:latest
```

## 📖 Documentation

For detailed installation and usage guides, visit our full documentation on Read the Docs:

> https://abrat.readthedocs.io/

Included topics:

- **Quick Guide** (Workflow, Folder Structure, Data Format, …)  
- **Installation & Setup** (incl. macOS Terminal instructions)  
- **Clustering Algorithms** (Iterative, Matrix-based, Hierarchical)

## 📜 Citation

If you use this software, please cite:

**Kreer, C. (2025).** AbRAT. DOI: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

## 🛡️ License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](LICENSE) file for details.

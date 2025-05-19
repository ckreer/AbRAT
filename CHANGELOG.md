# Changelog

All notable changes to **AbRAT** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), 
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [Unreleased]
### Fixed
- Ignore unknown amino acids (X) for gravy score calculation
- Save <NA> as pd.nan after clonal clustering in excel sheet
- Convert V Identity to numeric, when importing for repertoire statistics

### Docs
- Update instructions in `README.md` and documentation
- Add `CHANGELOG.md` file

### DevOps
- Add Dockerfile labels for GitHub workflows
- Add GitHub workflows for automated Docker builds

### Removed
- Remove duplicated database files

---

## [1.0.1] - 2025-05-01
### Fixed
- Folder mount issues in Docker when using prebuilt images (missing ab1files and output)

### Improved
- README and documentation consistency across `README.md` and documentation

---

## [1.0.0] - 2025-04-30
### Added
- First stable release of AbRAT
- Streamlit GUI for interactive use
- Repertoire analysis pipeline including clonal clustering
- Integrated analysis dashboard
- Docker support (multi-arch)
- Full documentation with ReadTheDocs and citation metadata (with DOI)

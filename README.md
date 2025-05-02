# AbRAT

[![Docs](https://readthedocs.org/projects/abrat/badge/?version=latest)](https://abrat.readthedocs.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15311638.svg)](https://doi.org/10.5281/zenodo.15311638)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

## 🌟 What is AbRAT?

AbRAT, the Antibody Repertoire Analysis Toolkit, is a containerized Streamlit app and Python package for analyzing single-cell BCR sequencing data
generated with the high-throughput protocol by [Gieselmann and Kreer et al.](https://www.nature.com/articles/s41596-021-00554-w).
It covers data preparation, clonal assignment, and exploratory & comparative analysis.

## 🚀 Recommended Usage: Dockerized App

AbRAT is designed to be used as a containerized application for consistent environment setup and best compatibility.  
I strongly recommend using Docker to run the app.

---

## Option A (Recommended): Use Docker with Local Clone

1. **Clone this repository:**

```bash
git clone https://github.com/ckreer/AbRAT.git
cd AbRAT
```

2. **Start the containerized app:**

```bash
docker-compose up
```

⚠️ Note: On first startup, the container will be built. This may take several minutes depending on your system and internet connection.

Your local `data/` folder will automatically be mounted to the container.
Place your `.ab1` files into:

```
data/userdata/ab1files
```

Results will appear in:

```
data/userdata/output
```
🔓 Once the container is running, open your browser and navigate to [http://localhost:8501](http://localhost:8501) to use the app.

---

## Option B: Use Prebuilt Image from Docker Hub

This option is best if you don’t want to build the image locally.  
**However, you must manually set up mount points on your machine.**

1. **Create local folders for data exchange and database access:**

```bash
mkdir -p ~/abrat_data/userdata/ab1files
mkdir -p ~/abrat_data/userdata/output
```

2. **Copy igblast and blast databases from the repository:**

Copy the full `database` folder from this repository to your `~/abrat_data/` folder with the terminal:

```bash
cp -r ./data/database ~/abrat_data/
```

Or simply use Finder to copy the `databse` folder into your `abrat_data` folder.'

See the [Quick Guide – Folder Structure](https://abrat.readthedocs.io/en/latest/quickguide/quickguide_folder_structure.html) for more details.

3. **Pull and run the image with mounted *abrat_data*-folder:**

```bash
docker pull ckreer/abrat:latest
docker run -d -p 8501:8501 \
  -v ~/abrat_data:/app/data \
  ckreer/abrat:latest
```

4. **Access the app:**  
Go to http://localhost:8501 in your browser.

> ℹ️ Tip: You can also start/stop the container via the Docker Desktop interface.
> Be aware that if you close the container or reload the browser, unsaved session data (e.g., temporary analysis steps) will be lost.

---

## Option C (Fallback): Run the App from Source (Not Recommended)

**⚠️ Only use this method if you can't or don't want to use Docker.**

1. **Clone the repo:**

```bash
git clone https://github.com/ckreer/AbRAT.git
cd AbRAT
```

2. **Create and activate the environment:**

```bash
conda env create -f environment.yml
conda activate abrat_env
```

3. **Install AbRAT as package:**

```bash
pip install -e .
```

4. **Run the Streamlit app:**

```bash
streamlit run abrat/gui/abrat_app.py
```

🧠 Note: On **Apple Silicon (M1/M2/M3/M4)**, `igblast` may not install properly.
If you face errors, launch your terminal with **Rosetta** and rerun:

```bash
arch -x86_64 conda create -n abrat_env python=3.12
arch -x86_64 conda activate abrat_env
arch -x86_64 conda env update -f environment.yml
```

---

## 📖 Documentation

Visit our full documentation on Read the Docs:

> https://abrat.readthedocs.io/

Topics include:

- **Quick Guide** (Workflow, Folder Structure, Data Format, …)  
- **Installation & Setup**  
- **Clonal Clustering Algorithms**

---

## 📜 Citation

If you use this software, please cite:

>Kreer, C. (2025). AbRAT: The Antibody Repertoire Analysis Toolkit (Version 1.0.1) [Computer software]. https://doi.org/10.5281/zenodo.15311638

---

## 🛡️ License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](LICENSE) file for details.
.. AbRAT documentation master file
.. include:: shared.rst

|AbRAT|
=======

Overview
--------
|AbRAT|, the Antibody Repertoire Analysis Toolkit, is a containerized Streamlit app and Python package for analyzing single-cell BCR sequencing data
generated with the high-throughput protocol by `Gieselmann and Kreer et al. <https://www.nature.com/articles/s41596-021-00554-w>`_.
It covers data preparation, clonal assignment, and exploratory & comparative analysis.

This documentation covers installation, setup, usage, API details, and more.
For a quick start, refer to the **Quick Guide** section below.

.. note::
   The Quick Guide, with a brief introduction and usage instructions for AbRAT is also provided
   within the app (accessible after installation).


Installation Quick Guide
^^^^^^^^^^^^^^^^^^^^^^^^
1. Create a destination folder (e.g., ``~/Applications``).
2. Obtain the app:

   - Option 1: Clone with Git.
   - Option 2: Download as ZIP.
   - Option 3: Pull the prebuilt image from Docker Hub.
   - Option 4: Pull the prebuilt image from GitHub Container Registry.

3. Navigate to the project folder (if applicable) and create exchange folders (ab1files and output) if working \
   with images (options 3 & 4)

4. Ensure Docker Desktop is running.
5. Start the app:
   - Use ``docker-compose up`` for Options 1 & 2 (build the docker image locally).
   - Use ``docker run`` for Options 3 & 4 (make sure to add mounted volumes).
6. Access the app via http://localhost:8501.
7. Exchange files by placing your `ab1-files` in the designated folder.
8. To stop and update the app, use ``docker-compose down`` and ``git pull`` (if using Git).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickguide
   clustering_algorithms

Indices and Tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
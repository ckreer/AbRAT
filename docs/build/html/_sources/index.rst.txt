.. AbRAT documentation master file
.. include:: shared.rst

|AbRAT|
=======

Overview
--------
|AbRAT|, the Antibody Repertoire Analysis Toolkit, runs as a containerized Streamlit app.
This documentation covers installation, setup, usage, API details, and more.
For a quick start, refer to the **Quick Guide** section below.

.. note::
   The Quick Guide, with a brief introduction and usage instructions for AbRAT is also provided
   within the app (accessible after installation).


Installation Quick Guide
^^^^^^^^^^^^^^^^^^^^^^^^
1. Create a destination folder (e.g., ``~/Applications``).
2. Obtain the app:
   - Option A: Clone with Git.
   - Option B: Download as ZIP.
   - Option C: Pull the prebuilt image from GitHub Container Registry.
   - Option D: Pull the prebuilt image from Docker Hub.
3. Navigate to the project folder (if applicable).
4. Ensure Docker Desktop is running.
5. Start the app:
   - Use ``docker-compose up`` for Options A & B.
   - Use ``docker run`` for Options C & D.
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
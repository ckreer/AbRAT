.. AbRAT installation documentation file
.. include:: shared.rst

Installation and Setup Guide
============================

.. note::

   This guide is tailored for **macOS**, but the installation instructions and Docker commands also work on **Linux** and **Windows**,
   provided that folder paths and system-specific commands are adjusted accordingly.

Hello and welcome! This guide will help you install and run my containerized |AbRAT| app on your Mac in a designated folder—even if you have little experience with GitHub, Docker, or the macOS Terminal.

Prerequisites
^^^^^^^^^^^^^

1. **Install Docker Desktop**

   - In Finder, open your **Applications** folder and double-click **Docker Desktop** (if already installed).
   - If not installed, download `Docker Desktop for Mac <https://www.docker.com/products/docker-desktop>`_ and double-click the downloaded file to install.
   - Launch Docker Desktop from the Applications folder to ensure it is running in the background.

2. **Git (optional but highly recommended)**

   - If you have Homebrew installed, open **Terminal** (see below) and run:

     .. code-block:: bash

        brew install git

   - Alternatively, download the installer from the `Git website <https://git-scm.com/download/mac>`_ and follow the on-screen instructions.
   - You can also download this repository as a ZIP file via your web browser if you don't want to use git.

Opening the Terminal
^^^^^^^^^^^^^^^^^^^^

If you are new to the macOS Terminal, here are two ways to open it:

   - **Via Finder:**

      1. Open Finder.
      2. Go to **Applications** > **Utilities**.
      3. Double-click **Terminal** to launch it.

   - **Via Spotlight Search:**

      1. Press **⌘ + Spacebar** to open Spotlight.
      2. Type **Terminal** and press **Enter**.

Step-by-Step Guide
^^^^^^^^^^^^^^^^^^

1. **Create a Destination Folder**

   You can use Finder or Terminal to create the folder where AbRAT will be installed.

   - **In Finder:**

     1. Open Finder and navigate to your home directory (click **Go** > **Home** in the menu bar) or
        any other directory, where you want to install |AbRAT|.
     2. Right-click in the folder and select **New Folder**.
     3. Name the folder, e.g. **Applications**.

   - **In Terminal:**

     .. code-block:: bash

        mkdir ~/Applications

     (the '~' symbol represents your home directory)

2. **Obtaining the App**

   You have four options to get started:

     - **Option 1/2:** Clone with Git or download and build the image locally
     - **Option 3/4:** Use a Prebuilt Image from Docker Hub or GitHub Container Registry


   **Option 1: Clone with Git (recommended)**

   - Open **Terminal**.
   - Change directory to your destination folder:

     .. code-block:: bash

        cd ~/Applications

   - Clone the repository:

     .. code-block:: bash

        git clone https://github.com/ckreer/AbRAT.git

   **Option 2: Download as ZIP**

   - In your web browser, visit the repository’s GitHub page.
   - Click the green **Code** button and select **Download ZIP**.
   - Use Finder to move the downloaded ZIP into `~/Applications` and double-click to unzip.

   .. note::

       For Options 1 & 2 skip Options 3 & 4 and continue with building the image locally (Step 3)

  **Option 3: Use Prebuilt Image from Docker Hub**

   .. note::

     This method is good for quickly running AbRAT without building the image locally.
     However, you must **manually create and mount all necessary data folders**,
     including subfolders for data exchange (ab1files, output) and required databases (blastdb, igblastdb).

   - Create local folders for ab1files and output on your machine with Finder or the Terminal:

    .. code-block::

        mkdir -p ~/abrat_data/userdata/ab1files
        mkdir -p ~/abrat_data/userdata/output

   - Copy the `database` folder (including `igblastdb` and `blastdb`) into the same location:

    If you cloned or downloaded the repository, you can copy the entire `database` folder with the terminal:

    .. code-block::

        cp -r .data/database ~/abrat_data/

    Or simply use Finder to copy the `databse` folder into your `abrat_data` folder.

   - Ensure Docker Desktop is running.

   - Open **Terminal** and run:

     .. code-block:: bash

        docker pull ckreer/abrat:latest

   - Start the container with mounted *abrat_data*-folder:

     .. code-block:: bash

        docker run -d -p 8501:8501 \
        -v ~/abrat_data:/app/data \
        ckreer/abrat:latest

   **Option 4: Use Prebuilt Image from GitHub Container Registry**

    As Option 3, but you use another source to pull the image:

   - Create local folders as above (Option 3)

   - Ensure Docker Desktop is running.

   - Open **Terminal** and run:

     .. code-block:: bash

        docker pull ghcr.io/ckreer/abrat:latest

   - Start the container with mounted *abrat_data*-folder:

     .. code-block:: bash

        docker run -d -p 8501:8501 \
        -v ~/abrat_data:/app/data \
        ghcr.io/ckreer/abrat:latest

3. **If using Options 1 & 2: Change to the Project Folder**

   If you cloned or downloaded the repository, open **Terminal** and navigate to the project folder:

   .. code-block:: bash

      cd ~/Applications/AbRAT

4. **Start the Application with Docker Compose (for Options 1 & 2)**

   - Ensure Docker Desktop is running.
   - In **Terminal**, run:

     .. code-block:: bash

        docker-compose up

   Wait until the logs indicate that the app is running (look for a message about Streamlit listening on port 8501).

5. **Access the App and Exchange Files**

   - Open your web browser and navigate to `http://localhost:8501` (or the URL shown in the Docker logs).
   - Use Finder to place your *ab1 files* in the `data/userdata/ab1files` folder inside the `AbRAT` directory, or via Terminal:

6. **Stop the Application**

   - In the **Terminal** window running `docker-compose up`, press **Ctrl + C** to stop.
   - Alternatively, in a new **Terminal** window, navigate to the project folder and run:

     .. code-block:: bash

        docker-compose down

7. **Applying Updates (when using Git)**

   - Open **Terminal** and make sure you are in the `AbRAT` folder:

     .. code-block:: bash

        cd ~/Applications/AbRAT

   - Pull the latest changes:

     .. code-block:: bash

        git pull

   - Restart the app:

     .. code-block:: bash

        docker-compose down
        docker-compose up

.. note::
    Instead of using the Terminal to start/shut down your container, it might be
    more convenient for you to use Docker Desktop.

    Note that all changes (except for the data in the mounted exchange folders) are
    lost, when you shut down the container or when you click reload in the browser.

.. warning::
   If you close the container (via Docker Desktop or Terminal), any unsaved data inside the container will be lost.
   Be sure to always use mounted volumes for input/output files to retain your data.

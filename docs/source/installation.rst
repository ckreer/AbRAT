.. AbRAT installation documentation file
.. include:: shared.rst

Installation and Setup Guide
============================
Hello and welcome! This guide will help you install and run our containerized |AbRAT| app on your Mac in a designated folder—even if you have little experience with GitHub, Docker, or the macOS Terminal.

Prerequisites
^^^^^^^^^^^^^

1. **Install Docker Desktop**
   - In Finder, open your **Applications** folder and double-click **Docker Desktop** (if already installed).
   - If not installed, download `Docker Desktop for Mac <https://www.docker.com/products/docker-desktop>`_ and double-click the downloaded file to install.
   - Launch Docker Desktop from the Applications folder to ensure it is running in the background.

2. **Git (optional but recommended for updates)**
   - If you have Homebrew installed, open **Terminal** (see below) and run:

     .. code-block:: bash

        brew install git

   - Alternatively, download the installer from the `Git website <https://git-scm.com/download/mac>`_ and follow the on-screen instructions.
   - You can also download the repository as a ZIP file via your web browser.

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
     1. Open Finder and navigate to your home directory (click **Go** > **Home** in the menu bar).
     2. Right-click in the folder and select **New Folder**.
     3. Name the folder **Applications**.

   - **In Terminal:**
     .. code-block:: bash

        mkdir ~/Applications

2. **Obtaining the App**

   You have four options to get started:

   **Option A: Clone with Git (recommended)**

   - Open **Terminal**.
   - Change directory to your destination folder:

     .. code-block:: bash

        cd ~/Applications

   - Clone the repository:

     .. code-block:: bash

        git clone https://github.com/ckreer/AbRAT.git

   **Option B: Download as ZIP**

   - In your web browser, visit the repository’s GitHub page.
   - Click the green **Code** button and select **Download ZIP**.
   - Use Finder to move the downloaded ZIP into `~/Applications` and double-click to unzip.

   **Option C: Use Prebuilt Image from GitHub Container Registry**

   - Ensure Docker Desktop is running.
   - Open **Terminal** and run:

     .. code-block:: bash

        docker pull ghcr.io/ckreer/abrat:latest

   - Start the container:

     .. code-block:: bash

        docker run -d -p 8501:8501 ghcr.io/ckreer/abrat:latest

   **Option D: Use Prebuilt Image from Docker Hub**

   - Ensure Docker Desktop is running.
   - Open **Terminal** and run:

     .. code-block:: bash

        docker pull ckreer/abrat:latest

   - Start the container:

     .. code-block:: bash

        docker run -d -p 8501:8501 ckreer/abrat:latest

3. **Change to the Project Folder (if using Options A or B)**

   If you cloned or downloaded the repository, open **Terminal** and navigate to the project folder:

   .. code-block:: bash

      cd ~/Applications/AbRAT

4. **Start the Application with Docker Compose (for Options A & B)**

   - Ensure Docker Desktop is running.
   - In **Terminal**, run:

     .. code-block:: bash

        docker-compose up

   Wait until the logs indicate that the app is running (look for a message about Streamlit listening on port 8501).

5. **Access the App and Exchange Files**

   - Open your web browser and navigate to `http://localhost:8501` (or the URL shown in the Docker logs).
   - Use Finder to place your *ab1 files* in the `data/userdata/ab1files` folder inside the `AbRAT` directory, or via Terminal:

     .. code-block:: bash

        mkdir -p data/userdata/ab1files
        cp /path/to/your/ab1files/*.ab1 data/userdata/ab1files/

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


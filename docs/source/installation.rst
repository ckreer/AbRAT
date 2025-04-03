.. AbRAT installation documentation file

Installation and Setup Guide
============================
Hello and welcome! This guide will help you install and run our containerized |AbRAT| app on your Mac
in a designated folder—even if you have little experience with GitHub and Docker.
You can choose from several options to obtain the app’s code or its prebuilt container image.

Prerequisites
^^^^^^^^^^^^^
1. **Install Docker Desktop**
   - Download `Docker Desktop for Mac <https://www.docker.com/products/docker-desktop>`_ and install it.
   - Launch Docker Desktop to ensure Docker is running in the background.

2. **Git (optional but recommended for updates)**
   - Install Git via Homebrew with ``brew install git`` or download it from the `Git website <https://git-scm.com/download/mac>`_.
   - Alternatively, download the repository as a ZIP file.

Step-by-Step Guide
^^^^^^^^^^^^^^^^^^
1. **Create a Destination Folder**

   Create a folder on your Mac where you want to install the software. For example:

   .. code-block:: bash

      mkdir ~/Applications

2. **Obtaining the App**

   You have four options to get started:

   **Option A: Clone with Git (recommended)**

   - Open Terminal.
   - Change directory to your destination folder:

     .. code-block:: bash

        cd ~/Applications

   - Clone the repository:

     .. code-block:: bash

        git clone https://github.com/ckreer/AbRAT.git

   **Option B: Download as ZIP**

   - Go to the repository’s GitHub page.
   - Click the green **“Code”** button and select **“Download ZIP”**.
   - Save and unzip the file in your destination folder.

   **Option C: Use Prebuilt Image from GitHub Container Registry**

   - Ensure Docker Desktop is running.
   - Open Terminal and run:

     .. code-block:: bash

        docker pull ghcr.io/ckreer/abrat:latest

   - Start the container:

     .. code-block:: bash

        docker run -d -p 8501:8501 ghcr.io/ckreer/abrat:latest

   **Option D: Use Prebuilt Image from Docker Hub**

   - Ensure Docker Desktop is running.
   - Open Terminal and run:

     .. code-block:: bash

        docker pull ckreer/abrat:latest

   - Start the container:

     .. code-block:: bash

        docker run -d -p 8501:8501 ckreer/abrat:latest

3. **Change to the Project Folder (if using Options A or B)**

   If you cloned or downloaded the repository, open Terminal and navigate to the project folder:

   .. code-block:: bash

      cd ~/Applications/AbRAT

4. **Start the Application with Docker Compose (for Options A & B)**

   - Ensure Docker Desktop is running.
   - In Terminal, run:

     .. code-block:: bash

        docker-compose up

   Wait until the logs indicate that the app is running.

5. **Access the App and Exchange Files**

   - Open your browser and navigate to http://localhost:8501 (or the URL shown in the logs).
   - Place your _ab1-files_ in the local ``data/userdata/ab1files`` folder (as defined in the docker-compose.yml).

6. **Stop the Application**

   - Stop it by pressing ``Ctrl + C`` in the Terminal running ``docker-compose up``.
   - Alternatively, open a new Terminal window, navigate to the project folder, and run:

     .. code-block:: bash

        docker-compose down

7. **Applying Updates (when using Git)**

   - Navigate to the project folder:

     .. code-block:: bash

        cd ~/Applications/AbRAT

   - Pull the latest changes:

     .. code-block:: bash

        git pull

   - Restart the app:
     - Stop it using:

       .. code-block:: bash

          docker-compose down

     - Start it again with:

       .. code-block:: bash

          docker-compose up
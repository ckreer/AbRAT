# Installation and Setup Guide

Hello and welcome! This guide will help you install and run our containerized app on your Mac in a designated 
folder—even if you have little experience with GitHub and Docker. You can choose from several options to 
obtain the app’s code or its prebuilt container image.

---

## Prerequisites

1. **Install Docker Desktop**  
   - Download [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop) and install it.  
   - Launch Docker Desktop to ensure Docker is running in the background.

2. **Git (optional but recommended for updates)**  
   - If Git is not installed, you can install it via Homebrew with `brew install git` or download it from the [Git website](https://git-scm.com/download/mac).  
   - If you prefer not to use Git, you can alternatively download the repository as a ZIP file.

---

## Step-by-Step Guide

### 1. Create a Destination Folder

Create a folder on your Mac where you want to install the software.  
Example:
```bash
mkdir ~/Applications
```
### 2. Obtaining the App

You have four options to get started:

**Option A: Clone with Git (recommended)**
1. Open Terminal. 
2. Change directory to your destination folder:

```bash
cd ~/Applications
```
3. Clone the repository:

git clone https://github.com/ckreer/AbRAT.git

This will create a subfolder yourRepository containing all the files.

**Option B: Download as ZIP**
1. Go to the repository’s GitHub page.
2. Click the green **“Code”** button and select **“Download ZIP”**.
3. Save the ZIP file in your destination folder.
4. Unzip the file (double-click the ZIP file in Finder).

**Option C: Use Prebuilt Image from GitHub Container Registry**

If you prefer not to build the image locally, you can pull the prebuilt image directly from the GitHub Container Registry.
1. Ensure Docker Desktop is running.
2. Open Terminal and run:

```bash
docker pull ghcr.io/ckreer/abrat:latest
```

3. Once the image is downloaded, you can start the container with:

```bash
docker run -d -p 8501:8501 ghcr.io/ckreer/abrat:latest
```
Adjust volume mounts or environment variables as needed (consult the docker-compose.yml for additional configuration details).

**Option D: Use Prebuilt Image from Docker Hub**

Alternatively, a prebuilt image is available on Docker Hub.
1. Make sure Docker Desktop is running.
2. Open Terminal and pull the image:

```bash
docker pull ckreer/abrat:latest
```

3. Start the container with:

```bash
docker run -d -p 8501:8501 ckreer/abrat:latest
```

Again, if your container requires specific volume mounts or other settings, adjust the docker run command accordingly.


### 3. Change to the Project Folder (if using Options A or B)

If you chose to clone or download the repository, open Terminal and navigate to the folder containing the source code. 
For example (if cloned via Git):

```bash
cd ~/Applications/AbRAT
```

### 4. Start the Application with Docker Compose (for Options A & B)
- Ensure Docker Desktop is running.
Check that the Docker icon is visible in the menu bar. 
- Start the App:
In Terminal, run:

```bash
docker-compose up
```

This command builds and starts all required containers. 

- Wait for Startup:
Once the logs indicate that the app is running, you can proceed to the next step.

### 5. Access the App and Exchange Files
- Open your browser and navigate to http://localhost:8501 (or the URL shown in the logs) to use the app.
- File Exchange:
Place ab1-files in the local data/userdata/ab1files folder (defined as a volume in the docker-compose.yml).
Any files you add to this folder will be accessible within the app.

### 6. Stop the Application
- Stop via Terminal:
Press Ctrl + C in the Terminal running docker-compose up to stop the process.
- Alternatively:
Open a new Terminal window, navigate to the project folder, and run:

```bash
docker-compose down
```

This command stops and removes all running containers.

### 7. Applying Updates (when using Git)

If there are updates to the repository, follow these steps:
1. Open Terminal and navigate to the project folder:

```bash
cd ~/Applications/AbRAT
```

2. Pull the latest changes:

```bash
git pull
```

3. Restart the app:
- First, stop it using:
```bash
docker-compose down
```
- Then, start it again with:
```bash
docker-compose up
```

### Summary
1. Create a Destination Folder:
E.g., ~/Applications.
2. Download the Repository or Image:
- Option A: Clone with Git
- Option B: Download as ZIP
- Option C: Pull the prebuilt image from GitHub Container Registry
- Option D: Pull the prebuilt image from Docker Hub
3. Navigate to the Project Folder: (if applicable).
4. Launch Docker Desktop:
Ensure Docker is running.
5. Start the App:
- Use docker-compose up for Options A & B.
- Use docker run for Options C & D.
6. Access the App:
Open http://localhost:8501 in your browser.
7. Exchange Files:
Place files in the designated local folder (e.g., data folder) to share with the app.
8. Stop and Update the App:
- Use docker-compose down to stop the app.
- Use git pull for updates (if using Git).

With these steps, you should be able to quickly and easily run our app on your Mac using your preferred method of 
obtaining the code or image. Enjoy!


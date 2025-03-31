# Installation and Setup Guide

Hello and welcome! This guide will help you install and run our containerized app on your Mac in a designated 
folder—even if you have little experience with GitHub and Docker.

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
mkdir ~/MyApp
```
### 2. Download the Repository

You have two options to get the source code and all necessary files (Dockerfile, docker-compose.yml, etc.):

**Option A: Clone with Git (recommended)**
1. Open Terminal. 
2. Change directory to your destination folder:

```bash
cd ~/MyApp
```
3. Clone the repository:

git clone https://github.com/yourUsername/yourRepository.git

This will create a subfolder yourRepository containing all the files.

**Option B: Download as ZIP**
1. Go to the repository’s GitHub page.
2. Click the green **“Code”** button and select **“Download ZIP”**.
3. Save the ZIP file in your destination folder.
4. Unzip the file (double-click the ZIP file in Finder).

### 3. Change to the Project Folder

Open Terminal and navigate to the folder containing the source code.
Example (if cloned via Git):

```bash
cd ~/MyApp/yourRepository
```

### 4. Start the Application with Docker Compose
- Ensure Docker Desktop is running.
Check that the Docker icon is visible in the menu bar. 
- Start the App:
In Terminal, run:

```bash
docker-compose up
```

This command builds and starts all required containers. Your docker-compose.yml file is likely configured with a volume for sharing local files (e.g., ./data:/app/data).

- Wait until the containers are running.
Once the logs indicate that the app is running, you should be able to access it, e.g., at http://localhost:8501.

### 5. Access the App and Exchange Files
	•	Open your browser and navigate to http://localhost:8501 (or the URL shown in the logs) to use the app.
	•	File Exchange:
Place files in the local folder defined as a volume in the docker-compose.yml (e.g., the data folder).
Any files you add to this folder will be accessible within the app.

### 6. Stop the Application
	•	Press Ctrl + C in the Terminal running docker-compose up to stop the process.
	•	Alternatively, open a new Terminal window, navigate to the project folder, and run:

```bash
docker-compose down
```

This command stops and removes all running containers.

### 7. Applying Updates (when using Git)

If there are updates to the repository, follow these steps:
- Open Terminal and navigate to the project folder:

```bash
cd ~/MyApp/yourRepository
```

- Pull the latest changes:

```bash
git pull
```

- Restart the app (first stop it using docker-compose down, then start it again with docker-compose up).

### Summary
	1.	Create a Destination Folder: E.g., ~/MyApp.
	2.	Download the Repository: Either clone with Git or download and unzip the ZIP file.
	3.	Navigate to the Project Folder: Open Terminal and change to the folder containing the source code.
	4.	Launch Docker Desktop: Ensure Docker is running.
	5.	Start the App: Run docker-compose up in Terminal.
	6.	Access the App: Open http://localhost:8501 in your browser.
	7.	Exchange Files: Place files in the designated local folder (e.g., data) to share with the app.
	8.	Stop and Update the App: Use docker-compose down to stop and git pull for updates.

With these steps, you should be able to quickly and easily run our app on your Mac. If you have any questions or need support, feel free to reach out!


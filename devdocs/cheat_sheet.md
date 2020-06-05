To activate your conda environment:
1. Open your Anaconda prompt or terminal
2. Run `conda activate YourEnvironment`, substituting in the name of your environment

To clone the repository from GitHub:
1. Open the GitHub page, click on the green button that says "Clone or Download", and copy the URL in the box
2. Open your Git Bash or terminal, and navigate to the directory you want to store the repository in
3. Run this command: `git clone https://the.url.provided/goes/here`, subbing in the URL you copied

To get any updates from the GitHub repository:
1. Open your Git Bash or terminal
2. Navigate to the directory where the respository is stored on your computer
3. Run this command: `git pull`

To install the GitHub version of cptac on your machine:
1. Open your Anaconda prompt or terminal
2. Activate your development environment
3. Navigate to the cptac directory (the upper one, not the lower one--should contain a file called `setup.py`)
4. Run this command: `pip install .` (don't forget the dot--it's a reference to your current directory, telling pip to build the package based on the setup.py file it finds in the current directory)
5. Alternatively, if you're in a different directory, you can run `pip install /path/to/cptac/directory/with/setup/py/file` (subbing in the proper path, and replacing / with \ if you're on Windows)

To install the PyPI version of cptac on your machine (less up to date):
1. Open your Anaconda prompt or terminal
2. Activate your development environment
3. Run this command: `pip install cptac`

To run Jupyter Notebooks to look at Python Notebooks:
1. Open your Anaconda prompt or terminal
2. Activate your development environment
3. Navigate to the directory with the Python notebooks you want to look at (for example, the cptac/docs folder if you want to look at the use cases and tutorials)
4. Run this command: `jupyter notebook` (singular `notebook`, not `notebooks`)
5. This will start a local Jupyter server in that terminal, and open the Jupyter interface in a tab of your default web browser.

To close Jupyter notebooks:
1. On any open notebooks, click "File", then "Close and halt" in the drop down box
2. On the main Jupyter interface page, click "Quit" in the top right corner. This will close the page, and stop the server process that was running in the terminal you started the Jupyter server from.

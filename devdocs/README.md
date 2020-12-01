# devdocs

These documents give instructions for developers working on the cptac package. They are also useful for others who are working on forks of the repository, or adapting the code for their own uses. **NOTE: Whenever you make any changes to the repository, follow the instructions in [05_HOW_TO_ADD_CODE.md](05_HOW_TO_ADD_CODE.md)**. This will avoid merge conflicts.

**Important:** Whenever you're testing changes to your code, make sure to locally install the package using `pip`, using the following instructions. These instructions will take the local copy of the package that you've been editing and install it in your Anaconda environment's package installation directory. This will make it so that when you've opened a Python prompt or a Jupyter Notebook from that Anaconda environment and then import the package, you'll be importing your edited version of the package. This allows you to test the edits you've made, without having to push them to PyPI. So, to install your locally edited version of the package:
1. Open your Anaconda prompt or terminal
2. Activate your development environment (`conda activate MyEnvironment`, subbing in the name of your environment)
3. Navigate to the cptac directory that contains the `setup.py` file (which is the upper cptac directory, not the lower one). `pip` reads this file to know how to install the package.
4. In that directory, run this command: `pip install .` (don't forget the dot--it's a reference to your current directory, telling pip to build the package based on the `setup.py` file it finds in the current directory)
5. Alternatively, if you're in a different directory, you could run `pip install /path/to/cptac/directory/with/setup/py/file`, subbing in the proper path to the cptac `setup.py` file, and replacing / with \ if you're on Windows. `pip` will follow that path, find the `setup.py` file, and then install the package based off of it.

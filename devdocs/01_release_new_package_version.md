**How to release a new version of the cptac package**



1. Make sure to update:
    1. cptac/version.py
    2. Development status in setup.py
    3. Dependency requirements in setup.py--match whatever versions you have in your current dev environment, especially for pandas. We also list a few dependencies in README.md for the user's convenience.
2. Make sure that if you have updated any datasets, the package can still load the old versions.
3. Use "grep -rn set_trace ." to make sure there are no files with breakpoints in them
4. After committing the last commit for the new version, tag it with the version number: 
    1. Use "git log" to list the commits, and copy the checksum of the commit you want (e.g. 964f16d36dfccde844893cac5b347e7b3d44abbc)
    2. Run this: git tag -a [tag] [checksum], as in 

        git tag -a v0.0 964f16d36dfccde844893cac5b347e7b3d44abbc

    3. By default, git push will not push tags, so after creating a tag, you need to push it to the repository by running "git push origin [tagname]", as in 

        git push origin v0.0

5. After tagging a version, follow these steps to create a release on GitHub: [Creating releases](https://help.github.com/en/articles/creating-releases)
6. After releasing a version on GitHub, follow these steps to release it on PyPI: 
    7. In the parent cptac directory (the one containing setup.py), run this command: "python setup.py bdist_wheel". This will generate several directories, including a "dist" directory, which contains the info we'll upload to PyPI.
    8. Then, in the same cptac directory, run "twine upload dist/*" and enter the lab account username and password as prompted.
        1. You'll need to have installed the 'twine' package to do this. If you don't have it, you can install it through from [Anaconda](https://anaconda.org/conda-forge/twine) or [PyPI](https://pypi.org/project/twine/).
7. After releasing on PyPI, **immediately** update the version.txt file on Box with the latest version number of the package. Do this by creating a new text file on your computer, and then using the "Upload New Version" button on Box to update the file. Make sure to use the "Upload New Version"--if you did something else, it would change the shared URL, which would make it so the package couldn't access the file.

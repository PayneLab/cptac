# How to release a new version of the cptac package

The cptac package is available on two different distribution services, PyPI and Bioconda. If you make a versioned release, you must update it with *both* services. Below are the general instructions for how to do this. Note that PyPI release is fairly simple, but Bioconda relase involves outside parties and some verification. So read up and follow these directions.

1. Make sure to update:
    1. cptac/version.py
    2. Development status in setup.py
    3. CITATION.cff version number
    4. Dependency requirements in setup.py.
2. Make sure that if you have updated any datasets, the package can still load the old versions.
3. Use `grep -rn set_trace .` to make sure there are no files with breakpoints in them
4. Make sure that any finished edits on the dev branch have been merged into the master branch--see [03_HOW_TO_ADD_CODE.md](03_HOW_TO_ADD_CODE.md) for details. (But if there are edits on dev that aren't ready to be released, then make sure to not merge them in.)
5. After committing the last commit for the new version, tag it with the version number: 
    1. Use `git log` to list the commits, and copy the checksum of the commit you want (e.g. 964f16d36dfccde844893cac5b347e7b3d44abbc)
    2. Run this: `git tag -a [tag] [checksum]`, as in 

        `git tag -a v0.0 964f16d36dfccde844893cac5b347e7b3d44abbc`

    3. By default, git push will not push tags, so after creating a tag, you need to push it to the repository by running `git push origin [tagname]`, as in 

        `git push origin v0.0`

6. After tagging a version, follow these steps to create a release on GitHub: [Creating releases](https://help.github.com/en/articles/creating-releases)
7. After releasing a version on GitHub, follow these steps to release it on PyPI: 
    1. In the parent cptac directory (the one containing setup.py), run this command: `python setup.py sdist bdist_wheel`. This will generate several directories, including a "dist" directory, which contains the info we'll upload to PyPI.
    2. Then, in the same cptac directory, run `twine upload dist/*` and enter the lab account username and password as prompted.
        1. You'll need to have installed the `twine` package to do this. If you don't have it, you can install it from [Anaconda](https://anaconda.org/conda-forge/twine) or [PyPI](https://pypi.org/project/twine/).
8. After releasing on PyPI, **immediately** update the version.txt file on Box with the latest version number of the package. Do this by creating a new text file on your computer, and then using the "Upload New Version" button on Box to update the file. Make sure to use the "Upload New Version"--if you did something else, it would change the shared URL, which would make it so the package couldn't access the file.
9. Once the package has been uploaded to PyPI, it can be used to make a recipe for Bioconda. This process should happen automatically! But if it doesn't, instructions for manually updating the package on Bioconda can be found at https://bioconda.github.io/contributor/. Note: It usually takes about half an hour to a day for it to update, and it can be seen at https://anaconda.org/bioconda/cptac

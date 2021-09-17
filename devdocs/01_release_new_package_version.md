# How to release a new version of the cptac package

The cptac package is available on two different distribution services, PyPI and Bioconda. If you make a versioned release, you must update it with *both* services. Below are the general instructions for how to do this. Note that PyPI release is fairly simple, but Bioconda relase involves outside parties and some verification. So read up and follow these directions.

1. Make sure to update:
    1. cptac/version.py
    2. Development status in setup.py
    3. CITATION.cff version number
    4. Dependency requirements in setup.py. For now, Google Colabs doesn't support Python 3.7, but pandas stopped officially supporting Python 3.6 as of 1.2.0. So far the only issue we've seen though is that if you have xlrd>=2.0.0 but pandas<=1.1.5, you won't be able to read .xlsx files. So, for now we're capping xlrd at 1.2.0, and everything should be fine.
2. Make sure that if you have updated any datasets, the package can still load the old versions.
3. Use `grep -rn set_trace .` to make sure there are no files with breakpoints in them
4. Make sure that any finished edits on the dev branch have been merged into the master branch--see [05_HOW_TO_ADD_CODE.md](05_HOW_TO_ADD_CODE.md) for details. (But if there are edits on dev that aren't ready to be released, then make sure to not merge them in.)
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
9. Once the package has been uploaded to PyPI, it can be used to make a recipe for Bioconda. The official bioconda instructions for this can be found at https://bioconda.github.io/contributor/, so go there if you run into problems with this process.

The process is fairly straightforward when you know what's happening, but there are a lot of moving parts. Here are the steps:

1. All bioconda packages are made from the bioconda github, so the first step is to [find and fork](https://github.com/bioconda/bioconda-recipes/fork) the bioconda repository. Great, now you can add your recipe without ruining everything if it doesn't work at first. We still need the files on your machine though, so clone it by going wherever you keep you github things and running these commands:

    1. `git clone https://github.com/<USERNAME>/bioconda-recipes.git`

        Which will create the bioconda-recipes folder, then:
    2. `cd bioconda-recipes`
    3. `git remote add upstream https://github.com/bioconda/bioconda-recipes.git`

2. The next step is to make your own branch to keep things clean and comply with the bioconda workflow. In this case probably something like "cptac" or "update-cptac"
        
      Make sure our master is up to date with Bioconda
    1. `git checkout master`
    2. `git pull upstream master`
    3. `git push origin master`

        Create and checkout a new branch for our work
    4. `git checkout -b update_my_recipe`

3. Now we are going to make a new recipe based on the package uploaded to PyPI. This would normally be done with conda skeleton like so:
`cd recipes`
`conda skeleton pypi pyaml`

   **However**, there is a tool that works better for cptac called grayskull, it is what has been used in the past. You can install it with:
`conda install -c conda-forge grayskull`

    Then simply run `grayskull pypi cptac` in the bioconda recipes folder

   Everything should be taken care of, but make sure the licensing file is included correctly. When everything in the meta.yaml file is looking correct we can go ahead and request our recipe be added to the main bioconda github.

4. The bioconda instructions say:
"Now head back to GitHub. Go to your fork of the bioconda-recipes repository and find the branch you have just created in the Branch: drop down. You should now see a message saying This branch is 1 commit ahead [...] bioconda:master. To the right of that line you will find a button Pull Request. Click this and follow the instructions to open a new Pull Request."

    After that the changes will be automatically tested. If they fail you can fix whatever broke it, then push your changes and it will run again.

5. Once it passes the tests someone associated with bioconda should review it, and if they don't you can go to the gitter page (link is on the bioconda repo readme) and ask for someone to review and merge it. Then you can celebrate, the latest version of cptac will be available on bioconda after half an hour, and can be seen at https://anaconda.org/bioconda/cptac

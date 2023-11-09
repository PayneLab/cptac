# How to Add Code to `cptac`

This document establishes the proper workflow for adding new code to the cptac package. This could be a use case, a new utils function, a new dataset, or anything else. In order to avoid merge conflicts, **no one should ever directly edit the master branch**. Instead, follow this workflow, which will have you make edits on the develop branch, and then merge the develop branch to master.

You can perform the git-related operations using either GitHub desktop, or git bash. As a help, I have provided the git bash commands needed to perform each step, in parentheses at the end of the step.

If at any point you encounter merge conflicts, then it's probably because someone didn't follow these steps, and there were conflicting edits made to both the master and develop branches. Or, you forgot to follow step 2 above and pull from remote before starting your edits. Or, it just so happened that between the time you pulled from latest and the time you committed your changes, someone else also pulled the same section of code and made and pushed their own edits to it. In any case, don't push the changes, and contact whoever's in charge of the package software for instructions on what to do. If you are that person, read up on resolving git merge conflicts if you aren't familiar with it, and fix the merge conflicts. Ask Dr. Payne or someone else for help if you need it--there will be problems and lost work if you resolve it improperly.

1. Switch to the develop branch of the repository (`git checkout develop`).
2. Pull the latest version of the develop branch. **You should always pull latest from remote before making any edits to anything**. (`git pull`)
3. Create and add your code and changes, whether in a new file or an existing file.
5. If you created any new files, first move them to the exact place you want them to be located in the repository, then add them to be tracked by git. (`git add [PUT_IN_YOUR_FILENAME_HERE]`).
5. Commit your changes (`git commit -a`).
6. Push the commit(s) to remote (`git push`).
7. Switch to the master branch (`git checkout master`).
8. Pull the latest version of the master branch (`git pull`).
9. Merge the develop branch into master (`git merge develop`).
10. Push the changes (`git push`).
11. Checkout the develop branch (`git checkout develop`).
12. Merge the master branch into the develop branch (this is just to be safe--you should get a message that says "Already up to date"). (`git merge master`).
13. Push the develop branch. The develop and master branches should now be even with (equivalent to) each other. (`git push`)

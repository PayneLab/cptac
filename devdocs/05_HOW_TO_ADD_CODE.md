# How to add code to `cptac`

This document establishes the proper workflow for adding new code to the cptac package. This could be a use case, a new utils function, a new dataset, or anything else. In brief, you pull the latest version of the dev branch, make your edits on that branch, push them, merge the dev branch into master, and then merge the master branch into dev to make sure everything is equivalent across branches.

You can do this using either GitHub desktop, or git bash. As a help, I have provided the git bash commands needed to perform each step, in parentheses at the end of the step.

If at any point you encounter merge conflicts, then it's probably because someone didn't follow these steps, and there were conflicting edits made to both the master and dev branches. Or, you forgot to follow step 2 above and pull from remote before starting your edits. In either case, don't push the changes, and contact whoever's in charge of the package software for instructions on what to do. If you are that person, read up on resolving git merge conflicts if you aren't familiar with it, and fix the merge conflicts. Ask Dr. Payne or someone else for help if you need it--there will be problems and lost work if you resolve it improperly.

1. Switch to the dev branch of the repository (`git checkout dev`).
2. Pull the latest version of the dev branch. You should always pull latest from remote before making any edits to any branch. (`git pull`).
3. Navigate to the docs/ directory (not devdocs/).
4. Create your use case, either copying in an existing notebook, or creating a new notebook.
5. Add the new file and commit your changes (`git add [PUT_IN_YOUR_FILENAME_HERE]` then `git commit -a`).
6. Push the commit(s) to remote (`git push`).
7. Switch to the master branch (`git checkout master`).
8. Pull the latest version of the master branch (`git pull`).
9. Merge the dev branch into master (`git merge dev`).
10. Push the changes (`git push`).
11. Checkout the dev branch (`git checkout dev`).
12. Merge the master branch into the dev branch (`git merge master`).
13. Push the changes. The dev and master branches should now be even with (equivalent to) each other. 

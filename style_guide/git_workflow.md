git & GitHub Workflow for Contributing
========================================
> **WARNING**: No further development of this project is expected. This section remains solely for completeness.

The best method for contributing software to the `jwst_magic` project is a workflow that involves making a branch off of the `jwst_magic` repository, developing changes on those branches, and opening pull requests through GitHub.

The first question you will have to figure out is whether you should open a JIRA issue. (Aside: At this time all contributors to `jwst_magic` are local to STScI, we ask that issues be made through STScI's JIRA in order to maintain the project management flow that we have established. If you are unsure of how to create and maintain a JIRA issue, please see the *JIRA Issue Creation* section below.) If you think that this change will be solving a significant problem or add a significant enhancement to the project then it would be advantageous to open an issue ticket here. This will allow both individuals and the team as a whole to keep track of the project and our progress as we go. Any appropriate individuals should be assigned to the issue, and a label(s) should be tagged in JIRA.

Following that, any changes that you want to eventually make to the master branch should be done through the workflow where you create a branch and work on your own branch before submitting those changes to be reviewed through a pull request. Instructions on how to do those things can be found below:

1. Make a local copy of the repository by cloning the `jwst_magic` repository by visiting the repository and clone via SSH or HTTPS (e.g. `git clone git@github.com:spacetelescope/jwst_magic.git` or by using a Git GUI such as SourceTree - our developers use SourceTree so we have included the instructions for SourceTree here, as well as the commands necessary to complete these steps in a terminal).  Unless you explicitly delete your clone, this only has to be done once. Note on SSH vs HTTPS: With HTTPS, you will be prompted for your AD credentials any time you push/pull code. SSH, on the other hand,  will require that you set up an SSH key, though this will only need to be done once, and you will never be prompted for your AD credentials when pushing or pulling code.  

2. Create a branch on the cloned repository to develop software changes on. Branch names should start with the name of the JIRA issue, so that they are linked, with a short description of issue; e.g. `JWSTFGS-375-fix-fgs-image-conversion`. Consistent use of hyphens is encouraged.

    *In a terminal:*

   1. `git branch <branchname>` - you only need to do this when you first create your branch.

   2. `git checkout <branchname>` - you can use this command to switch back and forth between existing branches.

   3. Perform local software changes using the nominal git add/git commit -m cycle.

      a. `git status` - allows you to see which files have changed.

      b. `git add <new or changed files you want to commit>`

      c. `git commit -m 'Explanation of changes you've done with these files'`

  *In SourceTree:*

   1. Click on the Branch icon in the top bar

   2. Give your branch a name and choose if you want to commit to the working copy parent (standard) or if there is specific commit you would like to commit to. If you want to immediately start working in this branch, make sure that the "Checkout new branch" box is checked. Click "Create Branch". Note: The branch that you are working in will be bolded under "Branches" in the menu on the left. Double click on different branches to switch branches (make sure you have committed your work before doing so).

   3. Perform local software changes using the nominal commit cycle

      a. Click on "File Status" at the top of the menu on the left (there there should be number next to this that indicates how many files have been changed) *or* by clicking on "Uncommitted changes" in the Description column showing the workflow in the repository.

      b. Select the files or changes to be committed by selecting the check boxes next to each changed or added file in the box at the bottom of the GUI.

      c. Once selected, click Commit icon in the top menu.

      d. In the dialog that opens, include your commit message (as you would with `git commit -m` in the terminal)

      e. Click "Commit"

3. Once the changes to your branch have been made, be sure to run tests locally using the command line command `pytest jwst_magic` before pushing your code and creating a pull request. This will allow GUI tests to be run, which cannot be run with GHA.


4. Push the branch to the GitHub repository - this will deliver all committed changes to the branch version on the web which makes it accessible to other team members. The following are the commands to do this:

  In a terminal:

   1. `git push origin <branchname>` for your first push, or

   2. `git push <branchname>` will also work after the first push of your branch.

  In SourceTree:

   1. Click the Push icon in the top menu

5. On the `jwst_magic` [GitHub repository](https://github.com/spacetelescope/jwst_magic), create a pull request - there will be a button ("Create Pull Request") for this after you have pushed your changes at the top of the page. Note that if the branch is still under heavy development, you can put WIP: at the beginning of the pull request title to signify that the pull request is still a work in progress (e.g. WIP: Example Pull Request Title). Not until the WIP: tag is explicitly removed will the pull request be deemed 'mergable'.

Assign the pull request a reviewer, selecting a member of the `jwst_magic` team. They will review your pull request and either accept the request and merge, or ask for additional changes.

Iterate with your reviewer(s) on additional changes if necessary. This will involve addressing any comments on your pull request which can be found on [this](https://github.com/spacetelescope/jwst_magic/pulls) webpage. You may end up iterating over steps 4 and 5 several times while working with your reviewer - do not despair.

Once the pull request has been accepted and merged, switch to master and delete your local branch with git branch -d <branchname> (from a terminal) or by right clicking on the branch name in the menu on the left and select "Delete <branchname>" (in SourceTree).

Expected Additional Changes to be Made
---------------------------
When you make a code change to `jwst_magic` you are requested to make sure that accompanying documentation is also up-to-date with that change and that you have included comprehensive tests for these changes. For example, if your change alters the GUIs or the steps that a user will take, please update the [JWST MAGIC User's Guide](../docs/magic_user_guide/README.md) and write a test for that GUI change. Your pull request will not be accepted by the reviewer until these updates have been made.

JIRA Issue Creation
-------------------
All issues for these tools exists in the [JWSTFGS project](https://jira.stsci.edu/projects/JWSTFGS/issues/JWSTFGS-76?filter=allopenissues) in JIRA. To create an issue:

1. Click the light blue "Create" button at the top of the page. In the form that pops up, please make sure to fill out the following cells;
* *Issue Type*
* *Summary* - this should be clear and short
* *Assignee*
* *Description* - this should clearly explain the need for the issue and the current plan for executing the solution. Any figures or other attachments should be added here. Also include any JIRA users that need to be watchers on this issue by tagging them with `@` in the description. The JWSTFGS admin will add them as official watchers.
* *Labels* - any MAGIC-related issues require the `JWSTMAGIC` label. Additional labels for each section of the tool or documentation can, and should be added as needed.

Additional cells should be filled out if that information is known.
2. Click the dark blue "Create" button

To start progress on this issue, go to the issue page (click on the issue name from the project home page or from a board), make sure that the issue has been assigned to you (click "Assign to me" under People on the right-hand side), and click the "Start Progress" button under the issue summary.To change where you are in the workflow, you can now click on the "Workflow" button and choose "Send To Backlog" or "Done". Once this issue has been reviewed and merged in GitHub, you are expected to close the issue on JIRA but selecting the "Done" option under "Workflow".


Attribution
------------
This workflow is adapted from the [`spacetelescope` `jwql` git & GitHub workflow for contributing page](https://github.com/spacetelescope/jwql/wiki/git-&-GitHub-workflow-for-contributing).

git & GitLab Workflow for Contributing
========================================
Note on language: In GitHub, "pull request" is used to indicate a request for code to be merged. In GitLab "merge request" is used. Since this repo exists in GitLab we will use the appropriate language however it is important to note this slight difference if you use both GitHub and GitLab.

The best method for contributing software to the `fgs-commissioning` project is a workflow that involves making a branch off of the `fgs-commissioning` repository, developing changes on those branches, and opening merge requests through GitLab.

The first question you will have to figure out is whether you should open a JIRA issue. (Aside: At this time all contributors to `fgs-commissioning` are local to STScI, we ask that issues be made through STScI's JIRA in order to maintain the project management flow that we have established. If you are unsure of how to create and maintain a JIRA issue, please see the JIRA Issue Creation section below.) If you think that this change will be solving a significant problem or add a significant enhancement to the project then it would be advantageous to open an issue ticket here. This will allow both individuals and the team as a whole to keep track of the project and our progress as we go. Any appropriate individuals should be assigned to the issue, and a label(s) should be tagged in JIRA.

Following that, any changes that you want to eventually make to the master branch should be done through the workflow where you create a branch and work on your own branch before submitting those changes to be reviewed through a merge request. Instructions on how to do those things can be found below:

1. Make a local copy of the repository by cloning the `tools` repository by visiting the repository and clone via SSH or HTTPS (e.g. `git clone https://grit.stsci.edu/wfsc/tools.git` or by using a Git GUI such as SourceTree - our developers use SourceTree so we have included the instructions for SourceTree here, as well as the commands necessary to complete these steps in a terminal).  Note that, unless you explicitly delete your clone, this only has to be done once.

2. Create a branch on the cloned repository to develop software changes on. Branch names should start with the name of the JIRA issue (e.g. so that they are linked) with a short description of issue; e.g. `JWSTFGS-375-fix-fgs-image-conversion`. Consistent use of hyphens is encouraged.
#####In a terminal:
   i. git branch <branchname> - you only need to do this when you first create your branch.
   ii. git checkout <branchname> - you can use this command to switch back and forth between existing branches.
   iii. Perform local software changes using the nominal git add/git commit -m cycle.
        a. `git status` - allows you to see which files have changed.
        b. `git add` <new or changed files you want to commit>
        c. `git commit -m` 'Explanation of changes you've done with these files'

#####In SourceTree:
   i. Click on the Branch icon in the top bar
   ii. Give your branch a name and choose if you want to commit to the working copy parent (standard) or if there is specific commit you would like to commit to. If you want to immediately start working in this branch, make sure that the "Checkout new branch" box is checked. Click "Create Branch". Note: The branch that you are working in will be bolded under "Branches" in the menu on the left. Double click on different branches to switch branches (make sure you have committed your work before doing so).
   iii. Perform local software changes using the nominal commit cycle
        a. Click on "File Status" at the top of the menu on the left (there there should be number next to this that indicates how many files have been changed) *or* by clicking on "Uncommitted changes" in the Description column showing the workflow in the repository.
        b. Select the files or changes to be committed by selecting the check boxes next to each changed or added file in the box at the bottom of the GUI.
        c. Once selected, click Commit icon in the top menu.
        d. In the dialog that opens, include your commit message (as you would with `git commit -m` in the terminal)
        e. Click "Commit"

3. Push the branch to the GitHub repository - this will deliver all committed changes to the branch version on the web which makes it accessible to other team members. The following are the commands to do this:
#####In a terminal:
   i. git push origin <branchname> for your first push, or
   ii. git push <branchname> will also work after the first push of your branch.
#####In SourceTree:
   i. Click the Push icon in the top menu
4. On the `tools` `fgs-commissioning` [GitLab repository](https://grit.stsci.edu/wfsc/tools/tree/master/fgs-commissioning), create a merge request - there will be a button ("Create Merge Request") for this after you have pushed your changes at the top of the page. Note that if the branch is still under heavy development, you can put WIP: at the beginning of the merge request title to signify that the merge request is still a work in progress (i.e. WIP: Example Merge Request Title). Not until the WIP: tag is explicitly removed will the merge request be deemed 'mergable'.

Assign the merge request a reviewer, selecting a member of the `fgs-commissioning` team. They will review your merge request and either accept the request and merge, or ask for additional changes.

Iterate with your reviewer(s) on additional changes if necessary. This will involve addressing any comments on your merge request which can be found on [this](https://grit.stsci.edu/wfsc/tools/merge_requests) webpage. You may end up iterating over steps 4.ii, 4.iii and 5.ii several times while working with your reviewer - do not despair.

Once the merge request has been accepted and merged, switch to master and delete your local branch with git branch -d <branchname> (from a terminal) or by right clicking on the branch name in the menu on the left and select "Delete <branchname>" (in SourceTree).

Expected Additional Changes to be Made
---------------------------
When you make a code change to `fgs-commissioning` you are requested to make sure that accompanying documentation is also up-to-date with that change. For example, if your change alters the GUIs or the steps that a user will take, please update the [JWST MAGIC User's Guide](./documentation/JWST_MaGIC_User_Guide.docx). Your merge request will not be accepted by the reviewer until these updates have been made.

JIRA Issue Creation
-------------------
All issues for these tools exists in the [JWSTFGS project](https://jira.stsci.edu/projects/JWSTFGS/issues/JWSTFGS-76?filter=allopenissues) in JIRA. (If you do not have access to this page, contact kbrooks@stsci.edu) To create an issue:

1. Click the light blue "Create" button at the top of the page. In the form that pops up, please make sure to fill out the following cells;
* Issue Type
* Summary - this should be clear and short
* Assignee
* Description - this should clearly explain the need for the issue and the current plan for executing the solution. Any figures or other attachments should be added here. Also include any JIRA users that need to be watchers on this issue by tagging them with `@` in the description. The JWSTFGS admin will add them as official watchers.
* Labels - any MAGIC-related issues require the `JWSTMAGIC` label. Additional labels for each section of the tool or documentation can, and should be added as needed.

Additional cells should be filled out if that information is known.
2. Click the dark blue "Create" button

To start progress on this issue, go to the issue page (click on the issue name from the project home page or from a board), make sure that the issue has been assigned to you (click "Assign to me" under People on the right-hand side), and click the "Start Progress" button under the issue summary.To change where you are in the workflow, you can now click on the "Workflow" button and choose "Send To Backlog" or "Done". Once this issue has been reviewed and merged in GitLab, you are expected to close the issue on JIRA but selecting the "Done" option under "Workflow".


Attribution
------------
This workflow is adapted from the [`spacetelescope` `jwql` git & GitHub workflow for contributing page](https://github.com/spacetelescope/jwql/wiki/git-&-GitHub-workflow-for-contributing).

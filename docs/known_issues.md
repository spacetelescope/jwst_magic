# MAGIC Package Known Issues

As with all software packages, MAGIC has a list of known issues ranging from issues due to computer setup or other packages, to features in the package itself.

## Downloading and Installing

There have been some recent issues reported when trying to download MAGIC from GitLab.

One of the following might occur:

### `xcrun` error
Comand line input:

    git clone https://grit.stsci.edu/JWST-FGS/jwst-magic.git

#### Error
    xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun

#### Solution
In this case, you will need to download Xcode.


### Permission denied (publickey) error

Please note that this must be solved differently depending on if you are cloning the repository from GitLab or GitHub. Please use the section that matches with your situation.

#### GitLab
Command line input:

    git clone git@grit.stsci.edu:JWST-FGS/jwst-magic.git

##### Error

    Cloning into 'jwst-magic'...

    git@grit.stsci.edu: Permission denied (publickey).

    fatal: Could not read from remote repository.


    Please make sure you have the correct access rights and the repository exists.

##### What this error means

This error is telling you that either A) you don't have SSH keys set up, or B) the SSH keys you set up previously no longer work for some reason.

##### Solution

You can do one of the following:

A. Make sure you have SSH keys working: https://docs.gitlab.com/ee/ssh/#rsa-ssh-keys

B. Use https instead:

    git clone https://grit.stsci.edu/JWST-FGS/jwst-magic.git

Note that this method will cause the following (non-fatal) error:

    Fontconfig warning: ignoring UTF-8: not a valid region tag

    Using backend:  Qt5Agg

    Downloading https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat

    |==========================================| 1.3k/1.3k (100.00%)         0s

    MAGIC version status cannot be checked due to the following issue: Warning: Permanently added the RSA host key for IP address '140.82.114.4' to the list of known hosts.

    FGS COUNTRATE version status cannot be checked due to the following issue: git@github.com: Permission denied (publickey).`

### GitHub

Note: If you are cloning the GitHub repository, it is assumed that you are making changes to the repository since you cannot make changes to the repository from GitLab.

Command line input:

    git clone git@github.com:spacetelescope/jwst_magic.git

##### Error

    Cloning into 'jwst-magic'...

    git@github.com: Permission denied (publickey).

    fatal: Could not read from remote repository.


    Please make sure you have the correct access rights and the repository exists.


##### What this error means

This error is telling you that either A) you don't have SSH keys set up, or B) the SSH keys you set up previously no longer work for some reason.

##### Solution

Because STScI requires 2 factor authentication for GitHub, you must use SSH keys for cloning repositories from GitHub.

To solve this error you will have to make sure you have SSH keys working: https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent



### lxml package error
#### Error

When installing MAGIC for the first time, some users have found that they get an error that traces back to the lxml package. At this moment, we have not been able to determine a permanent fix for this.

#### Solution

In your MAGIC environment, uninstall and reinstall `lxml` using `pip`:

    $ pip uninstall lxml

    $ pip install lxml


## MAGIC GUI Known Issues

### MAGIC running slowly
#### Issue
If MAGIC has been open for awhile, it will start to run very slowly.

#### Workaround
In order to make sure MAGIC is still running fast, exit the GUI and exit out of IPython every 5-10 runs/configs. If you are running the GUI with a specific path to an image or out directory (i.e. when running a lot of tests), I recommend copying those locations to a Note or some other place so that you can quickly copy them back into the GUI once it has been started back up.

### MAGIC Main GUI
### "Unable to configure formatter 'myFormatter'"

#### Issue
This is an issue that crops up every once in a while and the cause is not clear from the error message. Essentially you are getting this message because some wires have gotten crossed in the backend. This might happen if you switch branches without closing IPython or for other, less obvious reasons. 

#### Solution
Quit the MAGIC GUI and close the IPython session by typing 'exit' in the terminal. Restart MAGIC:

    $ ipython
    In[1]: import jwst_magic
    In[2]: jwst_magic.run_tool_GUI()


### GUI Crashing
#### Issue
Some individuals have experienced the GUI crashing whenever a cancel button is pressed.

#### Workaround
This behavior has been limited so no clear solution exists at this time. It is suggested that you quit your IPython session and restart the GUI.


### Loading the image: typing in the input image path
#### Issue
When loading an input image in the main GUI, if you start writing the path before clicking "Open", it will not open to the path you have typed, and if you close the file window, the path will be deleted from the line input.

#### Workaround
First click the Open button, when the dialog box appears, type "/" this will open the "Go To" window and you can type in your destination here.

### Loading the image: input FGS RAW image
#### Issue
For the case of an input image being in the FGS RAW frame (this is only expected in a testing case), MAGIC requires that the guider for the input image is the same as the output image. 

#### Workaround
Make any conversions ahead of time or simulate the image for the guider that you want. 

### MAGIC Version 1.1 ( check version number):
#### Issue
Cannot create files for both guiders under the same "root" provided in the GUI. This is not allowed by the GUI.

#### Workaround
If you need to switch guiders, remove the files in the directory in question or rename it if you anticipate needing to use those files for some reason. You can also run the GUI again using a different "root" by using the manual naming section.

### Loading selections from a file
#### Issue
When loading in a guiding selection file in the segment guiding section we currently get the locations and count rate, but for some CARs, we only really want to load the locations and not the count rates.

#### Workaround
We have two different ways to load previous selections. One will override the count rates, the other will only use the locations that are loaded and the count rates of the image in question are used: A) Loading guiding selections files through the main GUI in the star selector section will load the PSF positions and count rates, B) Loading guiding selections files through the pop-out star selector GUI will load only the PSF positions

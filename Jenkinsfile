// JWST MAGIC Jenkinsfile for automatic testing at ssbjenkins.stsci.edu
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
// Authors:
// --------
// - Lauren Chambers
// - Keira Brooks
//
// Notes:
// ------
// - Copied from spacetelescope/jwql GitHub repo on 5/6/19
// - More info here: https://github.com/spacetelescope/jenkinsfile_ci_examples
// - Syntax defined here: https://github.com/spacetelescope/jenkins_shared_ci_utils
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// scm_checkout() does the following:
//     1. Disables pipeline execution if [skip ci] or [ci skip] is present in the
//          commit message, letting users exclude individual commits from CI
//     2. Clones the Git repository
//     3. Creates a local cache of the repository to avoid commit drift between tasks
//        (i.e. Each stage is guaranteed to receive the same source code regardless of
//          commits taking place after the pipeline has started.)
if (utils.scm_checkout()) return

// Define helpful variables
CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda"
CONDA_INSTALL = "conda install -y -q --file=requirements.txt"

// Establish variables for the matrix
matrix_os = ["linux-stable"] // (Note that Jenkins can only be run with Linux, not MacOSX/Windows)
matrix_python = ["3.5", "3.6", "3.7"]

// Set up the matrix of builds
matrix = []

// Iterate over the above variables to define the build matrix.
for (os in matrix_os) {
    for (python_ver in matrix_python) {
        // Define each build configuration, copying and overriding values as necessary.

        // Create a new build configuration
        bc = new BuildConfig()

        // Define the OS (only "linux-stable" used here)
        bc.nodetype = os

        // Give the build configuration a name. This string becomes the
        // stage header on Jenkins' UI. Keep it short!
        bc.name = "linux-py${python_ver}"

        // (Required) Define what packages to include in the base conda environment.
        // This specification also tells Jenkins to spin up a new conda environment for
        // your build, rather than using the default environment.
        bc.conda_packages = ["python=${python_ver}"]

        // Execute a series of commands to set up the build, including
        // any packages that have to be installed with pip
        bc.build_cmds = [
            // Add astroconda as conda channel
            "conda config --add channels '${CONDA_CHANNEL}' ",
            // Install package requirements for given python version
            "${CONDA_INSTALL} python=${python_ver}",
            // Install fgscountrate package
            "git clone git://github.com/spacetelescope/jwst-fgs-countrate.git",
            "pip install jwst-fgs-countrate/",
            // Install jwst_magic
            "pip install ."
            ]

        // Execute a series of test commands
        bc.test_cmds = [
            // Run pytest
            "pytest --junitxml=results.xml",
            // Add a truly magical command that makes Jenkins work for Python 3.5
            "sed -i 's/file=\"[^\"]*\"//g;s/line=\"[^\"]*\"//g;s/skips=\"[^\"]*\"//g' results.xml"
        ]

        // Add the build to the matrix
        matrix += bc
    }
}

// Submit the build configurations and execute them in parallel
utils.run(matrix)

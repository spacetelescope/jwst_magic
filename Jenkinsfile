// Config file for automatic testing at ssbjenkins.stsci.edu

// - Copied from spacetelescope/jwql GitHub repo on 5/6/19
// - More info here: https://github.com/spacetelescope/jenkinsfile_ci_examples
// - Syntax defined here: https://github.com/spacetelescope/jenkins_shared_ci_utils

// Obtain files from source control system.
if (utils.scm_checkout()) return

CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda"
// CONDA_CREATE = "conda create -y -q -n magic --file=requirements.txt"
CONDA_INSTALL = "conda install -y -q --file=requirements.txt"

// Establish variables for the matrix
matrix_os = ["linux-stable"]
matrix_python = ["3.5", "3.6", "3.7"]
matrix = []

for (os in matrix_os) {
    for (python_ver in matrix_python) {
        // Define each build configuration, copying and overriding values as necessary.
        bc = new BuildConfig()
        bc.nodetype = "linux-stable"
        bc.name = "debug-linux-py${python_ver}"
        bc.conda_packages = ["python=${python_ver}"]
        bc.build_cmds = ["conda info --envs",
                         "conda list",
                         "conda config --add channels '${CONDA_CHANNEL}' ",
                         "${CONDA_INSTALL} python=${python_ver}",
                         "conda list",
                         "pip install .",
                         "conda list"]
        bc.test_cmds = ["pytest --junitxml=result.xml",
                        "sed -i 's/file=\"[^\"]*\"//g;s/line=\"[^\"]*\"//g;s/skips=\"[^\"]*\"//g' results.xml"]

        matrix += bc
    }
}

// Iterate over configurations that define the (distributed) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)

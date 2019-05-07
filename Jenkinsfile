// Config file for automatic testing at ssbjenkins.stsci.edu

// - Copied from spacetelescope/jwql GitHub repo on 5/6/19
// - More info here: https://github.com/spacetelescope/jenkinsfile_ci_examples
// - Syntax defined here: https://github.com/spacetelescope/jenkins_shared_ci_utils

// Obtain files from source control system.
if (utils.scm_checkout()) return

CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda"
CONDA_CREATE = "conda create -y -q -n magic --file=requirements.txt"

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
        bc.build_cmds = ["conda config --add channels '${CONDA_CHANNEL}' ",
                         "${CONDA_CREATE} python=${python_ver}",
                         "with_env -n magic pip install ."]
        bc.test_cmds = ["with_env -n magic pytest --junitxml=result.xml"]

        matrix += bc
    }
}

// Iterate over configurations that define the (distributed) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)

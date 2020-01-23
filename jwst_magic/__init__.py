import os
import pathlib
import pkg_resources
import socket
import subprocess

JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from .masterGUI import run_MasterGui as run_tool_GUI
    from .run_magic import run_all as run_tool

module_path = pkg_resources.resource_filename('jwst_magic', '')
setup_path = os.path.normpath(os.path.join(module_path, '../setup.py'))

try:
    with open(setup_path) as f:
        data = f.readlines()

    for line in data:
        if 'VERSION =' in line:
            __version__ = line.split(' ')[-1].replace("'", "").strip()

except FileNotFoundError:
    print('Could not determine jwst_magic version')
    __version__ = '0.0.0'

# Check is MAGIC version matches online version and print results
if "sogs" not in socket.gethostname():
    pathlib.Path(__file__).parent.absolute()
    path = str(pathlib.Path(__file__).parent.absolute()).split('/jwst_magic')[0] + '/jwst_magic'

    cmd = 'git --git-dir={}.git ls-remote --tags git@github.com:spacetelescope/jwst_magic.git'.format(path)
    p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = p.communicate()

    # Check for error in calling version
    if p[0].decode('utf-8') != '':
        tag_list = p[0].decode('utf-8').strip()
        remote_max_version = tag_list.split('\n')[-1].split('tags/')[-1]
        if __version__ == remote_max_version:
            print("Your MAGIC package is up to date")
        else:
            print("**WARNING**: LOCAL MAGIC VERSION {} IS BEHIND THE CURRENT ONLINE VERSION {}\nPlease "
                  "update MAGIC, e.g. run `git pull` and `pip install .` for both\nthe "
                  "jwst_magic and jwst-fgs-countrate packages".format(
                __version__, remote_max_version))
    else:  # no network access
        print('MAGIC version status cannot be checked due to the following issue: {}'.format(
            p[1].decode('utf-8').strip().split('\n')[0]))
else:
    # The SOGS case cannot be checked due to network access on SOGS machines
    print("Your MAGIC package is up to date")

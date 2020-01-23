import importlib
import os
import pkg_resources
import socket
import subprocess

JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from .masterGUI import run_MasterGui as run_tool_GUI
    from .run_magic import run_all as run_tool

module_path = pkg_resources.resource_filename('jwst_magic', '')
setup_path = os.path.normpath(os.path.join(module_path, '../setup.py'))

# Set version number from setup.py
try:
    with open(setup_path) as f:
        data = f.readlines()

    for line in data:
        if 'VERSION =' in line:
            __version__ = line.split(' ')[-1].replace("'", "").strip()

except FileNotFoundError:
    print('Could not determine jwst_magic version')
    __version__ = '0.0.0'

# Check is MAGIC version and the FGS Countrate version matches online versions
if "sogs" not in socket.gethostname():
    # JWST_MAGIC
    magic_path = module_path.split('/jwst_magic')[0] + '/jwst_magic'
    cmd_magic = 'git --git-dir={}.git ls-remote --tags git@github.com:spacetelescope/jwst_magic.git'.format(magic_path)
    p_magic = subprocess.Popen(cmd_magic, shell=True, executable='/bin/bash',
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p_magic = p_magic.communicate()

    # JWST-FGS-COUNTRATE
    cr_path = importlib.util.find_spec('fgscountrate').origin.split('fgscountrate/')[0]
    cmd_cr = 'git --git-dir={}.git ls-remote --tags git@github.com:spacetelescope/jwst-fgs-countrate.git'.format(cr_path)
    p_cr = subprocess.Popen(cmd_cr, shell=True, executable='/bin/bash',
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p_cr = p_cr.communicate()

    # Check MAGIC version
    if p_magic[0].decode('utf-8') != '':
        tag_list = p_magic[0].decode('utf-8').strip()
        latest_remote_version_magic = tag_list.split('\n')[-1].split('tags/')[-1]
        if __version__ == latest_remote_version_magic:
            print("Your MAGIC package is up to date")
        else:
            print("**WARNING**: LOCAL MAGIC VERSION {} IS BEHIND THE CURRENT ONLINE VERSION {}\nPlease "
                  "update MAGIC, e.g. run `git pull` and `pip install .`".format(__version__,
                                                                                 latest_remote_version_magic))
    else:  # no network access
        print('MAGIC version status cannot be checked due to the following issue: {}'.format(
            p_magic[1].decode('utf-8').strip().split('\n')[0]))

    # Check FGSCountrate version
    cr_version = pkg_resources.get_distribution("fgscountrate").version
    if p_cr[0].decode('utf-8') != '':
        tag_list = p_cr[0].decode('utf-8').strip()
        latest_remote_version_cr = tag_list.split('\n')[-1].split('tags/')[-1]
        if cr_version == latest_remote_version_cr:
            print("Your FGS Countrate package is up to date")
        else:
            print("**WARNING**: LOCAL FGS COUNTRATE VERSION {} IS BEHIND THE CURRENT ONLINE VERSION {}\nPlease "
                  "update FGS Countrate, e.g. run `git pull` and `pip install .`".format(cr_version,
                                                                                         latest_remote_version_cr))
    else:  # no network access
        print('FGS Countrate version status cannot be checked due to the following issue: {}'.format(
            p_cr[1].decode('utf-8').strip().split('\n')[0]))

else:
    # The SOGS case cannot be checked due to network access on SOGS machines
    print("Your MAGIC package is up to date")
    print("Your FGS Countrate package is up to date")

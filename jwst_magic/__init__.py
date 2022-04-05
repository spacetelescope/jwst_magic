import importlib
import os
import socket
import subprocess
import pkg_resources

from .utils.utils import check_reference_files

JENKINS = '/home/developer/workspace/' in os.getcwd()
if not JENKINS:
    from jwst_magic.mainGUI import run_MainGui as run_tool_GUI
    from jwst_magic.run_magic import run_all as run_tool

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
    magic_path = module_path.split('/jwst_magic')[0] + '/jwst_magic'
    cr_path = importlib.util.find_spec('fgscountrate').origin.split('fgscountrate/')[0]
    cr_version = pkg_resources.get_distribution("fgscountrate").version
    ssh_magic = 'git@github.com:spacetelescope/jwst_magic.git'
    ssh_cr = 'git@github.com:spacetelescope/jwst-fgs-countrate.git'

    for name, path, ssh, version in [('MAGIC', magic_path, ssh_magic, __version__),
                                     ('FGS COUNTRATE', cr_path, ssh_cr, cr_version)]:
        cmd = f'git --git-dir={path}.git ls-remote --tags {ssh}'
        p = subprocess.Popen(cmd, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p = p.communicate()

        if p[0].decode('utf-8') != '':
            tag_list = p[0].decode('utf-8').strip()
            latest_remote_version = tag_list.split('\n')[-1].split('tags/')[-1]
            latest_remote_version = latest_remote_version.replace('^{}', '')  # remove ^{} from tag if tag was updated

            if version == latest_remote_version:
                print("Your {} package is up to date".format(name))
            else:
                print("**WARNING**: LOCAL {} VERSION {} IS BEHIND THE CURRENT ONLINE VERSION {}\nPlease "
                      "update {}, e.g. run `git pull` and `pip install .`".format(name, version,
                                                                                  latest_remote_version, name))
        else:  # no network access
            print('{} version status cannot be checked due to the following issue: {}'.format(name, p[1].decode(
                                                                                              'utf-8').strip().split(
                                                                                              '\n')[0]))
else:
    # The SOGS case cannot be checked due to network access on SOGS machines
    print("Your MAGIC package is up to date")
    print("Your FGS Countrate package is up to date")

# Make sure that all of our references are up to date and saved locally
try:
    check_reference_files()
except FileNotFoundError:
    print('Warning: Cannot check for newest reference files. See solutions '
          'https://github.com/spacetelescope/jwst_magic#running-the-tools')

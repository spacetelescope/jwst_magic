import importlib
import os
import socket
import subprocess
import pkg_resources

from crds import CrdsLookupError, CrdsNetworkError, CrdsDownloadError

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

# Make sure that all of our references are up to date and saved locally
if not JENKINS:
    try:
        check_reference_files()
    except (FileNotFoundError, CrdsLookupError, CrdsNetworkError, CrdsDownloadError, ValueError):
        print('Warning: Cannot check for newest reference files. See solutions '
              'https://github.com/spacetelescope/jwst_magic#running-the-tools')

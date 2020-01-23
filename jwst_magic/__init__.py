import os
import pkg_resources

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

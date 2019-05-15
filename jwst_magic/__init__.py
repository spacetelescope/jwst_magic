import os
JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from .masterGUI import run_MasterGui as run_tool_GUI
    from .run_magic import run_all as run_tool
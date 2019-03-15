"""Load the config.ini file that contains FGS step parameters

Authors
-------
    - Keira Brooks

Use
---
    This module can be used as such:
    ::
        from jwst_magic.fsw_file_writer import config
        config_ini = config.load_config_ini(config_file_name)

    Required arguments:
        ``config_file_name`` - Path to the config.ini file

"""

# Standard Library Imports
import configparser
import os


def get_config_ini_path(config_file_name):
    """Define path to config file within jwst_magic/fsw_file_writer

    Parameters
    ----------
    config_file_name : str
        Path to config file

    Returns
    -------
    str
        Path to config file
    """
    # return os.path.join('data', config_file_name)
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), config_file_name)


def load_config_ini(config_file_name):
    """Read and parse the config file.

    Parameters
    ----------
    config_file_name : str
        Path to config file

    Returns
    -------
    config : obj
        configparser object containing parameters read from config.ini
    """
    # Read config file once here.
    config = configparser.ConfigParser()
    config._interpolation = configparser.ExtendedInterpolation()

    filename = get_config_ini_path(config_file_name)
    config.read(filename)

    return config

import configparser
import os


def get_config_ini_path(config_file_name):
    #return os.path.join('data', config_file_name)
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), config_file_name)


def load_config_ini(config_file_name):
    # Read config file once here.
    config = configparser.ConfigParser()
    config._interpolation = configparser.ExtendedInterpolation()

    filename = get_config_ini_path(config_file_name)
    config.read(filename)

    return config

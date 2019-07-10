""" Common utilities for the JWST MaGIC test suite.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    ::
        from jwst_magic.tests import utils
        utils.ensure_dir_exists(dir)
"""

import os
import yaml

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def parametrized_data():
    """Load parametrized data from file.

    Returns
    -------
    test_data : dict
        Dictionary containing parametrized test data
    """
    parametrized_data_file = os.path.join(__location__, 'data', 'parametrized_test_data.yml')
    with open(parametrized_data_file) as f:
        test_data = yaml.safe_load(f.read())

    return test_data
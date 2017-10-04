"""Some simple console message functions with counting features for
errors, warnings, and info.  Also error exception raising and
tracebacks.

Adapted from log.py from the HSI package: https://github.com/STScI-SSB/metrology

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# STDLIB
import logging
import optparse
import os
import pprint
import sys
from functools import wraps

# LOCAL
import utils


DEFAULT_VERBOSITY_LEVEL = 50


class Logger(object):
    """Base class for logging."""
    def __init__(self, name, enable_console=True, level=logging.DEBUG,
                 format='%(asctime)-15s %(levelname)-8s %(message)s'):
        self.name = name
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        self.formatter = logging.Formatter(format)
        self.console = None
        if enable_console:
            self.add_console_handler(level)
        self.errors = 0
        self.warnings = 0
        self.infos = 0
        self.eol_pending = False
        try:
            self.verbose_level =  int(os.environ.get(
                self.name.replace(".","_").upper() + "_VERBOSITY", 0))
        except Exception:
            self.verbose_level = DEFAULT_VERBOSITY_LEVEL

    def format(self, *args, **keys):
        end = keys.get("end", "\n")
        sep = keys.get("sep", " ")
        output = sep.join([str(arg) for arg in args]) + end
        return output

    def eformat(self, *args, **keys):
        keys["end"] = ""
        if self.eol_pending:
            self.write()
        return self.format(*args, **keys)

    def info(self, *args, **keys):
        self.infos += 1
        self.logger.info(self.eformat(*args, **keys))

    def error(self, *args, **keys):
        self.errors += 1
        self.logger.error(self.eformat(*args, **keys))

    def exception(self, e):
        """Full traceback."""
        self.errors += 1
        self.logger.exception(e)

    def warn(self, *args, **keys):
        self.warnings += 1
        self.logger.warn(self.eformat(*args, **keys))

    def debug(self, *args, **keys):
        self.logger.debug(self.eformat(*args, **keys))

    def verbose(self, *args, **keys):
        verbosity = keys.get("verbosity", DEFAULT_VERBOSITY_LEVEL)
        if self.verbose_level < verbosity:
            return
        self.debug(*args, **keys)

    def verbose_warning(self, *args, **keys):
        verbosity = keys.get("verbosity", DEFAULT_VERBOSITY_LEVEL)
        if self.verbose_level < verbosity:
            return
        self.warn(*args, **keys)

    def write(self, *args, **keys):
        """Output a message to stdout, formatting each positional parameter
        as a string.
        """
        output = self.format(*args, **keys)
        self.eol_pending = not output.endswith("\n")
        sys.stderr.write(output)
        sys.stderr.flush()

    def status(self):
        return self.errors, self.warnings, self.infos

    def reset(self):
        self.errors = self.warnings = self.infos = 0

    def set_verbose(self, level=True):
        assert 0 <= level <= 100,  "verbosity level must be in range 0..100"
        old_verbose = self.verbose_level
        if level == True:
            level = DEFAULT_VERBOSITY_LEVEL
        elif level == False:
            level = 0
        self.verbose_level = level
        return old_verbose

    def get_verbose(self):
        return self.verbose_level

    def add_console_handler(self, level=logging.DEBUG):
        if self.console is None:
            self.console = self.add_stream_handler(sys.stdout)

    def remove_console_handler(self):
        if self.console is not None:
            self.console = self.remove_stream_handler(self.console)

    def add_stream_handler(self, filelike, level=logging.DEBUG):
        handler = logging.StreamHandler(filelike)
        handler.setLevel(level)
        handler.setFormatter(self.formatter)
        self.logger.addHandler(handler)
        return handler

    def remove_stream_handler(self, handler):
        self.logger.removeHandler(handler)


class FgsLogger(Logger):
    """FGS logger."""
    pass


THE_LOGGER = FgsLogger("FGS")

info = THE_LOGGER.info
error = THE_LOGGER.error
exception = THE_LOGGER.exception
warning = THE_LOGGER.warn
verbose_warning = THE_LOGGER.verbose_warning
verbose = THE_LOGGER.verbose
debug = THE_LOGGER.debug
status = THE_LOGGER.status
reset = THE_LOGGER.reset
write = THE_LOGGER.write
set_verbose = THE_LOGGER.set_verbose
get_verbose = THE_LOGGER.get_verbose
add_console_handler = THE_LOGGER.add_console_handler
remove_console_handler = THE_LOGGER.remove_console_handler
add_stream_handler = THE_LOGGER.add_stream_handler
remove_stream_handler = THE_LOGGER.remove_stream_handler


def set_log_file(filename, level=logging.DEBUG, mode="w+"):
    """Output log info to `filename`."""
    utils.ensure_dir_exists(os.path.dirname(filename))
    f = open(filename, mode)
    h = THE_LOGGER.add_stream_handler(f, level=level)
    return f, h  # So we can close it gracefully


def close_log_file(file_stream, log_handler):
    """See :func:`set_log_file`."""
    THE_LOGGER.remove_stream_handler(log_handler)
    file_stream.close()


def logtofile(logfile):
    """Decorator to log everything to a given file.
    Log is also displayed on screen.

    Parameters
    ----------
    logfile : str
        Full path to the output log file.

    function
        Function to log.

    args, kwargs
        Arguments to the function.

    Returns
    -------
    result
        Output(s) from the function, if any.

    Examples
    --------
    To save log to 'myfunc.log' for ``myfunc()``::

        from hsi import log

        @log.logtofile('myfunc.log')
        def myfunc():
            log.info('My info message.')
            raise ValueError('My traceback.')

    To run the decorated function above:

    >>> myfunc()

    The file 'myfunc.log' should look like this::

        2014-01-24 12:52:02,688 FGS: INFO     My info message.
        2014-01-24 12:52:02,688 FGS: ERROR    My traceback.
        Traceback (most recent call last):
          File ".../log.py", line 218, in wrapper
            result = function(*args, **kwargs)
          File "<ipython-input-4-95d19561f284>", line 4, in myfunc
            raise ValueError('My traceback.')
        ValueError: My traceback.

    """
    def real_decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            # This is not going into the file, just screen.
            THE_LOGGER.info('Started logging to {0}'.format(logfile))
            # Start logging.
            log_stream, log_handler = set_log_file(logfile)
            # Catch traceback, if any.
            try:
                result = function(*args, **kwargs)
            except Exception as e:
                THE_LOGGER.exception(e)
                result = None
            # End logging.
            close_log_file(log_stream, log_handler)
            # Original function output
            return result
        return wrapper
    return real_decorator


def errors():
    """Return the global count of errors."""
    return THE_LOGGER.errors


# ===========================================================================


class PP(object):
    """A wrapper to defer pretty printing until after it's known a verbose
    message will definitely be output.
    """
    def __init__(self, ppobj):
        self.ppobj = ppobj

    def __str__(self):
        return pprint.pformat(self.ppobj)


class Deferred(object):
    """A wrapper to delay calling a callable until after it's known a verbose
    message will definitely be output.
    """
    def __init__(self, ppobj):
        self.ppobj = ppobj

    def __str__(self):
        return str(self.ppobj())


# ===========================================================================


def handle_standard_options(args, parser=None,
                            usage="usage: %prog [options] <inpaths...>"):
    """Set some standard options on an optparser object,  many
    of which interplay with the log module itself.
    """
    if parser is None:
        parser = optparse.OptionParser(usage)

    parser.add_option("-V", "--verbose", dest="verbose",
                      help="Set verbosity level.",
                      metavar="VERBOSITY", default=0)

    options, args = parser.parse_args(args)

    set_verbose(int(options.verbose))

    return options, args


def standard_run(run_str, options, globals_dict, locals_dict):
    """Use options to step run_str, profile run_str,  or just run it."""
    exec(run_str in globals_dict, locals_dict)


def standard_status():
    """Print out errors, warnings, and infos."""
    errors, warnings, infos = THE_LOGGER.status()
    info(errors, "errors")
    info(warnings, "warnings")
    info(infos, "infos")

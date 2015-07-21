import argparse
import cython
from warnings import warn

from operator import attrgetter as oag
from ._bmftools_helper import defaultConfig
from .HTSUtils import printlog as pl, parseSketchConfig, to_bool


class ArgumentSketcher(object):
    """
    Class made for working equally easily with config files
    and command-line arguments.
    Loads the defaultConfig file, then overrides it with the parsed
    config file, and then overrides that with the command-line arguments.

    __init__:
    :param Namespace args - argparse's Namespace output
    :param cystr configPath - path to config file.

    Note: config style is key|value|typechar, where typechar
    is 'b' for bool, 'f' for float, 's' for string, and 'i' for int.
    Anything after a # character is ignored.
    """
    @cython.locals(key=cystr, oag=object)
    def arbitrate(self, oag=oag):
        for key in self.keys:
            value = getattr(self.args, key)
            if(value is not None):
                self.config[key] = value

    @cython.locals(args=object, config=cystr,
                   key=cystr, value=object)
    def __init__(self, args, configPath):
        # assert isinstance(args, argparse.Namespace)
        self.config = {key: value for key, value in
                       defaultConfig.iteritems()}
        # Load in config by copying each item so as to not give
        # ArgumentSketcher a dict subclass.
        try:
            tmpConfDict = parseSketchConfig(configPath)
        except TypeError:
            pl("No config path provided - skipping.")
            tmpConfDict = {}
        for key, value in tmpConfDict.iteritems():
            self.config[key] = value
        self.args = args
        self.keys = [i for i in dir(args) if i[0] != "_"]
        self.arbitrate()

    def __getitem__(self, item):
        return self.config[item]

    def __setitem__(self, key, value):
        self.config[key] = value

    @cython.locals(key=cystr)
    def __call__(self, key):
        return self.config[key]

    def values(self, *args, **kwargs):
        return self.config.values(*args, **kwargs)

    def itervalues(self, *args, **kwargs):
        return self.config.itervalues(*args, **kwargs)

    def iteritems(self, *args, **kwargs):
        return self.config.iteritems(*args, **kwargs)

    def items(self, *args, **kwargs):
        return self.config.items(*args, **kwargs)

    def getkeys(self, *args, **kwargs):
        return self.config.keys(*args, **kwargs)

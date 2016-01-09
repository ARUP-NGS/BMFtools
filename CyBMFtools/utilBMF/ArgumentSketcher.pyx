import argparse
import cython
from warnings import warn
import logging
from copy import copy

from operator import attrgetter as oag
from ._bmftools_helper import DefaultConfig, GetArgumentType
from .HTSUtils import printlog as pl


@cython.returns(bint)
@cython.locals(input_str=cystr)
def to_bool(cystr input_str):
    return (input_str.lower() == "true")


TypeConversionDict = {"s": str, "i": int, "f": float, "b": to_bool,
                      "o": lambda x: x, 'e': lambda x: eval(x)}


@cython.locals(path=cystr, parsedLines=list)
@cython.returns(dict)
def parseSketchConfig(path):
    """
    Parses in a file into a dictionary of key value pairs.

    Note: config style is key|value|typechar, where typechar
    is 'b' for bool, 'f' for float, 's' for string, and 'i' for int.
    Anything after a # character is ignored.
    """
    # Note that the key is mangled to make the key match up with
    # argparse's name
    return {line[0].replace(" ", "_"): TypeConversionDict[GetArgumentType(
                line[0].replace(' ', '_'))](line[1]) for
            line in [l.strip().split("#")[0].split("|") for l in
                     open(path, "r").readlines()
                     if l[0] != "#"]}


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
    'o' for identity (leave unchanged), 'e' to eval() the string to
    make the type.
    Anything after a # character is ignored.
    """
    def arbitrate(self):
        for key in self.keys:
            value = getattr(self.args, key)
            if(value is not None):
                self.config[key] = value

    @cython.locals(args=object, config=cystr,
                   key=cystr, value=object)
    def __init__(self, args, configPath):
        '''
        :param args: [argparse.Namespace/arg] Returned value from argparse.
        :param configPath: [cystr/arg] Path to config file.
        :return: None, like all __init__'s.
        '''

        # Retrieves the default values from DefaultConfig.
        self.config = {key: value[0] for key, value in
                       DefaultConfig.iteritems()}

        # Load in config by copying each item so as to not give
        # ArgumentSketcher a dict subclass.
        try:
            tmpConfDict = parseSketchConfig(configPath)
        except TypeError:
            pl("No config path provided - skipping.", level=logging.DEBUG)
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

# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import logging
import operator
from operator import add as oadd
from operator import sub as osub
from operator import div as odiv
from operator import mul as omul
import math
from math import pow as mpow
from math import sqrt as msqrt
import cython
import numpy as np
from scipy.stats import binom
from scipy.misc import comb
from utilBMF.HTSUtils import printlog as pl, ThisIsMadness

cimport numpy as np

"""
Contains tools relating to calling Copy-Number Alterations.
"""

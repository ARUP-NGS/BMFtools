import numpy as np
import operator
import os
import os.path
import pysam
import shlex
import subprocess
import sys
from Cython.Build import cythonize
# from setuptools import setup
from distutils.core import setup

#  BMFTools tries to install two binaries into your bin folder.
#  If you want to change the install directory, edit this variable
#  or copy it manually.
installDir = "/mounts/bin"

#  Find the ideal -march argument for the system.
try:
    print("Retrieving optimal -march flag.")
    marchValue = subprocess.Popen("gcc -c -Q -march=native --help=target | "
                                  "grep 'march' | awk '{print $NF}'",
                                  shell=True,
                                  stdout=subprocess.PIPE).stdout.read().strip()
    marchFlag = "-march=%s" % marchValue
    if(os.path.isfile("help-dummy.o")):
        os.remove("help-dummy.o")
except ImportError:
    print("Error retrieving optimal -march flag. Give up!")
    marchFlag = ""
compilerList = ["-O3", "-pipe", marchFlag, "-funroll-loops", "-floop-block",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-fivopts", "-ftree-loop-im", "-floop-nest-optimize",
                "-fprefetch-loop-arrays", "-floop-strip-mine", "-flto"]
print("Removing all .c files - this is "
      "important for making sure things get rebuilt.")
subprocess.check_call(shlex.split("find . -name \"*.c\" -exec rm \{\} \\;"))
"""
compilerList = ["-O3", "-pipe", marchFlag, "-funroll-loops", "-floop-block"]
compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse", "-funroll-loops",
                "-floop-strip-mine", "-flto"]
print("Removing all .c files - this is "
      "important for making sure things get rebuilt.")
subprocess.check_call(shlex.split("find . -name \"*.c\" -exec rm \{\} \\;"))

compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse", "-funroll-loops",
                "-floop-strip-mine"]
compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse", "-funroll-loops",
                "-floop-strip-mine"]
compilerList = ["-O3", "-pipe", marchFlag, "-funroll-loops", "-floop-block",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-fivopts", "-ftree-loop-im", "-floop-nest-optimize",
                "-fprefetch-loop-arrays", "-floop-strip-mine", "-flto"]



compilerList = ["marchFlag, "-pipe", "-msse2",
                "-funroll-loops", "-floop-block",
                "-floop-strip-mine", "-floop-nest-optimize", "-ftracer",
                "-fbranch-target-load-optimize2",
                "-ftree-loop-distribution", "-ftree-loop-im", "-fivopts",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-funswitch-loops", "-funsafe-math-optimizations",
                "-fprefetch-loop-arrays", "-fmodulo-sched",
                "-fmodulo-sched-allow-regmoves", "-fgcse",
                "-floop-unroll-and-jam",
                "--mfpmath=sse", "-fomit-frame-pointer"]


compilerList = [marchFlag, "-pipe", "-msse2",
                "-funroll-loops", "-floop-block",
                "-floop-strip-mine", "-floop-nest-optimize", "-ftracer",
                "-fbranch-target-load-optimize2",
                "-ftree-loop-distribution", "-ftree-loop-im", "-fivopts",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-funswitch-loops",
                "-fprefetch-loop-arrays", "-fmodulo-sched",
                "-fmodulo-sched-allow-regmoves", "-fgcse",
                "-floop-unroll-and-jam",
                "-fomit-frame-pointer", "-Ofast"]
"""
ext = cythonize('*/*.pyx') + cythonize("*/*.py") + cythonize("*/*.pxd")
# Insist on -O3 optimization
# If more complex optimizations fail, fall back from line 31 to line 30.
for x in ext:
    x.extra_compile_args += compilerList

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.1.0.0beta',
    'install_requires': ['pysam', 'biopython', 'cytoolz', 'matplotlib',
                         'cython', 'cutadapt', 'lxml', 'scipy', 'entropy',
                         'statsmodels'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster', 'SecC'],
    'ext_modules': ext,
    'include_dirs': [np.get_include()] + pysam.get_include(),
    'scripts': [],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}


setup(**config)

try:
    subprocess.check_call(["cp", "BMFMain/main.py", installDir + "/BMFMain"])
    subprocess.check_call(["cp", "utilBMF/bmftools.py",
                           installDir + "/bmftools"])
except subprocess.CalledProcessError:
    raise ValueError("You don't seem to have permissions to install BMFTools"
                     " executables. You'll have to do that manually.")

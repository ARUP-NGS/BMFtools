import numpy as np
import pysam
import subprocess
import operator
import os
import os.path
import sys

#BMFTools tries to install two binaries into your bin folder.
#If you want to change it, copy 
installDir = "/usr/local/bin"

# Find the ideal -march argument for the system.
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

"""
compilerList = ["-flto", marchFlag, "-pipe", "-msse2",
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
"""
compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse",
                "-funroll-loops", "-floop-unroll-and-jam",
                "-floop-nest-optimize", "-fvariable-expansion-in-unroller"]
# compilerList = ["-O3", "-pipe", marchFlag]


try:
    from setuptools import setup, Extension
except ImportError:
    print("setuptools not available. Trying distutils.")
    from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = cythonize('*/*.pyx') + cythonize("*/*.py")
# Insist on -O3 optimization
# If more complex optimizations fail, fall back from line 31 to line 30.
for x in map(operator.attrgetter("extra_compile_args"), ext):
    # x += ["-Ofast", "-flto", marchFlag, "-pipe",
    x += compilerList


config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.0.7.2',
    'install_requires': ['pysam', 'biopython',
                         'cython', 'cutadapt', 'lxml', 'scipy'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster', 'SECC'],
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
    subprocess.check_call(["cp", "utilBMF/bmftools.py", installDir + "/bmftools"])
except CalledProcessError:
    raise ValueError("You don't seem to have permissions to install BMFTools"
                     " executables. You'll have to do that manually.")

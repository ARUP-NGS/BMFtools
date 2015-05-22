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
compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse", "-funroll-loops"]
"""
compilerList = ["-O3", "-pipe", marchFlag, "-funroll-loops", "-floop-block",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-fivopts", "-ftree-loop-im", "-floop-nest-optimize",
                "-fprefetch-loop-arrays", "-floop-strip-mine", "-flto"]
compilerList = ["-O3", "-pipe", marchFlag, "-funroll-loops", "-floop-block",
                "-fvariable-expansion-in-unroller", "-fsplit-ivs-in-unroller",
                "-fivopts", "-ftree-loop-im", "-floop-nest-optimize",
                "-fprefetch-loop-arrays", "-floop-strip-mine", "-flto"]
compilerList = ["-O3", "-pipe", marchFlag, "-mfpmath=sse", "-funroll-loops",
                "-floop-strip-mine", "-flto"]
print("Removing all .c files - this is "
      "important for making sure things get rebuilt.")
subprocess.check_call(shlex.split("find . -name \"*.c\" -exec rm \{\} \\;"))

"""
ext = cythonize('*/*.pyx') + cythonize("*/*.py") + cythonize("*/*.pxd")
# Insist on -O3 optimization
# If more complex optimizations fail, fall back to -O2
for x in ext:
    x.extra_compile_args += compilerList

install_requires = ['pysam', 'cytoolz', 'matplotlib', 'cython', 'cutadapt',
                    'lxml', 'scipy', 'entropy', 'statsmodels', 'pudb']

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.1.0.0beta',
    'install_requires': install_requires,
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
    print("Could not install the bmftools executables - have you set installDir?")
    sys.exit(1)
for requirement in install_requires:
    if(requirement == "pysam"):
        continue
    try:
        print("Now importing %s" % requirement)
        exec("import %s" % requirement)
    except ImportError:
        raise ImportError("Unable to import requirement %s" % requirement)
sys.exit(0)

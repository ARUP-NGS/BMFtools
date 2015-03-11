import numpy as np
import pysam
import subprocess

try:
    from setuptools import setup, Extension
except ImportError:
    print("setuptools not available. Trying distutils.")
    from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = cythonize('*/*.pyx') + cythonize('*/*.py')

# Find the ideal -march argument for the system.
try:
    print("Retrieving optimal -march flag.")
    marchValue = subprocess.Popen("gcc -c -Q -march=native --help=target | "
                                  "grep 'march' | awk '{print $NF}'",
                                  shell=True,
                                  stdout=subprocess.PIPE).stdout.read().strip()
    marchFlag = "-march=%s" % marchValue
except ImportError:
    print("Error retrieving optimal -march flag. Give up!")
    marchFlag = ""

# Insist on -O3 optimization
for x in ext:
    if(x.extra_compile_args == []):
        x.extra_compile_args = ["-O3", "-flto", marchFlag, "-pipe", "-funroll-loops", "-ftree-vectorize", "-msse2"]
    else:
        x.extra_compile_args += ["-O3", "-flto", marchFlag, "-pipe", "-ftree-vectorize", "-msse2", "funroll-loops"]


config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.0.6.1',
    'install_requires': ['pysam', 'biopython', 'pudb',
                         'cython', 'numconv', 'cutadapt'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster'],
    'ext_modules': ext,
    'include_dirs': [np.get_include()] + pysam.get_include(),
    'scripts': [],
    'name': 'BMFTools',
    'license': 'GPLv3',
    'include': 'README.md',
    'package_data': {'': ['LICENSE']}
}


setup(**config)

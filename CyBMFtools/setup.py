import numpy as np
import operator
import os
import os.path
import pysam
import shlex
import subprocess
import sys
from sys import stderr
from itertools import chain
from Cython.Build import cythonize
# from setuptools import setup
from distutils.core import setup

compilerList = ["-O2", "-pipe", "-march=native", "-mfpmath=sse", "-DSAMTOOLS=1", "-std=c99"]

ext = list(chain.from_iterable(map(cythonize, ['*/*.pyx'])))

# If more complex optimizations fail, fall back to -O2
for x in ext:
    x.extra_link_args += pysam.get_libraries()
    x.define_macros += pysam.get_defines()
    x.include_dirs += ['../include', 'lib', '../']
    x.extra_compile_args += compilerList
    #if x.name ["MawCluster.Math":
    x.sources += ["../include/igamc_cephes.c"]

install_requires = ['pysam>=0.8.3', 'cytoolz', 'matplotlib', 'cython>=0.22',
                    'cutadapt>=1.8']

includes = [np.get_include(), os.path.abspath("src/include"), os.path.abspath("src/include/cephes")] + pysam.get_include()


config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.1.1',
    'install_requires': install_requires,
    'packages': ["BMFMain", "utilBMF", "MawCluster",
                 "SecC", "analyscripts"],
    'ext_modules': ext,
    'include_dirs': includes,
    'scripts': 'scripts/workflow',
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}

config['scripts'] = []


setup(**config)

stderr.write("Installation successful!\n")

sys.exit(0)

# distutils: language = c++
import numpy as np
import operator
import os
import os.path
import pysam
import shlex
import subprocess
import sys
from itertools import chain
from Cython.Build import cythonize
# from setuptools import setup
from distutils.core import setup

marchFlag = "-march=native"

compilerList = ["-O2", "-pipe", marchFlag, "-mfpmath=sse", "-std=c99"]

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
ext = list(chain.from_iterable(map(cythonize, ['*/*.pyx', '*/*.py'])))

# Insist on -O3 optimization
# If more complex optimizations fail, fall back to -O2
for x in ext:
    if(x.name in ['MawCluster.BCFastq', 'utilBMF.MPA']):
        x.sources += ["include/cephes/igam.c", "include/cephes/const.c",
                      "include/cephes/gamma.c", "include/cephes/mtherr.c",
                      "include/cephes/sf_error.c"]
    x.extra_compile_args += compilerList

install_requires = ['pysam>=0.8.2', 'cytoolz', 'matplotlib', 'cython>=0.22',
                    'cutadapt>=1.5', 'lxml', 'scipy', 'entropy', 'statsmodels',
                    're2']

includes = [np.get_include(), os.path.abspath("include"), os.path.abspath("include/cephes")] + pysam.get_include()

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.1.0.2beta',
    'install_requires': install_requires,
    'packages': ["BMFMain", "utilBMF", "MawCluster",
                 "SecC", "analyscripts"],
    'ext_modules': ext,
    'include_dirs': includes,
    'scripts': ['utilBMF/bmftools', 'include/dnbtools'],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']},
}


setup(**config)

print("Installation successful!")

sys.exit(0)

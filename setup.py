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

marchFlag = "-march=native"

compilerList = ["-O2", "-pipe", marchFlag, "-mfpmath=sse", "-std=c99", "-DSAMTOOLS=1",]

"""
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
ext = list(chain.from_iterable(map(cythonize, ['*/*.pyx'])))

# If more complex optimizations fail, fall back to -O2
for x in ext:
    x.extra_link_args += pysam.get_libraries()
    x.define_macros += pysam.get_defines()
    if(x.name in ['MawCluster.BCFastq', 'utilBMF.MPA', 'MawCluster.BCBam', 'MawCluster.Math']):
        x.sources += ["src/include/cephes/igam.c", "src/include/cephes/const.c",
                      "src/include/cephes/gamma.c", "src/include/cephes/mtherr.c",
                      "src/include/cephes/sf_error.c"]
    x.extra_compile_args += compilerList

install_requires = ['pysam>=0.8.3', 'cytoolz', 'matplotlib', 'cython>=0.22',
                    'cutadapt>=1.8']

includes = [np.get_include(), os.path.abspath("src/include"), os.path.abspath("src/include/cephes")] + pysam.get_include()


subprocess.check_call(['make'])

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
    'scripts': ['utilBMF/bmftools', 'src/dmp/fqmarksplit', 'src/sort/lh3sort', 'src/dmp/fqmarksplit_p',
                'src/dmp/crms', 'src/dmp/dmp', 'src/dmp/uthash_dmp', 'src/dmp/bam_rescue',
                'src/dmp/bmfsort', 'src/dmp/crms_db',  'src/dmp/crms_p'],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}


setup(**config)

stderr.write("Installation successful!\n")

sys.exit(0)

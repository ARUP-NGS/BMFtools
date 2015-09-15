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

install_requires = ['pysam>=0.8.2', 'cytoolz', 'matplotlib', 'cython>=0.22',
                    'cutadapt>=1.5', 'lxml', 'scipy', 'entropy', 'statsmodels',
                    're2']

includes = [np.get_include(), os.path.abspath("src/include"), os.path.abspath("src/include/cephes")] + pysam.get_include()

subprocess.check_call('gcc src/dmp/fqmarksplit.c src/dmp/igamc_cephes.c -I src/dmp/ -I src/dmp/include -lm -lz -o src/dmp/fqmarksplit -fopenmp -std=gnu11', shell=True)
subprocess.check_call('cd src/sort;make;cd ../..;', shell=True)
subprocess.check_call('cd src/dmp;gcc fqmarksplit_inline.c -o fqmarksplit_inline -lm -lz -std=gnu11;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc crms.c -I . -o bmf_crms -fopenmp -lm -lz -std=gnu11;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc dmp.c -I. igamc_cephes.c -lm -lz -o dmp -fopenmp  -std=gnu11;cd ../..', shell=True)
#subprocess.check_call('cd src/dmp; gcc -g -Wall -O2 dmp.c igamc_cephes.c isnanl.c -o igamc -fopenmp -lm -std=c99; cd ../..', shell=True)
#gcc  -g -Wall -O2 dmp.c igamc_cephes.c isnanl.c -o omgz -fopenmp -lm -lz -std=gnu99
#gcc -I. -o hash_dmp hash_dmp.c igamc_cephes.c -lz -lm -std=gnu11
#gcc -I. -o dmp dmp.c igamc_cephes.c -lz -lm -std=gnu11

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
    'scripts': ['utilBMF/bmftools', 'src/dmp/fqmarksplit', 'src/dmp/fqmarksplit_inline', 'src/sort/lh3sort',
                'src/dmp/bmf_crms'],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}


setup(**config)

stderr.write("Installation successful!\n")

sys.exit(0)

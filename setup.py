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

subprocess.check_call('gcc src/dmp/fqmarksplit.c src/dmp/igamc_cephes.c -I src/dmp/ -I src/dmp/include -lm -lz -o src/dmp/fqmarksplit -fopenmp -std=gnu11', shell=True)
subprocess.check_call('cd src/sort;make;cd ../..;', shell=True)
subprocess.check_call('cd src/dmp;gcc fqmarksplit_inline.c -o fqmarksplit_inline -lm -lz -DNDEBUG -std=gnu11;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc crms.c -I . -o bmf_crms -fopenmp -lm -lz -DNDEBUG -std=gnu11 ;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc dmp.c -I. igamc_cephes.c -lm -lz -o dmp -fopenmp -DNDEBUG=1 -std=gnu11 -fivopts;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc -I. -o uthash_dmp uthash_dmp.c -lz -lm -std=gnu11 -DNDEBUG;cd ../..', shell=True)
subprocess.check_call('cd src/holloway;gcc htslib/faidx.c htslib/bgzf.c htslib/hfile.c htslib/hfile_net.c htslib/kfunc.c htslib/md5.c htslib/hts.c htslib/kstring.c htslib/sam.c htslib/knetfile.c htslib/cram/mFILE.c htslib/cram/thread_pool.c htslib/cram/pooled_alloc.c htslib/cram/cram_external.c htslib/cram/cram_encode.c htslib/cram/cram_codecs.c htslib/cram/cram_io.c htslib/cram/sam_header.c htslib/cram/files.c htslib/cram/vlen.c htslib/cram/string_alloc.c htslib/cram/cram_decode.c htslib/cram/cram_samtools.c htslib/cram/rANS_static.c htslib/cram/open_trace_file.c htslib/cram/cram_index.c htslib/cram/cram_stats.c htslib/cram/zfio.c sam_opts.c bam_rescue.c ../dmp/igamc_cephes.c -I. -I htslib/htslib -I htslib/ -I ../dmp/include/ -o frmsi -fopenmp -lm -lz -std=gnu11 -DNDEBUG;cd ../..', shell=True)
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
                'src/dmp/bmf_crms', 'src/dmp/dmp', 'src/dmp/uthash_dmp', 'src/dmp/hash_dmp'],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}


setup(**config)

stderr.write("Installation successful!\n")

sys.exit(0)

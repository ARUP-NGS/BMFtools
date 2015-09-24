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

subprocess.check_call('cd src/dmp/include/htslib; make; cd ../../../../', shell=True)
subprocess.check_call('gcc src/dmp/fqmarksplit.c src/dmp/igamc_cephes.c -I src/dmp/ -I src/dmp/include -lm -lz -o src/dmp/fqmarksplit -fopenmp -std=gnu11', shell=True)
subprocess.check_call('cd src/sort;make;cd ../..;', shell=True)
subprocess.check_call('cd src/dmp;gcc fqmarksplit_inline.c -o fqmarksplit_inline -lm -lz -DNDEBUG -std=gnu11;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc crms.c -I . -o bmf_crms -fopenmp -lm -lz -DNDEBUG -std=gnu11 ;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc dmp.c -I. igamc_cephes.c -lm -lz -o dmp -fopenmp -DNDEBUG -std=gnu11 -fivopts;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc -I. -o uthash_dmp uthash_dmp.c -lz -lm -std=gnu11 -DNDEBUG;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc -o vlrcms vlcrms.c -lz -lm -std=gnu11 -DNDEBUG;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc include/htslib/faidx.c include/htslib/bgzf.c include/htslib/hfile.c include/htslib/hfile_net.c include/htslib/kfunc.c include/htslib/md5.c include/htslib/hts.c include/htslib/kstring.c include/htslib/sam.c include/htslib/knetfile.c include/htslib/cram/mFILE.c include/htslib/cram/thread_pool.c include/htslib/cram/pooled_alloc.c include/htslib/cram/cram_external.c include/htslib/cram/cram_encode.c include/htslib/cram/cram_codecs.c include/htslib/cram/cram_io.c include/htslib/cram/sam_header.c include/htslib/cram/files.c include/htslib/cram/vlen.c include/htslib/cram/string_alloc.c include/htslib/cram/cram_decode.c include/htslib/cram/cram_samtools.c include/htslib/cram/rANS_static.c include/htslib/cram/open_trace_file.c include/htslib/cram/cram_index.c include/htslib/cram/cram_stats.c include/htslib/cram/zfio.c sam_opts.c bam_rescue.c ../dmp/igamc_cephes.c -I. -I include/htslib/include/htslib -I include/htslib/ -I ../dmp/include/ -o bam_rescue -fopenmp -lm -lz -std=gnu11 -DNDEBUG;cd ../..', shell=True)
subprocess.check_call('cd src/dmp;gcc include/htslib/faidx.c include/htslib/bgzf.c include/htslib/hfile.c include/htslib/hfile_net.c include/htslib/kfunc.c include/htslib/md5.c include/htslib/hts.c include/htslib/kstring.c include/htslib/sam.c include/htslib/knetfile.c include/htslib/cram/mFILE.c include/htslib/cram/thread_pool.c include/htslib/cram/pooled_alloc.c include/htslib/cram/cram_external.c include/htslib/cram/cram_encode.c include/htslib/cram/cram_codecs.c include/htslib/cram/cram_io.c include/htslib/cram/sam_header.c include/htslib/cram/files.c include/htslib/cram/vlen.c include/htslib/cram/string_alloc.c include/htslib/cram/cram_decode.c include/htslib/cram/cram_samtools.c include/htslib/cram/rANS_static.c include/htslib/cram/open_trace_file.c include/htslib/cram/cram_index.c include/htslib/cram/cram_stats.c include/htslib/cram/zfio.c sam_opts.c bmf_bam_sort.c ../dmp/igamc_cephes.c -I. -I include/htslib/include/htslib -I include/htslib/ -I ../dmp/include/ -o bmf_bam_sort -fopenmp -lm -lz -std=gnu11 -DNDEBUG;cd ../..', shell=True)
#subprocess.check_call('cd src/dmp; gcc -g -Wall -O2 dmp.c igamc_cephes.c isnanl.c -o igamc -fopenmp -lm -std=c99; cd ../..', shell=True)
#gcc  -g -Wall -O2 dmp.c igamc_cephes.c isnanl.c -o omgz -fopenmp -lm -lz -std=gnu99
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
                'src/dmp/bmf_crms', 'src/dmp/dmp', 'src/dmp/uthash_dmp', 'src/dmp/bam_rescue', 'src/dmp/bmf_bam_sort'],
    'name': 'BMFTools',
    'license': 'GNU Affero General Public License, '
               'pending institutional approval',
    'include': 'README.md',
    'package_data': {'': ['README.md']}
}


setup(**config)

stderr.write("Installation successful!\n")

sys.exit(0)

import numpy as np
import pysam

try:
    from setuptools import setup, Extension
except ImportError:
    print("setuptools not available. Trying distutils.")
    from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = cythonize('*/*.pyx') + cythonize('*/*.py')

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.0.5.1',
    'install_requires': ['pysam', 'biopython', 'pudb', 'cython'],
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

import numpy as np

try:
    from setuptools import setup
except ImportError:
    print("setuptools not available. Trying distutils.")
    from distutils.core import setup
from Cython.Build import cythonize

ext = cythonize('*/*.pyx') + cythonize('*/*.py')

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.0.5',
    'install_requires': ['pysam', 'biopython', 'pudb', 'cython'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster'],
    'ext_modules': ext,
    'extra_compile_args': ["-O3"],
    'include_dirs': [np.get_include()],
    'scripts': [],
    'name': 'BMFTools',
    'license': 'GPLv3',
    'include': 'README.md',
    'package_data': {'': ['LICENSE']}
}


setup(**config)

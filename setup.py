try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from Cython.Build import cythonize

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.0.5',
    'install_requires': ['pysam', 'biopython', 'pudb', 'cython'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster'],
    'ext_modules': cythonize('MawCluster/*.py'),
    'extra_compile_args': ["-O3"],
    'scripts': [],
    'name': 'BMFTools',
    'license': 'GPLv3',
#    'long_description': open('dist/README.md').read(),
    'include': 'README.md',
    'package_data': {'': ['LICENSE']}
}


setup(**config)

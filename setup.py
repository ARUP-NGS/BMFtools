try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': '',
    'author': 'Daniel Baker',
    'url': 'https://github.com/ARUP-NGS/BMFTools',
    'author_email': 'daniel.baker@aruplab.com',
    'version': '0.4.0',
    'install_requires': ['pysam', 'biopython'],
    'packages': ['BMFMain', 'utilBMF', 'MawCluster'],
    'scripts': [],
    'name': 'BMFTools',
    'license': 'GPLv3',
#    'long_description': open('dist/README.md').read(),
    'include': 'README.md',
    'package_data': {'': ['LICENSE']}
}


setup(**config)

try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
	'description':'',
	'author': 'Daniel Baker',
	'url': 'https://github.com/ARUP-NGS/BMFTools',
	'author_email': 'daniel.baker@aruplab.com',
	'version':'0.2',
	'install_requires': ['pysam', 'biopython'],
	'packages': ['BMFMain', 'BMFUtils', 'MawCluster'],
	'scripts':[],
	'name': 'BMFTools',
    'license' : 'GPLv3',
    'long_description' : open('dist/README.md').read(),
    'include' : 'README.md',
    'package_data' : {'': ['LICENSE']}
}

setup(**config)

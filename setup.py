from setuptools import setup
from setuptools import find_packages 
import os.path as path

# Get the path to our current directory
path_here = path.abspath(path.dirname(__file__))

# Get the package version from its universal storage location, cptac/version.py
version = {}
version_path = path.join(path_here, "cptac", "version.py")
with open(version_path) as fp:
	exec(fp.read(), version)

# Get the long description from the README file
readme_path = path.join(path_here, "README.md")
with open(readme_path) as readme_file:
    readme_text = readme_file.read()

setup(name='cptac',
	version=version['__version__'],
	description='Python packaging for CPTAC data',
	long_description=readme_text,
	long_description_content_type='text/markdown',
	url='http://github.com/PayneLab/cptac',
	author='Dr. Samuel Payne',
	author_email='sam_payne@byu.edu',
	license='Apache 2.0',
	packages=find_packages(
        where='.',
        include=['cptac*'],  # ["*"] by default
        exclude=[],  # empty by default
    ),
	install_requires=[
		'numpy>=1.16.3',
		'pandas>=2.0.0',
		'requests>=2.21.0',
		'scipy>=1.10.0',
		'openpyxl>=2.6.0',
		'statsmodels>=0.10.0',
		'pyranges>=0.0.111',
        'tqdm>=4.65.0',
	],
	classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
	],
	keywords='bioinformatics cancer proteomics genomics open science open data',
	python_requires='>=3.6',
	zip_safe=False,
	include_package_data=True,
	project_urls={
	   'Documentation': 'https://paynelab.github.io/cptac/'
	},
)

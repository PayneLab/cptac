from setuptools import setup
import os.path

version = {}
version_path = os.path.join("cptac", "version.py") 
with open(version_path) as fp:
	exec(fp.read(), version)

setup(name='cptac',
	version=version['__version__'],
	description='Python packaging for CPTAC data',
	url='http://github.com/PayneLab/cptac',
	author='Dr. Samuel Payne',
	author_email='sam_payne@byu.edu',
	license='Apache 2.0',
	packages=['cptac','cptac.endometrial','cptac.algorithms','cptac.ovarian','cptac.colon'],
	install_requires=[
		'numpy',
		'pandas'
	],
	zip_safe=False,
	include_package_data=True
	)

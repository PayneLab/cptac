from setuptools import setup

version = {}
with open('version.py') as fp:
	exec(fp.read(), version)

setup(name='CPTAC',
	version=version['__version__'],
	description='Python packaging for CPTAC data',
	url='http://github.com/PayneLab/CPTAC',
	author='Dr. Samuel Payne',
	author_email='sam_payne@byu.edu',
	license='MIT',
	packages=['CPTAC'],
	zip_safe=False,
	include_package_data=True
	)

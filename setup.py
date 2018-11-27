from setuptools import setup

exec(open('CPTAC/version.py').read())

setup(name='CPTAC',
	version=__version__,
	description='Python packaging for CPTAC data',
	url='http://github.com/PayneLab/CPTAC',
	author='Dr. Samuel Payne',
	author_email='sam_payne@byu.edu',
	license='MIT',
	packages=['CPTAC'],
	zip_safe=False,
	include_package_data=True
	)

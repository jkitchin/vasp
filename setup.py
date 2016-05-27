from setuptools import setup

setup(name = 'vasp',
      version='0.1',
      description='ase-compliant calculator for Vasp',
      url='http://github.com/jkitchin/vasp',
      author='John Kitchin',
      author_email='jkitchin@andrew.cmu.edu',
      license='GPL',
      platforms=['linux'],
      packages=['vasp'],
      scripts=['bin/runvasp.py', 'bin/vaspsum'],
      test_suite='nose.collector',
      long_description='''Python module for setting up, running and analysing VASP calculations.''',
      install_requires=[
          'nose',
          "numpy",
          "matplotlib",
          "scipy",
          ],)

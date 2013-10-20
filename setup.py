#!/usr/bin/env python
import sys

if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup
else:
    from distutils.core import setup

with open('README.rst') as file:
    long_description = file.read()

#with open('CHANGES') as file:
#    long_description += file.read()


from pyradmc3d import __version__ as version

setup(name='pyradmc3d',
      version=version,
      description='Python-radmc3d',
      long_description=long_description,
      author='Adam Ginsburg',
      author_email='adam.g.ginsburg@gmail.com',
      url='http://github.com/keflavich/pyradmc3d/',
      packages=['pyradmc3d'],
      )

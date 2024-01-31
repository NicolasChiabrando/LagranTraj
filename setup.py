#! /usr/bin/env python

from setuptools import setup


setup(name='lagrantraj',
      install_requires=['numpy','scipy','netCDF4','xarray','pandas'],
      packages=['lagrantraj']
      )

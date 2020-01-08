#!/usr/bin/env python

''' Setup script, only usable via pip. '''

from setuptools import setup

setup(
    name='radical.facts',
    version='0.1.0',
    description='Framework for Assessing Changes To Sea-level (FACTS)',
    url='https://github.com/radical-collaboration/facts.git',
    author='Gregory Garner',
    author_email='gregory.garner@rutgers.edu',
    license='unlicense',
    packages=['facts'],
    zip_safe=False
)

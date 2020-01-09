#!/usr/bin/env python

''' Setup script, only usable via pip. '''

import os
import sys
import glob

from setuptools import setup, find_packages

# ------------------------------------------------------------------------------
# check python version. we need >= 2.7, <3.x
if  sys.hexversion < 0x02070000 or sys.hexversion >= 0x03000000:
    raise RuntimeError('%s requires Python 2.x (2.7 or higher)' % 'facts')


# ------------------------------------------------------------------------------
#
# This copies the contents like examples/ dir under sys.prefix/share/$name
# It needs the MANIFEST.in entries to work.
base = 'share/%s' % 'facts'
df = [('%s/'                      		% base, ['FACTS.py']),
     ('%s/experiments/temp_exp'                 % base, glob.glob('experiments/temp_exp/*.yml')),
     ('%s/experiments/temp_exp/input'           % base, glob.glob('experiments/temp_exp/input/*')),
     ('%s/modules/genmod/directsample'  	% base, glob.glob('modules/genmod/directsample/*')),
     ('%s/modules/gilford/icesheets'     	% base, glob.glob('modules/gilford/icesheets/*')),
     ('%s/modules/kopp14/glaciers'              % base, glob.glob('modules/kopp14/glaciers/*')),
     ('%s/modules/kopp14/icesheets'             % base, glob.glob('modules/kopp14/icesheets/*')),
     ('%s/modules/kopp14/landwaterstorage'      % base, glob.glob('modules/kopp14/landwaterstorage/*')),
     ('%s/modules/kopp14/lib'                   % base, glob.glob('modules/kopp14/lib/*')),
     ('%s/modules/kopp14/oceandynamics'         % base, glob.glob('modules/kopp14/oceandynamics/*')),
     ('%s/modules/kopp14/thermalexpansion'      % base, glob.glob('modules/kopp14/thermalexpansion/*')),
     ('%s/modules/kopp14/verticallandmotion'    % base, glob.glob('modules/kopp14/verticallandmotion/*')),
     ('%s/modules/kopp14SROCC/icesheets'        % base, glob.glob('modules/kopp14SROCC/icesheets/*')),
     ('%s/modules/ssp/landwaterstorage'         % base, glob.glob('modules/ssp/landwaterstorage/*')),
     ('%s/other'                  		% base, glob.glob('other/*')),
]


# ------------------------------------------------------------------------------
#
setup(
    name='radical.facts',
    version='0.2.0',
    description='Framework for Assessing Changes To Sea-level (FACTS)',
    url='https://github.com/radical-collaboration/facts.git',
    author='Gregory Garner',
    author_email='gregory.garner@rutgers.edu',
    license='unlicensed',
    keywords='radical facts',
    classifiers=[
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Utilities',
        'Topic :: System :: Distributed Computing',
        'Topic :: Scientific/Engineering',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix'
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    zip_safe=False,
    data_files=df,
)

# ------------------------------------------------------------------------------
#
os.system('rm -rf %s.egg-info' % 'radical.facts')


# ------------------------------------------------------------------------------

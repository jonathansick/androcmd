#!/usr/bin/env python

import os
import sys
from glob import glob

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://androcmd.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='androcmd',
    version='0.1.0',
    description='ANDROIDS star formation history inference from resolved stars',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Jonathan Sick',
    author_email='jonathansick@mac.com',
    url='https://github.com/jonathansick/androcmd',
    packages=[
        'androcmd',
    ],
    package_dir={'androcmd': 'androcmd'},
    include_package_data=True,
    package_data={'andrcomd': ['data/*.json']},
    install_requires=[
        'numpy',
        'astropy',
        'matplotlib',
        'astroml',
    ],
    license='MIT',
    zip_safe=False,
    keywords='androcmd',
    scripts=glob('scripts/*.py') + glob('scripts/starfish_magphys/*.py'),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)

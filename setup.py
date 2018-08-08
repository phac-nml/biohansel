# -*- coding: utf-8 -*-

import sys
from distutils.core import setup
from glob import glob

import io
from os.path import basename, splitext
from os.path import dirname
from os.path import join
from setuptools import find_packages
from setuptools.command.test import test as TestCommand

from biohansel import __version__, program_name, program_summary


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: Apache Software License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.6
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


# The following taken from https://docs.pytest.org/en/latest/goodpractices.html#manual-integration
class PyTest(TestCommand):
    """Required to be able to run tests with `python setup.py test`"""
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest

        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)


setup(
    name=program_name,
    version=__version__,
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/{}'.format(program_name),
    license='Apache v2.0',
    author='Peter Kruczkiewicz',
    author_email='peter.kruczkiewicz@gmail.com',
    description=program_summary,
    long_description=read('README.rst'),
    keywords=[
        'Salmonella',
        'enterica',
        'Heidelberg',
        'Enteritidis',
        'SNP',
        'kmer',
        'subtyping',
        'Aho-Corasick'
    ],
    classifiers=classifiers,
    package_dir={program_name: program_name},
    package_data={program_name: ['subtype/data/*/*.fasta', ]},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'numpy>=1.12.1',
        'pandas>=0.20.1',
        'pyahocorasick>=1.1.6',
        'attrs>=17.2.0',
        'click>=6.7',
        'scipy>=1.1.0',
        'biopython>=1.72',
        'ete3>=3.1.1'
    ],
    tests_require=['pytest', ],
    cmdclass={'test': PyTest},  # run tests with `python setup.py test`
    entry_points={
        'console_scripts': [
            f'{program_name}={program_name}.cli:cli',
        ],
    }
)

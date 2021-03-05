#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    'numpy>=1.12.1',
    'pandas>=0.20.1',
    'pyahocorasick>=1.1.6',
    'attrs',
    'rich'
]

setup(
    author='Peter Kruczkiewicz',
    author_email='peter.kruczkiewicz@gmail.com',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering'
    ],
    description='Subtype microbial whole-genome sequencing (WGS) '
                'data using SNV targeting k-mer subtyping schemes.',
    entry_points={'console_scripts': [
        'hansel=bio_hansel.main:main',
        'biohansel=bio_hansel.main:main',
    ]},
    install_requires=requirements,
    keywords='Salmonella enterica Heidelberg Enteritidis SNP kmer subtyping Aho-Corasick',
    license='Apache Software License 2.0',
    long_description=readme,
    name='bio_hansel',
    package_data={'bio_hansel': ['data/*/*.fasta', 'data/*/*.tsv',]},
    packages=find_packages(exclude=['test_*.py', 'tests']),
    test_suite='tests',
    url='https://github.com/phac-nml/biohansel',
    version='2.6.1',
    zip_safe=False,
)

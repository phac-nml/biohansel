|logo|

|conda| |nbsp| |pypi| |nbsp| |license| |nbsp| |nbsp| Master:|citest-master| |nbsp| Development:|citest-dev|



.. |logo| image:: logo.png
    :target: https://github.com/phac-nml/biohansel
.. |pypi| image:: https://badge.fury.io/py/bio-hansel.svg
    :target: https://pypi.python.org/pypi/bio_hansel/
.. |license| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
    :target: http://www.apache.org/licenses/LICENSE-2.0
.. |citest-dev|  image:: https://travis-ci.org/phac-nml/biohansel.svg?branch=development
    :target: https://travis-ci.org/phac-nml/biohansel
.. |citest-master| image:: https://travis-ci.org/phac-nml/biohansel.svg?branch=master
    :target: https://travis-ci.org/phac-nml/biohansel
.. |conda|   image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: https://bioconda.github.io/recipes/bio_hansel/README.html
.. |nbsp| unicode:: 0xA0
    :trim:

Subtype microbial whole-genome sequencing (WGS) data using SNV targeting k-mer subtyping schemes.

Includes 33 bp k-mer SNV subtyping schemes for *Salmonella enterica* subsp. enterica serovar Heidelberg and Enteritidis genomes developed by Genevieve Labbe et al.

Works on genome assemblies (FASTA files) or reads (FASTQ files)! Accepts Gzipped FASTA/FASTQ files as input!


Citation
========

If you find the ``biohansel`` tool useful, please cite as:

.. epigraph::

    A robust genotyping scheme for *Salmonella enterica* serovar Heidelberg clones circulating in North America.
    Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Marisa Rankin, Matthew Gopez, Chad R. Laing, Philip Mabon, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Joanne MacKinnon, Roger Johnson, John H.E. Nash.
    [Manuscript in preparation]


Requirements and Dependencies
=============================

Each new build of ``biohansel`` is automatically tested on Linux using `Continuous Integration <https://travis-ci.org/phac-nml/bio_hansel/branches>`_. ``biohansel`` has been confirmed to work on Mac OSX (versions 10.13.5 Beta and 10.12.6) when installed with Conda_.

These are the dependencies required for ``biohansel``:

- Python_ (>=v3.6)
    - numpy_ >=1.12.1
    - pandas_ >=0.20.1
    - pyahocorasick_ >=1.1.6
    - attrs_


Installation
============

With Conda_
-----------

Install ``biohansel`` from Bioconda_ with Conda_ (`Conda installation instructions <https://bioconda.github.io/#install-conda>`_):

.. code-block:: bash

    # setup Conda channels for Bioconda and Conda-Forge (https://bioconda.github.io/#set-up-channels)
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    # install biohansel
    conda install bio_hansel

With pip_ from PyPI_
---------------------

Install ``biohansel`` from PyPI_ with pip_:

.. code-block:: bash

    pip install bio_hansel

With pip_ from Github
---------------------

Or install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/biohansel.git@master

Install into Galaxy_ (version >= 17.01)
---------------------------------------

Install ``biohansel`` from the main Galaxy_ toolshed:

https://toolshed.g2.bx.psu.edu/repository?repository_id=59b90ef18cc5dbbc&changeset_revision=4654c51dae72


Usage
=====

If you run ``hansel -h``, you should see the following usage statement:

.. code-block::

    usage: hansel [-h] [-s SCHEME] [--scheme-name SCHEME_NAME]
                  [-p forward_reads reverse_reads] [-i fasta_path genome_name]
                  [-D INPUT_DIRECTORY] [-o OUTPUT_SUMMARY]
                  [-O OUTPUT_TILE_RESULTS] [-S OUTPUT_SIMPLE_SUMMARY] [--force]
                  [--json] [--min-kmer-freq MIN_KMER_FREQ]
                  [--max-kmer-freq MAX_KMER_FREQ]
                  [--low-cov-depth-freq LOW_COV_DEPTH_FREQ]
                  [--max-missing-tiles MAX_MISSING_TILES]
                  [--min-ambiguous-tiles MIN_AMBIGUOUS_TILES]
                  [--low-cov-warning LOW_COV_WARNING]
                  [--max-intermediate-tiles MAX_INTERMEDIATE_TILES] [-t THREADS]
                  [-v] [-V]
                  [F [F ...]]

    Subtype microbial genomes using SNV targeting k-mer subtyping schemes.
    Includes schemes for Salmonella enterica spp. enterica serovar Heidelberg and Enteritidis subtyping.
    Developed by Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Marisa Rankin, Matthew Gopez, Chad R. Laing, Philip Mabon, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Marisa Rankin, Joanne MacKinnon, Roger Johnson, John H.E. Nash.

    positional arguments:
      F                     Input genome FASTA/FASTQ files (can be Gzipped)

    optional arguments:
      -h, --help            show this help message and exit
      -s SCHEME, --scheme SCHEME
                            Scheme to use for subtyping (built-in: "heidelberg",
                            "enteritidis"; OR user-specified:
                            /path/to/user/scheme)
      --scheme-name SCHEME_NAME
                            Custom user-specified SNP substyping scheme name
      -p forward_reads reverse_reads, --paired-reads forward_reads reverse_reads
                            FASTQ paired-end reads
      -i fasta_path genome_name, --input-fasta-genome-name fasta_path genome_name
                            fasta file path to genome name pair
      -D INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                            directory of input fasta files (.fasta|.fa|.fna) or
                            FASTQ files (paired FASTQ should have same basename
                            with "_\d\.(fastq|fq)" postfix to be automatically
                            paired) (files can be Gzipped)
      -o OUTPUT_SUMMARY, --output-summary OUTPUT_SUMMARY
                            Subtyping summary output path (tab-delimited)
      -O OUTPUT_TILE_RESULTS, --output-tile-results OUTPUT_TILE_RESULTS
                            Subtyping tile matching output path (tab-delimited)
      -S OUTPUT_SIMPLE_SUMMARY, --output-simple-summary OUTPUT_SIMPLE_SUMMARY
                            Subtyping simple summary output path
      --force               Force existing output files to be overwritten
      --json                Output JSON representation of output files
      --min-kmer-freq MIN_KMER_FREQ
                            Min k-mer freq/coverage
      --max-kmer-freq MAX_KMER_FREQ
                            Max k-mer freq/coverage
      --low-cov-depth-freq LOW_COV_DEPTH_FREQ
                            Frequencies below this coverage are considered low
                            coverage
      --max-missing-tiles MAX_MISSING_TILES
                            Decimal proportion of maximum allowable missing tiles
                            before being considered an error. (0.0 - 1.0)
      --min-ambiguous-tiles MIN_AMBIGUOUS_TILES
                            Minimum number of missing tiles to be considered an
                            ambiguous result
      --low-cov-warning LOW_COV_WARNING
                            Overall tile coverage below this value will trigger a
                            low coverage warning
      --max-intermediate-tiles MAX_INTERMEDIATE_TILES
                            Decimal proportion of maximum allowable missing tiles
                            to be considered an intermediate subtype. (0.0 - 1.0)
      -t THREADS, --threads THREADS
                            Number of parallel threads to run analysis (default=1)
      -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                            show debug info)
      -V, --version         show program's version number and exit




Example Usage
=============

Analysis of a single FASTA file
-------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta


Contents of ``results.tab``:

.. code-block::

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR1002850  heidelberg  2.2.2.2.1.4 2; 2.2; 2.2.2; 2.2.2.2; 2.2.2.2.1; 2.2.2.2.1.4  1037658-2.2.2.2.1.4; 2154958-2.2.2.2.1.4; 3785187-2.2.2.2.1.4   True        202 202 17  17  3   3   SRR1002850.fasta


Contents of ``match_results.tab``:

.. code-block::

    tilename    stitle  pident  length  mismatch    gapopen qstart  qend    sstart  send    evalue  bitscore    qlen    slen    seq coverage    is_trunc    refposition subtype is_pos_tile sample  file_path   scheme
    775920-2.2.2.2  NODE_2_length_512016_cov_46.4737_ID_3   100.0   33  0   0   1   33  474875  474907  2.0000000000000002e-11  62.1    33  512016  GTTCAGGTGCTACCGAGGATCGTTTTTGGTGCG   1.0 False   775920  2.2.2.2 True    SRR1002850  SRR1002850.fasta   heidelberg
    negative3305400-2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  276235  276267  2.0000000000000002e-11  62.1    33  427905  CATCGTGAAGCAGAACAGACGCGCATTCTTGCT   1.0 False   negative3305400 2.1.1.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3200083-2.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  170918  170950  2.0000000000000002e-11  62.1    33  427905  ACCCGGTCTACCGCAAAATGGAAAGCGATATGC   1.0 False   negative3200083 2.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3204925-2.2.3.1.5   NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  175760  175792  2.0000000000000002e-11  62.1    33  427905  CTCGCTGGCAAGCAGTGCGGGTACTATCGGCGG   1.0 False   negative3204925 2.2.3.1.5   False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3230678-2.2.2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  201513  201545  2.0000000000000002e-11  62.1    33  427905  AGCGGTGCGCCAAACCACCCGGAATGATGAGTG   1.0 False   negative3230678 2.2.2.1.1.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3233869-2.1.1.1.1   NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  204704  204736  2.0000000000000002e-11  62.1    33  427905  CAGCGCTGGTATGTGGCTGCACCATCGTCATTA   1.0 False   
    [Next 196 lines omitted.]


Analysis of a single FASTQ readset
----------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_forward.fastqsanger SRR5646583_reverse.fastqsanger


Contents of ``results.tab``:

.. code-block::

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR5646583  heidelberg  2.2.1.1.1.1 2; 2.2; 2.2.1; 2.2.1.1; 2.2.1.1.1; 2.2.1.1.1.1  1983064-2.2.1.1.1.1; 4211912-2.2.1.1.1.1    True        202 202 20  20  2   2   SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger


Contents of ``match_results.tab``:

.. code-block::

    seq freq    sample  file_path   tilename    is_pos_tile subtype refposition is_kmer_freq_okay   scheme
    ACGGTAAAAGAGGACTTGACTGGCGCGATTTGC   68  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    21097-2.2.1.1.1 True    2.2.1.1.1   21097   True    heidelberg
    AACCGGCGGTATTGGCTGCGGTAAAAGTACCGT   77  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    157792-2.2.1.1.1    True    2.2.1.1.1   157792  True    heidelberg
    CCGCTGCTTTCTGAAATCGCGCGTCGTTTCAAC   67  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    293728-2.2.1.1  True    2.2.1.1 293728  True    heidelberg
    GAATAACAGCAAAGTGATCATGATGCCGCTGGA   91  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    607438-2.2.1    True    2.2.1   607438  True    heidelberg
    CAGTTTTACATCCTGCGAAATGCGCAGCGTCAA   87  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    691203-2.2.1.1  True    2.2.1.1 691203  True    heidelberg
    CAGGAGAAAGGATGCCAGGGTCAACACGTAAAC   33  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    944885-2.2.1.1.1    True    2.2.1.1.1   944885  True    heidelberg
    [Next 200 lines omitted.]

Analysis of all FASTA/FASTQ files in a directory
------------------------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/


``hansel`` will only attempt to analyze the FASTA/FASTQ files within the specified directory and will not descend into any subdirectories!


Development
===========


Get the latest development code using Git from GitHub:

.. code-block:: bash

    git clone https://github.com/phac-nml/biohansel.git
    cd biohansel/
    git checkout development
    # Create a virtual environment (virtualenv) for development
    virtualenv -p python3 .venv
    # Activate the newly created virtualenv
    source .venv/bin/activate
    # Install biohansel into the virtualenv in "editable" mode
    pip install -e .


Run tests with pytest_:

.. code-block:: bash

    # In the biohansel/ root directory, install pytest for running tests
    pip install pytest
    # Run all tests in tests/ directory
    pytest
    # Or run a specific test module
    pytest -s tests/test_qc.py



Legal
=====

Copyright Government of Canada 2017

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

Contact
=======

**Gary van Domselaar**: gary.vandomselaar@phac-aspc.gc.ca


.. _PyPI: https://pypi.org/project/bio-hansel/
.. _Conda: https://conda.io/docs/
.. _Bioconda: https://bioconda.github.io/
.. _pip: https://pip.pypa.io/en/stable/quickstart/
.. _numpy: http://www.numpy.org/
.. _pandas: http://pandas.pydata.org/
.. _pyahocorasick: http://pyahocorasick.readthedocs.io/en/latest/
.. _attrs: http://www.attrs.org/en/stable/
.. _Python: https://www.python.org/
.. _Galaxy: https://galaxyproject.org/
.. _pytest: https://docs.pytest.org/en/latest/

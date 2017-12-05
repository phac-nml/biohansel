|logo|
======

|pypi| |nbsp| |license| |citest| |conda| |nbsp|

.. |logo| image:: https://s2.postimg.org/5m5hiakax/Selection_044.png
    :target: https://github.com/phac-nml/bio_hansel
.. |pypi| image:: https://badge.fury.io/py/bio-hansel.svg
    :target: https://pypi.python.org/pypi/bio_hansel/
.. |license| image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0
.. |citest|  image:: https://travis-ci.org/phac-nml/bio_hansel.svg?branch=master
    :target: https://travis-ci.org/phac-nml/bio_hansel
.. |conda|   image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: https://bioconda.github.io/recipes/bio_hansel/README.html
.. |nbsp| unicode:: 0xA0
    :trim:

Subtype *Salmonella enterica* subsp. enterica serovar Heidelberg and Enteritidis genomes using *in-silico* 33 bp k-mer SNP subtyping schemes developed by Genevieve Labbe et al.

Subtype *Salmonella* genome assemblies (FASTA files) and/or whole-genome sequencing reads (FASTQ files)!

Citation
========

If you find this tool useful, please cite as:

.. epigraph::

    A robust genotyping scheme for *Salmonella enterica* serovar Heidelberg clones circulating in North America.
    Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Chad R. Laing, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Marisa Rankin, Joanne MacKinnon, Roger Johnson, John H.E. Nash.
    [Manuscript in preparation]


Requirements and Dependencies
=============================

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the external dependencies required for ``bio_hansel``:

- Python (>=v3.5)
- BLAST+ (>=v2.6)
    + for subtyping of FASTA input (assemblies)
- JELLYFISH (v2.0) (http://www.cbcb.umd.edu/software/jellyfish/)
    + for subtyping of FASTQ input (raw reads)


Installation
============

Ensure that BLAST+ and/or JELLYFISH are installed and accessible in your ``$PATH``.

Install ``bio_hansel`` from Conda:

.. code-block:: bash

    conda install bio_hansel

Install ``bio_hansel`` from PyPI:

.. code-block:: bash

    pip install bio_hansel

Or install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/bio_hansel.git@master


Usage
=====

If you run ``hansel -h``, you should see the following usage statement:

.. code-block:: none

          usage: hansel [-h] [-s SCHEME] [--scheme-name SCHEME_NAME]
                        [-p forward_reads reverse_reads] [-i fasta_path genome_name]
                        [-D INPUT_DIRECTORY] [-o OUTPUT_SUMMARY]
                        [-O OUTPUT_TILE_RESULTS] [-S OUTPUT_SIMPLE_SUMMARY] [--force]
                        [--min-kmer-freq MIN_KMER_FREQ] [--max-kmer-freq MAX_KMER_FREQ]
                        [--low-cov-depth-freq LOW_COV_DEPTH_FREQ]
                        [--max-missing-tiles MAX_MISSING_TILES]
                        [--min-ambiguous-tiles MIN_AMBIGUOUS_TILES]
                        [--max-intermediate-tiles MAX_INTERMEDIATE_TILES] [-t THREADS]
                        [-T TMP_DIR] [-K] [-v] [-V]
                        [F [F ...]]

          Subtype Salmonella enterica genomes using 33bp k-mer typing schemes.
          Includes schemes for Heidelberg and Enteritidis subtyping.
          Developed by Genevi�ve Labb�, James Robertson, Peter Kruczkiewicz, Chad R. Laing, Kim Ziebell, Marisa Rankin, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Eduardo N. Taboada, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Joanne MacKinnon, Benjamin M. Hetman, Roger Johnson, John H.E. Nash.

          positional arguments:
            F                     Input genome FASTA/FASTQ files

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
                                  paired)
            -o OUTPUT_SUMMARY, --output-summary OUTPUT_SUMMARY
                                  Subtyping summary output path (tab-delimited)
            -O OUTPUT_TILE_RESULTS, --output-tile-results OUTPUT_TILE_RESULTS
                                  Subtyping tile matching output path (tab-delimited)
            -S OUTPUT_SIMPLE_SUMMARY, --output-simple-summary OUTPUT_SIMPLE_SUMMARY
                                  Subtyping simple summary output path
            --force               Force existing output files to be overwritten
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
            --max-intermediate-tiles MAX_INTERMEDIATE_TILES
                                  Decimal proportion of maximum allowable missing tiles
                                  to be considered an intermediate subtype. (0.0 - 1.0)
            -t THREADS, --threads THREADS
                                  Number of parallel threads to run analysis (default=1)
            -T TMP_DIR, --tmp-dir TMP_DIR
                                  Base temporary working directory for intermediate
                                  analysis files
            -K, --keep-tmp        Keep temporary analysis files
            -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                                  show debug info)
            -V, --version         show program's version number and exit




Quality Checking
================

`bio_hansel` runs quality checking on files passed to it to provide feedback to the end user how the analysis went.
Each analysis will have either a `QC_STATUS`: PASS, WARNING, or FAIL followed by the error code if one occured.

Error Codes
-----------
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Type of Error         | Error Code | Description                                                                                     | How to Proceed                              |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Missing Tiles         | Error 1    | 5% of scheme tiles not found                                                                    | 1) Need more WGS data                       |
|                       |            |                                                                                                 | 2) Wrong serovar/species                    |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Mixed Subtype         | Error 2    | 1) 2 subtypes present in final subtype call                                                     | 1) Potential Contaminated Data              |
|                       |            | 2) There are conflicting tiles that match (+/-) for the same position                           | 2) Re-isolate Colony & Re-sequence          |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Ambiguous Results     | Error 3    | <5% scheme tiles are missing but there +/- tiles missing                                        | 1) Possible low geneome coverage            |
|                       |            | for target sites.                                                                               | 2) Possible recombination events            |
|                       |            | Ex. Final Subtype call = 2.1.2 but we're missing 2.1's subtyping tiles                          |                                             |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Non Confident Results | Error 4    | We have a final subtype call, but further downstream subtype's tiles are not present.           | 1) Need more WGS data                       |
|                       |            | Ex. Final subtype = 2.1.2 but we're missing 2.1.2.X                                             | 2) Re-sequence                              |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+
| Intermediate Subtypes | Warning    | Issue with subtyping scheme itself, where the analysis falls between two chairs of the scheme.  | 1) Requires further analysis of the scheme. |
|                       |            | Requires further analysis                                                                       |                                             |
+-----------------------+------------+-------------------------------------------------------------------------------------------------+---------------------------------------------+

Parameters
----------
The Quality Checking module within `bio_hansel` contains parameters which you can use to fine tune your quality checking results.

.. code-block:: none

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
            --max-intermediate-tiles MAX_INTERMEDIATE_TILES
                                  Decimal proportion of maximum allowable missing tiles
                                  to be considered an intermediate subtype. (0.0 - 1.0)


Example Usage
=============

Analysis of a single FASTA file
-------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta


Contents of ``results.tab``:

.. code-block:: none

    sample  scheme  scheme_version  subtype all_subtypes  tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes n_tiles_matching_all  n_tiles_matching_all_expected n_tiles_matching_positive n_tiles_matching_positive_expected  n_tiles_matching_subtype  n_tiles_matching_subtype_expected file_path qc_status qc_message  
    SRR1002850  heidelberg  0.5.0 2.2.2.2.1.4 2;  2.2;  2.2.2;  2.2.2.2;  2.2.2.2.1;  2.2.2.2.1.4 1037658-2.2.2.2.1.4;  2154958-2.2.2.2.1.4;  3785187-2.2.2.2.1.4 True  202 202 17  17  3 3 SRR1002850.fasta  PASS  


Contents of ``match_results.tab``:

.. code-block:: none

    tilename  stitle  refposition subtype is_pos_tile sample  file_path scheme  scheme_version  qc_status qc_message
    775920-2.2.2.2  NODE_2_length_512016_cov_46.4737_ID_3 775920  2.2.2.2 True  SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3305400-2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5 3305400 2.1.1.1 False SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3200083-2.1 NODE_3_length_427905_cov_48.1477_ID_5 3200083 2.1 False SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3204925-2.2.3.1.5 NODE_3_length_427905_cov_48.1477_ID_5 3204925 2.2.3.1.5 False SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3230678-2.2.2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5 3230678 2.2.2.1.1.1 False SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3233869-2.1.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5 3233869 2.1.1.1.1 False SRR1002850_smalltestdata  /home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    negative3254229-2.2.3.1.3 NODE_3_length_427905_cov_48.1477_ID_5 3254229 2.2.3.1.3 False SRR1002850_smalltestdata/home/mgopez/hansel_test_data/old_test_data/SRR1002850_smalltestdata.fasta  heidelberg  0.5.0 PASS  
    [Following lines ommitted.]


Analysis of a single FASTQ readset
----------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_forward.fastqsanger SRR5646583_reverse.fastqsanger


Contents of ``results.tab``:

.. code-block:: none

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR5646583  heidelberg  2.2.1.1.1.1 2; 2.2; 2.2.1; 2.2.1.1; 2.2.1.1.1; 2.2.1.1.1.1  1983064-2.2.1.1.1.1; 4211912-2.2.1.1.1.1    True        202 202 20  20  2   2   SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger


Contents of ``match_results.tab``:

.. code-block:: none

    seq freq    sample  file_path   tilename    is_pos_tile subtype refposition is_kmer_freq_okay   scheme
    ACGGTAAAAGAGGACTTGACTGGCGCGATTTGC   68  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    21097-2.2.1.1.1 True    2.2.1.1.1   21097   True    heidelberg
    AACCGGCGGTATTGGCTGCGGTAAAAGTACCGT   77  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    157792-2.2.1.1.1    True    2.2.1.1.1   157792  True    heidelberg
    CCGCTGCTTTCTGAAATCGCGCGTCGTTTCAAC   67  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    293728-2.2.1.1  True    2.2.1.1 293728  True    heidelberg
    GAATAACAGCAAAGTGATCATGATGCCGCTGGA   91  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    607438-2.2.1    True    2.2.1   607438  True    heidelberg
    CAGTTTTACATCCTGCGAAATGCGCAGCGTCAA   87  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    691203-2.2.1.1  True    2.2.1.1 691203  True    heidelberg
    CAGGAGAAAGGATGCCAGGGTCAACACGTAAAC   33  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    944885-2.2.1.1.1    True    2.2.1.1.1   944885  True    heidelberg
    [Next 200 lines omitted.]


Contents of ``tech_results.tab``:

..code-block:: none

    sample  subtype qc_status qc_message
    SRR1002850 2.2.2.2.1.4 PASS  


Analysis of all FASTA/FASTQ files in a directory
------------------------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/


``hansel`` will only attempt to analyze the FASTA/FASTQ files within the specified directory and will not descend into any subdirectories!


License
=======

Copyright 2017 Public Health Agency of Canada

Distributed under the GNU Public License version 3.0

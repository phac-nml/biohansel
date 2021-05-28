.. Bio_Hansel documentation master file, created by
   sphinx-quickstart on Wed May 23 15:05:31 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Biohansel's Read the Docs!
=====================================

**Currently debating on changing Subtype to Genotype in all outputs and as such, they are synonymous for biohansel in the current version (2.4.0) and in this manual**

Biohansel genotypes clonal microbial whole-genome sequencing (WGS) data using SNV targeting k-mer genotyping schemes.

This tool works on genome assemblies (FASTA files) or reads (FASTQ files)! Accepts Gzipped FASTA/FASTQ files as input!

Biohansel includes 33 base-pair k-mer SNV genotyping schemes focused on *Salmonella enterica* subsp. enterica serovars. K-mers can be any length as long as they are a odd number.

Currently, there are subtyping schemes for the following *Salmonella* serovars:

- Heidelberg

- Enteritidis

- Typhimurium

- Typhi

    - Typhi scheme adapted from Wong et al paper titled: `"An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid" <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5059462/>`_.

These genotyping schemes have been created, maintained, and/or adapted by Genevieve Labbe et al.

There is also an included *Mycobacterium tuberculosis* scheme that was modified from the Francesc Coll et al. paper titled:
`"A robust SNP barcode for typing Mycobacterium tuberculosis complex strains" <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4166679/>`_ 
for biohansel use by Daniel Kein.


Code is available on GitHub under https://github.com/phac-nml/biohansel. 

* :ref:`user-docs`
* :ref:`installation-docs`
* :ref:`legal`
* :ref:`contact`


.. _user-docs:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   user-docs/usage
   user-docs/Tutorial
   user-docs/input
   user-docs/genotyping_schemes
   user-docs/degenerate_base_expansion
   user-docs/output
   user-docs/parameters
   user-docs/command-line

.. _installation-docs:

.. toctree::
   :maxdepth: 2
   :caption: Installation

   installation-docs/home
   installation-docs/command-line
   installation-docs/galaxy
   installation-docs/versions

.. _legal:

.. toctree::
   :maxdepth: 2
   :caption: Legal

   legal

.. _contact:

.. toctree::
   :maxdepth: 2
   :caption: Contact

   contact

.. toctree::
   :maxdepth: 2
   :caption: poop

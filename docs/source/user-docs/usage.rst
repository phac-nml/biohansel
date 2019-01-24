Usage
=====

BioHansel subtypes microbial whole-genome sequencing (WGS) data using singe nucleotide variant (SNV) k-mer subtyping schemes.

SNV k-mer schemes can be created for use by the tool or the included 33 bp k-mer SNV subtyping schemes for Salmonella enterica subsp. enterica serovar Heidelberg and Enteritidis genomes developed by Genevieve Labbe et al. can be used for analysis of Salmonella WGS data. 

BioHansel can be installed with Conda, pip, or within an existing Galaxy infrastructure. View the `install guide <../installation-docs/home.html>`_ of your preference for additional details.

Requirements and Dependencies
-----------------------------

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the dependencies required for BioHansel

- Python_ (>=v3.5)
    - numpy_ >=1.12.1
    - pandas_ >=0.20.1
    - pyahocorasick_ >=1.1.6
    - attrs_

Quick Installation
------------------

With Conda_
-----------

Conda is the easiest way to install BioHansel and run it through the use of the command line.

First, install Conda_ (`Conda installation instructions <https://bioconda.github.io/#install-conda>`_).

Then, install ``bio_hansel``through Bioconda_ (64bit linux and MAC OSX) using the following commands:

.. code-block:: bash

    # OPTIONAL: To create a new Conda environment input this command on terminal:
    conda create -n "name of environment" python=3.6
    # Then to create/activate conda environment: (*note* name of environment is what user decides to name environment)
    source activate "name of environment"
    
    # You can then install bio_hansel in the new environment
    # To deactivate environment, input:
    source deactivate
    
    # Setup Conda channels for Bioconda and Conda-Forge (https://bioconda.github.io/#set-up-channels)
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    #Activate wanted Conda environment (base or user created)
    conda activate

    # Install bio_hansel
    conda install bio_hansel

    #Check installation with the following command; make sure to be in the correct environment
    hansel -h
    #This will display the usage statement

Remeber to activate the Conda environment that BioHansel is installed at to allow it to run.

With pip_ from PyPI_
---------------------

Install BioHansel from PyPI_ with pip_:

.. code-block:: bash

    pip install bio_hansel

This will install BioHansel along with the required dependencies.

Check that installation is correct with the command:

 .. code-block:: bash

    hansel -h
    #This will display the usage statement.

With pip_ from Github
---------------------

Install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/bio_hansel.git@master

Check that BioHansel is working with the command:

 .. code-block:: bash

    hansel -h
    #This will display the usage statement.

Install into Galaxy_ (version >= 17.01)
---------------------------------------

Galaxy admins install BioHansel from the main Galaxy toolshed (`tutorial <https://galaxyproject.org/admin/tools/add-tool-from-toolshed-tutorial/>`_):

https://toolshed.g2.bx.psu.edu/repository?repository_id=59b90ef18cc5dbbc&changeset_revision=4654c51dae72

Users can download and set up their own instance of Galaxy following the `get Galaxy tutorial <https://galaxyproject.org/admin/get-galaxy/>`_ and then install BioHansel from the toolshed as an admin using the admin instructions linked above.

Input Data
----------

BioHansel uses genome assemblies (FASTA files) or reads (FastQ files) from WGS data as an input. It also accepts these files as their Gzipped FASTA/FASTQ formats. Genomes can be fully assembled or a collection of contigs when analyzed without impacting the output.

SNV subtyping schemes have to be defined for BioHansel to run correctly. Two schemes are included and user created schemes can be used by creating SNV k-mers in the FASTA format. See `Creating schemes <subtyping_schemes.html>`_ for more details.

Subtype metadata schemes can be added to the analysis but are not required. Added metadata is joined with the subtype field of the final results. More detailed info on formatting of metadata schemes can be found in the `input section <input.html>`_. 

Output Results
--------------

Output of results will be found in three .tab files in the same folder as the analyzed data was stored or in the Galaxy histories window after analysis. The three files include:

- tech_results.tab
- results.tab
- match_results.tab

Detailed info about the results output can be found in the `output section <output.html>`_.

Parameters
----------

Parameters can be modified for users of both Galaxy and the command line. These can be changed based on the users need. Modifiable parameters include:

- K-mer Frequency Thresholds
    - Min k-mer frequency/coverage (default 8)
    - Max k-mer frequency/coverage (default 1000)

- Quality Checking Thresholds
    - QC: Frequency below this coverage are considered low coverage (default 20)
    - QC: Min number of tiles missing for Ambiguous Result (default 3)
    - QC: Decimal Proportion of max allowed missing tiles (default 0.05)
    - QC: Decimal Proportion of max allowed missing tiles for an intermediate subtype (default 0.05)
    - QC: Overall tile coverage below this value will trigger a low coverage warning (default 20)

Detailed info on BioHansels parameters and their functions can be found in the `parameter section <parameters.html>`_.


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


Usage
=====

Biohansel subtypes microbial whole-genome sequencing (WGS) data using single nucleotide variant (SNV) k-mer subtyping schemes.

SNV k-mer schemes can be created for use by the tool or the included 33 bp k-mer SNV subtyping schemes for Salmonella enterica subsp. enterica serovar Heidelberg and Enteritidis genomes developed by Genevieve Labbe et al. can be used for analysis of Salmonella WGS data. 

Biohansel can be installed with Conda, pip, or within an existing Galaxy infrastructure. View the `install guide <../installation-docs/home.html>`_ of your preference for additional details.

Requirements and Dependencies
-----------------------------

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the dependencies required for biohansel

- Python_ (>=v3.5)
    - numpy_ >=1.12.1
    - pandas_ >=0.20.1
    - pyahocorasick_ >=1.1.6
    - attrs_

Quick Installation
------------------

With Conda_
-----------

Conda is the easiest way to install and biohansel through the use of the command line.

First, install Conda_ (`Conda installation instructions <https://bioconda.github.io/#install-conda>`_).

Then, install ``biohansel`` through Bioconda_ (64bit linux and MAC OSX) using the following commands:

.. code-block:: bash

    # OPTIONAL: To create a new Conda environment input this command on terminal:
    conda create -n "name of environment" python=3.6
    # Then to create/activate conda environment: (*note* name of environment is what user decides to name environment)
    source activate "name of environment"
    
    # You can then install biohansel in the new environment
    # To deactivate environment, input:
    source deactivate
    
    # Setup Conda channels for Bioconda and Conda-Forge (https://bioconda.github.io/#set-up-channels)
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    #Activate wanted Conda environment (base or user created)
    conda activate

    # Install biohansel
    conda install bio_hansel

    #Check installation with the following command; make sure to be in the correct environment
    hansel -h
    #This will display the usage statement

Remember to activate the Conda environment that biohansel is installed into each time you want to run it after opening a new terminal window.

With pip_ from PyPI_
---------------------

Install biohansel from PyPI_ with pip_:

.. code-block:: bash

    pip install bio_hansel

This will install biohansel along with the required dependencies.

Check that installation is correct with the command:


 .. code-block:: bash

    hansel -h
    #This will display the usage statement.


With pip_ from Github
---------------------

Install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/biohansel.git@master

Check that biohansel is working with the command:

 .. code-block:: bash

    hansel -h
    #This will display the usage statement.

Install into Galaxy_ (version >= 17.01)
---------------------------------------

Galaxy admins install biohansel from the main Galaxy toolshed (`tutorial <https://galaxyproject.org/admin/tools/add-tool-from-toolshed-tutorial/>`_):

https://toolshed.g2.bx.psu.edu/view/nml/biohansel/ba6a0af656a6

Users can download and set up their own instance of Galaxy following the `get Galaxy tutorial <https://galaxyproject.org/admin/get-galaxy/>`_ and then install ``biohansel`` from the toolshed as an admin using the admin instructions linked above.

Input Data
----------

Biohansel uses genome assemblies (FASTA files) or reads (FastQ files) from WGS data as an input. It also accepts these files as their Gzipped FASTA/FASTQ formats. Genomes can be fully assembled or a collection of contigs when analyzed without impacting the output.

SNV subtyping schemes have to be defined for biohansel to run correctly. Two schemes are currently included in biohansel and user created schemes can be developed by creating SNV k-mer pairs in the FASTA format. See `Creating schemes <subtyping_schemes.html>`_ for more details.

Subtype metadata schemes can be added to the analysis but are not required. Added metadata is joined with the subtype field of the final results. More detailed info on formatting of metadata schemes can be found in the `Input section <input.html>`_ along with additional information on all of the other input files biohansel can use. 

Output Results
--------------

Output of the results generated through biohansel will be found in three .tab files in the directory that biohansel was run from or in the Galaxy histories window after analysis is complete. The three files include:

- tech_results.tab -> Most basic results file
- results.tab -> More advanced information on the results generated
- match_results.tab -> All k-mer information used to generate the subtype

Detailed info about the results output can be found in the `output section <output.html>`_.

Parameters
----------

Parameters can be modified for users of both Galaxy and the command line. These can be changed based on the users need. Modifiable parameters include:

- K-mer Frequency Thresholds - **only apply to raw reads/.fastq datasets**
    - Min k-mer frequency/coverage (default 8, cannot lower past 8 in current build)
    - Max k-mer frequency/coverage (default 1000)

- Quality Checking Thresholds - Important parameters for the final results of the QC columns
    - QC: Frequency below this coverage are considered low coverage (default 20)
    - QC: Min number of tiles missing for Ambiguous Result (default 3)
    - QC: Decimal Proportion of max allowed missing tiles (default 0.05)
    - QC: Decimal Proportion of max allowed missing tiles for an intermediate subtype (default 0.05)
    - QC: Overall tile coverage below this value will trigger a low coverage warning (default 20)

Detailed info on biohansels parameters and their functions can be found in the `parameter section <parameters.html>`_.

Running biohansel
-----------------

More detailed information is available under the `Tutorial section <tutorial.html>`_, the `input section <input.html>`_, or the `Command Line section <command-line.html>`_.

A basic command to run biohansel on an assembled Heidelberg fasta file  would be:

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab -S tech_results.tab </path/to/data_file>

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


Command Line  
============ 

There are three ways to install the latest version of BioHansel to run analyses through the use of the command line. These instructions have been tested for linux terminal and thus, no comments on the other operating system instructions can be found here at the moment. MAC should be extremely similar once you install Miniconda following the MAC installation steps. BioHansel commands have not been tested for the Windows OS terminal.

BioHansel on Miniconda Linux Installation Instructions
------------------------------------------------------

Miniconda is a mini version of `Anaconda <https://conda.io/projects/conda/en/latest/glossary.html#anaconda-glossary>`_ that includes only conda and its dependencies. If you wish to `install Anaconda <https://docs.continuum.io/anaconda/install.html>`_ then follow the steps found on the Anaconda installation instructions page and join back at step 5. 

The steps below will be the fully detailed for Linux operating systems. You can find installation guides for Miniconda for the other operating systems at: https://conda.io/projects/conda/en/latest/user-guide/install/index.html


1. Download Miniconda - Python 3.7 from `Conda's Download page <https://conda.io/en/latest/miniconda.html>`_

- Select the correct installer for your operating system (either 64-bit or 32-bit)



2. Path to where the miniconda package was downloaded and run the following command:

.. code-block:: bash

    bash Miniconda-latest-Linux-x86_64.sh

You can use the tab key after typing the first few letters of the file to finish the rest.


3. Follow the prompts on the installer screens


4. Close and re-open the terminal window to make the changes take effect

5. Add BioConda to MiniConda with the following commands in the following order:

.. code-block:: bash

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

6. Create a new environment to run BioHansel from. This environment is going to be called "bio_hansel" but can be whatever you choose. Command:

.. code-block:: bash

    conda create -n bio_hansel python=3.6

7. Activate the newly created conda environment:

.. code-block:: bash

    source activate bio_hansel

8. Check that the channels added earlier are properly included in the new environment by repeating the command in step 5. This can be skipped if step 5 went well:

.. code-block:: bash

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

The output will either be nothing if the channel is re-added, or it will say that the channel is already present. Either output is good.

9. Install bio_hansel and all of its dependencies into this environment:

.. code-block:: bash

    conda install bio_hansel

10. Check that the installation was successful. Running the following command to display the usage statement:

.. code-block:: bash

    hansel -h

11. When you open a new terminal window to run BioHansel, remember to activate the environment you set it to before running a job or it will not work:

.. code-block:: bash

    conda activate bio_hansel

    # Then run the analysis you want to do. Example:
    hansel -s heidelberg -o results.tab STR13341

If there are problems running/installing BioHansel, check to see if any of the following are are occuring:

1. Make sure that that all of the system requirements are met for Miniconda on this page: https://conda.io/projects/conda/en/latest/user-guide/install/index.html#system-requirements

2. Check that the right python version/installation is being used. It should be found under the /Miniconda3/bin/python3.7/ directory if installed with Miniconda:

.. code-block:: bash

    which python

- If the wrong directory or python version is being run by the terminal, then try the following:

		- Restart the terminal window and check again
		
		- Use the following commands

			.. code-block:: bash

			    alias python=python3
			    # This will set python 3 as the working python
			    # Check that this worked with the command:
   			    python --version
			    # Should print out the version as 3.x.x depending on which version is installed.

3. If using the Fish shell, make sure that you add the following line to your ``fish.config`` file if there are problems occuring:

.. code-block:: bash

    source (conda info --root)/etc/fish/conf.d/conda.fish

|
BioHansel installation with pip from PyPI
-----------------------------------------

If you have pip and python3 installed already onto your machine, then the following steps can be used to install BioHansel. If not, follow along and install them as prompted:

1. Make sure that python3 is the active python version and that it is installed onto the machine:

.. code-block:: bash

    python -V
    # This will print the python version used

    # If this doesn't output 3.X and instead outputs 2.X, then type:
    alias python=python3
    
    python -V
    # Now it should output python version as 3.X. 

    # If not, then python3 may need to be installed with the following:
    apt-install python3.7-minimal

BioHansel needs python3 to work correctly. If installed this way, you may need to use the alias command to get the correct version of python active before each run of BioHansel

2. Install BioHansel with pip. If pip is not installed on your current machine, then follow the `installing pip tutorial <https://pip.pypa.io/en/stable/installing/>`_:

.. code-block:: bash
    
    # You can check that pip is installed with the input:
    pip

    # If pip is installed, then install bio_hansel with it
    pip install bio_hansel

    # This will install bio_hansel along with all of its needed dependencies. 

3. Check that BioHansel has been correctly installed with:

.. code-block:: bash

    hansel -h

Common problems encountered:

1. pip installing BioHansel to the wrong python environment. Instead of installing to python 3, it installs to python 2.
	
- Set the correct path for pip/python to install files

2. Make sure the correct version of python is being installed to (v3.x)


BioHansel installation with pip from Github
-------------------------------------------

Use the following command:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/bio_hansel.git@master

If that doesn't work, look at the common problems encountered with pip from PyPI or try the PyPI installation instructions. Both installation methods are extremely similar.


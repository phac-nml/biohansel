Usage
=====

Requirements and Dependencies
-----------------------------

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the dependencies required for ``bio_hansel``:

- Python_ (>=v3.5)
    - numpy_ >=1.12.1
    - pandas_ >=0.20.1
    - pyahocorasick_ >=1.1.6
    - attrs_

Installation
------------

With Conda_
-----------

Install ``bio_hansel`` from Bioconda_ with Conda_ (`Conda installation instructions <https://bioconda.github.io/#install-conda>`_):

.. code-block:: bash

    # setup Conda channels for Bioconda and Conda-Forge (https://bioconda.github.io/#set-up-channels)
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    # install bio_hansel
    conda install bio_hansel

With pip_ from PyPI_
---------------------

Install ``bio_hansel`` from PyPI_ with pip_:

.. code-block:: bash

    pip install bio_hansel

With pip_ from Github
---------------------

Or install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/phac-nml/bio_hansel.git@master

Install into Galaxy_ (version >= 17.01)
---------------------------------------

Install ``bio_hansel`` from the main Galaxy_ toolshed:

https://toolshed.g2.bx.psu.edu/repository?repository_id=59b90ef18cc5dbbc&changeset_revision=4654c51dae72

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
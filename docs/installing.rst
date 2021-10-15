.. _installation:
############################
Requirements & Installation
############################


Installation for Users
#########################

Getting ``aplc-optimization`` up and running on your own computer requires four steps, detailed below:

 1. Cloning the GitHub repository.
 2. Installing the ``conda`` environment.
 3. Installing and licensing Gurobi.
 4. Installing the python package.

.. _installing-with-conda:

Prerequisites
=============
For ease of installation, it is highly recommended that users have a working
installation of either Anaconda or Miniconda for Python 3.7.
Download and installation instructions can be found here:

 - `Miniconda <https://conda.io/miniconda.html>`_
 - `Anaconda <https://www.continuum.io/downloads>`_

Clone the ``aplc_optimization`` repository
==========================================
First, you will need to clone the current version of ``aplc_optimization``. The simplest way to do this is to go to the
directory you want your copy of the repository to be located and clone the repository there.
Once you are in the directory you can do the following::

    $ git clone https://github.com/spacetelescope/aplc_optimization.git
    $ cd aplc_optimization

Environment Installation
========================
The next step is to install the ``aplc_optimization`` conda environment via
the environment yaml file, which contains all of the dependencies for the project::

    $ conda env create --file environment.yml

and then activate the environment::

    $ conda activate aplc_optimization

Package Installation
====================
Next, you need to install the ``aplc_optimization`` package. This can be accomplished using ``pip``::

    $ pip install .

or by running the ``setup.py`` script::

    $ python setup.py install

You can check if ``aplc_optimization`` is installed correctly by importing it in Python::

    >>> import aplc_optimization

--------------------------------------------------

Installation for Contributors
#############################
To instead install ``aplc_optimization`` in development mode (i.e. so that edits made to the source files in the
``aplc_optimization/`` directory will change the library), use ``pip install`` with the ``-e`` option when
installing the package::

    pip install -e .

This allows for easy updates by simply pulling the git repository::

    git pull
    git setup.py egg_info



.. _installing-gurobi:

----------------------------------------

Installing Gurobi
#################

The ``aplc_optimization`` toolkit relies on the `Gurobi solver <https://www.gurobi.com/>`_, which it calls directly from
Python using the ``gurobipy`` package.


Register for an Academic account
================================

In order to download the Gurobi Optimizer you will need to first register for an account.
If you already have an account, `Log In <https://www.gurobi.com/login>`_; otherwise,
`register <https://pages.gurobi.com/registration>`_ for a free *Academic* account.

.. _download-gurobi:

Download the Gurobi optimizer
=============================

After registering and logging in to your Gurobi account, go to the
`Gurobi software download page <https://www.gurobi.com/downloads/gurobi-software/>`_. Find your platform
(we'll assume Mac OS X in this document) and choose the corresponding file to download. Once downloaded, double-click
on the appropriate Gurobi installer (e.g., ``gurobi9.1.2_mac64.pkg`` for Gurobi 9.1.2) and follow the prompts.
By default, the installer will place the Gurobi files in ``/Library/gurobi911/mac64`` (note that this is the system
``/Library`` directory, not your personal ``~/Library`` directory).

.. _get-gurobi-license:

Obtain a Gurobi license
=======================
In order to use the Gurobi Optimizer you will require a Gurobi license. Once you have downloaded the Gurobi Optimizer,
as above, visit the the `Academic License page <https://www.gurobi.com/downloads/end-user-license-agreement-academic/>`_ to
request a free license (note: you will first need to read and agree to the End User License Agreement and the
Conditions for academic use; once you have done so, click on "Request License").

Install the Gurobi license
--------------------------
Your next step is to install this Gurobi license on your machine. Once your license is visible on the
`Current Gurobi Licenses page <https://www.gurobi.com/downloads/licenses/>`_, click on the *License ID*
to view the *License Detail* page.

To obtain a Gurobi license key you'll need to run the ``grbgetkey`` command on your machine. The exact ``grbgetkey`` command
to run for a specific license is indicated at the bottom of the License Detail page (e.g., ``grbgetkey 253e22f3-...``).
Copy the entire ``grbgetkey`` command and paste it into a Terminal window. Once run, the ``grpgetkey`` program will prompt you to store
the license key on your machine. You can store the license key file anywhere, but we strongly recommend that you accept
the default location by hitting Enter. Setting up a non-default location is error-prone and a frequent source of trouble.

.. note::

    If you would like to store your ``gurobi.key`` license file in a non-default location, you can do so by setting the **GRB_LICENSE_FILE** environment variable to point to the license key file location.

Test the Gurobi license
-----------------------
Once you have obtained a license key for your machine, you are ready to test your license using the Gurobi Interactive Shell.
To do this, type ``gurobi.sh`` in a Terminal window. The shell should produce the following output::

    Using license file /Library/gurobi/gurobi.lic
    Set parameter LogFile to value gurobi.log

    Gurobi Interactive Shell, Version 9.1.1
    Copyright (c) 2020, Gurobi Optimization, LLC
    Type "help()" for help

    gurobi>

If the Gurobi shell didn't produce the desired output, there's a problem with your license (see the Gurobi
`documentation <https://www.gurobi.com/documentation/9.1/quickstart_mac/testing_your_license.html#subsection:testlicense>`_ for more information).

------------------------------------------------------------

Software Requirements
#####################

See `the environment.yml specification file <https://github.com/spacetelescope/aplc_optimization/blob/scda_21/environment.yml>`_ for the required package dependencies.

**Required Python version**: ``aplc_optimization`` requires Python 3.7 or higher.

**Conda channels:**

 - AstroConda (http://ssb.stsci.edu/astroconda)
 - Gurobi (http://conda.anaconda.org/gurobi)
 - Conda-forge (https://anaconda.org/conda-forge)

**Major Python dependencies**

 - `hcipy <https://docs.hcipy.org/0.3.1/>`_ (for coronagraphic simulations)
 - `gurobi <https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html#section:Python>`_ (for building and running optimization models)
 - `numpy <>`_ (for all numerical calculations)
 - `matplotlib <http://matplotlib.org>`_ (for visualizations)
 - `Astropy <http://astropy.org>`_ (for fits file reading and writing)
 - `asdf <https://pypi.org/project/asdf/>`_ (for reading and writing of HCIpy objects)
 - `imageio <https://pypi.org/project/imageio/>`_ (for writing image data)
 - `SciPy <http://www.scipy.org/scipylib/download.html>`_ (for advanced linear algebra)



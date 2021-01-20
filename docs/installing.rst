.. _installation:

Requirements & Installation
===========================

.. warning::

   Pardon our dust, this space is under construction.

.. _installing-with-conda:

Installing with Conda
----------------------

For ease of installation, we recommend using Conda to install the ``aplc-optimization`` package as follows:

1. Clone the  ``aplc-optimization`` repository::

    $ git clone https://github.com/spacetelescope/aplc-optimization.git

2. Create a new ``conda`` environment using the ``environment.yml`` file::

    $ cd aplc_optimization
    $ conda env create --file environment.yml
    $ source activate aplc-optimization


.. _installing-gurobi:

Installing Gurobi
-------------------

The ``aplc-optimization`` toolkit relies on the `Gurobi solver <https://www.gurobi.com/>`, which it calls directly from Python using the ``gurobipy`` package.

Register for an academic account
'''''''''''''''''''''''''''''''''
In order to download the Gurobi Optimizer you will need to first register for an account.
If you already have an account, `Log In <https://www.gurobi.com/login>`_. Otherwise,
`Register <https://pages.gurobi.com/registration>`_ for a free *Academic* account.

.. _download-gurobi:

Download the Gurobi optimizer
'''''''''''''''''''''''''''''''
After registering and logging in, go to the `Gurobi software download page <https://www.gurobi.com/downloads/gurobi-software/>`_. Find your platform
(we'll assume Mac OS X in this document) and choose the corresponding file to download. Once downloaded, double-click on the appropriate Gurobi installer
(e.g., ``gurobi9.1.1_mac64.pkg`` for Gurobi 9.1.1) and follow the prompts. By default, the installer will place the Gurobi files
in /Library/gurobi911/mac64 (note that this is the system /Library directory, not your personal ~/Library directory).

.. _get-gurobi-license:

Obtain a Gurobi license
''''''''''''''''''''''''''''
In order to use the Gurobi Optimizer, a Gurobi license is required. Once you have downloaded the Gurobi Optimizer, as above,
visit the the `Academic License page <https://www.gurobi.com/downloads/end-user-license-agreement-academic/>`_ to
request a free license (note: you will first need to read and agree to the End User License Agreement and the Conditions for academic use;
once you have done so, click on "Request License"). Your new license will be visible immediately on the
`Current Gurobi Licenses page <https://www.gurobi.com/downloads/licenses/>`_ and you can create as many academic licenses as you like.

Install the Gurobi license
```````````````````````````
Your next step is to install this Gurobi license on your machine. Once your license is visible on the
`Current page <https://www.gurobi.com/downloads/licenses/>`_ of the Gurobi website, click on the *License ID*
to view the License Detail page.

To obtain a Gurobi license key you'll need to run the ``grbgetkey`` command on your machine. The exact ``grbgetkey`` command
to run for a specific license is indicated at the bottom of the License Detail page (e.g., ``grbgetkey 253e22f3-...``).
Copy the entire ``grbgetkey`` command and paste it into a Terminal window. Once run, the ``grpgetkey`` program will prompt you to store
the license key on your machine. You can store the license key file anywhere, but we strongly recommend that you accept
the default location by hitting Enter. Setting up a non-default location is error-prone and a frequent source of trouble.

.. note::

    If you would like to store your ``gurobi.key`` license file in a non-default location, you can do so by setting the **GRB_LICENSE_FILE** environment variable to point to the license key file location.

Test the Gurobi license
''''''''''''''''''''''''
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
----------------------

See `the environment.yml specification file <https://github.com/spacetelescope/aplc-optimization/blob/scda_21/environment.yml>`_ for the required package dependencies.

**Required Python version**: ``aplc-optimization`` requires Python 3.7 or higher.

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


Using the STScI Linux Servers
------------------------------

The Linux servers at STScI are available to optimize select high-resolution design cases, as needed. For instructions on how to
install and use ``aplc-optimization`` on one of these servers, see :ref:`servers`.
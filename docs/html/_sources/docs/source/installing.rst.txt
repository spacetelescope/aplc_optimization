.. _installation:

Requirements & Installation
===========================

.. warning::

   Pardon our dust, this space is under construction.

.. _installing-with-conda:

Installing with Conda
----------------------

For ease of installation, we recommend using Conda to install the `aplc-optimization`` package as follows:

1. Clone the  ``aplc-optimization`` repository::

    $ git clone -b develop https://github.com/spacetelescope/aplc-optimization.git

2. Create a new ``conda`` environment using the ``environment.yml`` file::

    $ cd aplc_optimization
    $ conda env create --file environment.yml
    $ source activate aplc-optimization


.. _installing-gurobi:

Installing the Optimization Solver (Gurobi)
--------------------------------------------

Register for an academic account
'''''''''''''''''''''''''''''''''
In order to download the Gurobi Optimizer you will need to register for an account.
If you already have an account, `log in <https://www.gurobi.com/login>`_. Otherwise,
`register <https://pages.gurobi.com/registration>`_ for a free *Academic* account.

.. _download-gurobi:

Download Gurobi Optimizer
'''''''''''''''''''''''''''''''
Before using the Gurobi Optimizer, you'll need to install the software to your computer.
Go to the `Gurobi software download page <https://www.gurobi.com/downloads/gurobi-software/>`_, find your platform (we'll assume Mac OS X in this document) and
choose the corresponding file to download. In addition to the software, a `README <https://packages.gurobi.com/9.0/README.txt>`_ file containing
installation instructions is also available.

The next step is to double-click on the appropriate Gurobi installer (e.g., `gurobi7.5.2_mac64.pkg` for Gurobi 7.5.2)
and follow the prompts.

.. note::

    If you would like an overview of the files included in the Gurobi distribution,
    see the `File Overview  <https://www.gurobi.com/documentation/7.5/quickstart_mac/file_overview.html#section:Overview>`_ page on the Gurobi website.

.. _get-gurobi-license:

Obtaining a Gurobi License
''''''''''''''''''''''''''''
In order to use the Gurobi Optimizer, a Gurobi license is required. Once you have downloaded the Gurobi Optimizer, as above,
visit the the `Academic License page <https://www.gurobi.com/downloads/end-user-license-agreement-academic/>`_ to
request a free license (note: you will first need to read and agree to the End User License Agreement and the Conditions for academic use;
once you have done so, click on "Request License"). Your new license will be visible immediately on the
`Current Gurobi Licenses page <https://www.gurobi.com/downloads/licenses/>`_ and you can create as many academic licenses as you like.

Install the Gurobi License
```````````````````````````
The next step is to install this Gurobi license on your machine. Once your license is visible on the
`Current page <https://www.gurobi.com/downloads/licenses/>`_ of the Gurobi website, click on the *License ID*
to view the License Detail page and run the ``grpgetkey`` argument provided therein
(e.g. ``grbgetkey ae36ac20-16e6-acd2-f242-4da6e765fa0a``) to install the license where
the Gurobi Optimizer is installed.

.. note::

    The ``grbgetkey`` program will prompt you to store the license key on your machine, as well as validate your
    eligibility by confirming your academic domain (e.g., any ‘.edu’ address). It is strongly recommend that
    you place your client ``gurobi.lic`` file in a default location for your platform (either your home directory or ``/Library/gurobi``).
    Setting up a non-default location is error-prone and a frequent source of trouble.




Software Requirements
----------------------

See `the environment.yml specification file <https://github.com/spacetelescope/aplc-optimization/blob/scda_21/environment.yml>`_ for the required package dependencies.

**Required Python version**: aplc-optimization requires Python 3.7 or higher.

**Conda channels:**

 - AstroConda (http://ssb.stsci.edu/astroconda)
 - Gurobi (http://conda.anaconda.org/gurobi)
 - Conda-forge (https://anaconda.org/conda-forge)

**Python dependencies**

 - `hcipy <https://docs.hcipy.org/0.3.1/>`_ (for coronagraphic simulations)
 - `gurobi <https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html#section:Python>`_ (for building and running optimization models)
 - numpy (for all numerical calculations)
 - `matplotlib <http://matplotlib.org>`_ (for visualizations)
 - `Astropy <http://astropy.org>`_ (for fits file reading and writing)
 - asdf (for reading and writing of HCIpy objects)
 - imageio (for image writing)
 - `SciPy <http://www.scipy.org/scipylib/download.html>`_ (for advanced linear algebra)

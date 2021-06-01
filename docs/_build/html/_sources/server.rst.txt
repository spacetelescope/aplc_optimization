.. _servers:

###########################
Using STScI Linux Servers
###########################

Due to the wide parameter space open to the SCDA study, apodizer optimization jobs are often computationally expensive, with some high-resolution design cases requiring
several days to complete. In order to scale up the feasible volume of optimization trails, we utilize back room Linux compute servers at STScI.
There are three shared Telescope servers that can be used for the purpose of coronagraph design work: telserv1, telserv2 and telserv3.
They are physical machines, not virtual machines, and run the Linux operating system.

The Linux virtual workstations can be accessed from any STScI Macintosh or Windows system, on-site or through VPN and
are available to anyone with an Active Directory account. Information on the Telescopes shared computers can be found in
the following STScI `Innerspace page <https://innerspace.stsci.edu/display/INSTEL/Telescopes+Shared+Computers>`_.


Getting Started
====================
Verifying network mounts
--------------------------

In order to access both Central Storage, and the Linux servers, you must be on the STScI internal network
(either physically or wirelessly connected, or via VPN). Before proceeding, it may be useful to verify your
network mounts.

Open a terminal terminal and check that you can list the contents of the ``/user`` directory as shown::

        $ ls /user

If instead, you see something like this...::

    $ ls /user
    ls: /user: No such file or directory

... then you will need to talk to ITSD.

.. _Login:

Logging into a server
----------------------

The command to remotely connect to a given server is ``ssh <SERVERNAME>`` (where ``<SERVERNAME>`` is the server hostname,
e.g. ``ssh telserv3``). To enable X forwarding, you can use the ``-XY`` option (e.g. ``ssh -XY science1``). This allows
for displays that use X (Linux's default backend) to be forwarded to your monitor.

Go ahead and verify that you are able to log in to the TELSERV3 server, as follows::

    $ ssh telserv3

By default, the STScI servers run on ``tsch`` (not ``bash``) when you log in. To switch to a ``bash`` shell use::

    $ bash


Software Installation
======================

.. _Storage:

An aside on storage
'''''''''''''''''''

When logged into a server (see :ref:`below <Login>`), any commands run in that SSH session will be executed on the server, not on your personal machine.
Thus, files and directories on your local machine will not be accessible when running the servers. Instead, all of the
servers have access to network-mounted home directories (i.e. ``/home/<username>``) and access to the ``/user`` space
on Central Store.

Home directories
``````````````````

Each user has a home directory (``/home/<username>``) that is mounted on each of the Linux machines. This home directory
is mounted and shared among all the machines, so that any file in that directory will be available on each machine.
However, space in this directory is limited (each user has a 10 GB quota) and installing software in the home directory
is usually discouraged.

``/user`` space
`````````````````
A user's Central Storage space is a better choice than the home directory for user-installed software,
as it has considerably more space (at least 850 GB). This directory is shared and mounted on each of the linux machines,
and can also be mounted on the user's individual computer (when connected to the STScI internal network). This allows for
a convenient way to share and move files among/between the Linux machines and the user's individual machine.

On the servers, this space is available under the path ``/user/<username>``. For Mac users, this directory is
auto-mounted at the same location (i.e. ``/user/<username>``). On windows machines, the directory is not auto-mounted,
but is still available at the path ``\\cs10d\user``. When installing in this location, it is important to remain aware
of the differences in the operating system among the machines.

Within a SSH session (see :ref:`Logging into a server <Login>`), verify that you can access your Central Storage directory::

    $ ls /user/myname


Installing Conda
-----------------

It is recommended that the ``aplc_optimization`` package be installed using Conda (see :ref:`installing`). While you may
have a version of Conda installed to your local machine, this version will not be accessible to programs running on
the servers at STScI (see :ref:`above <Storage>`). Therefore, in order to install ``aplc_optimization`` on the servers,
we will first need set up Conda on one of the servers. (Fortunately, you only have to go through the process on *one*
Linux server to make conda available on *all* of the Linux servers.)

Because of the storage limitations associated with the Linux home diretory (see :ref:`above <Storage>`),
we recommend installing Conda your Central Store directory instead.

First ssh into one of the servers::

    $ ssh telserv3
    telserv3>

Start a ``bash`` shell:

    telserv3> bash
    bash$

Change directories to your Central Store directory::

    bash$ cd "/user/$(logname)"

Because the servers at STScI run Linux, you will need to install a *Linux* version of Conda. Download and install
the latest Miniconda for linux using curl::

    bash$ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash$ bash Miniconda3-latest-Linux-x86_64.sh -b -p "/user/$(logname)/miniconda3"

Now add the miniconda folder to the $PATH on the server::

    bash$ echo "export PATH=\"/user/$(logname)/miniconda3/bin:\$PATH\"" >> ~/.bashrc

Log out and ssh back into the server again::

    bash$ exit
    telserv3> exit
    logout
    Connection to telserv3.stsci.edu closed.
    $ ssh telserv3
    telserv3>

The ``conda`` command works with ``bash``, not ``tcsh``; so before you continue, start a bash shell on the server::

    telserv3> bash
    bash$

Now check that the ``conda`` command is available::

    bash$ which conda
    /user/myname/miniconda3/bin/conda

Now that you have the conda command available, we can install ``aplc_optimization`` on the server.

Installing ``aplc_optimization``
---------------------------------

Because of the storage limitations associated with the Linux home diretory (see :ref:`above <Storage>`), we recommend
intalling ``aplc_optimization`` in your Central Store directory. If you haven't already, ``ssh`` into the server and open a
``bash`` shell::

    ssh telserv3
    telserv3> bash
    bash$

Go to your user directory on Central Store and clone the ``_optimization`` repository::

    bash$ cd /user/$(logname)/
    bash$ git clone https://github.com/spacetelescope/aplc_optimization.git

Create a new Conda environment on the Server using the ``environment.yml`` file::

    bash$ cd aplc_optimization
    bash$ conda env create --file environment.yml
    bash$ source activate aplc_optimization
    (aplc_optimization) bash$


Installing Gurobi
----------------------

Go back into your Central Store directory and download the latest 64-bit Linux version of the Gurobi optimizer, using ``curl``::

    (aplc_optimization) bash$ cd /user/bnickson
    (aplc_optimization) bash$ curl -OL https://packages.gurobi.com/9.0/gurobi9.0.2_linux64.tar.gz

Open and extract the ``.tar.gz`` file::

    (aplc_optimization) bash$ tar xvzf gurobi9.0.2_linux64.tar.gz


Obtain a new Gurobi license
'''''''''''''''''''''''''''''
Academic users are able to install and license Gurobi for their own use on more than one machine, however each
individual academic license can only be installed on a single physical machine. Consequently, in order to use the optimizer
on any of the Linux machines, a new individual academic license must be obtained (aside from the license previously
obtained for your personal machine).

Log into the Gurobi website and visit the `Academic License page <https://www.gurobi.com/downloads/end-user-license-agreement-academic/>`_
and request a new free license. Once obtained, this new license will be visible in your
`Current Gurobi Licenses <https://www.gurobi.com/downloads/licenses/>`_ library.

Go to your `Current Gurobi Licenses <https://www.gurobi.com/downloads/licenses/>`_ library and open
the newly created license by clicking on the **Licence ID** number.

Retrieve and set up the Gurobi license on the server
''''''''''''''''''''''''''''''''''''''''''''''''''''''
The next step is to install the new Gurobi license on the Linux machine. First, log in to the server
and go to the ``/user`` directory::

    ssh telserv3
    telserv3> bash
    bash$ cd "/user/$(logname)"

After navigating to the **License Detail** page on the Gurobi website, copy and paste the ``grbgetkey`` command provided to
install the ``gurobi.lic`` file on the Linux machine.

.. warning::

    While both the Linux machines and your personal machines have access to the ``/user`` space on Central Store,
    it is important that the Gurobi license only be installed whilst connected to the server. When the
    ``grbgetkey`` command is run, it passes identifying information about the machine it is run on to the Gurobi key server, and the server responds with a
    license key registered to that particular machine. As such, if the ``grbgetkey`` command is run whilst not connected to the server,
    the corresponding license will be registered to your personal machine and will be invalid for use while connected to the server.

Test the license on the server
''''''''''''''''''''''''''''''''

Once you have obtained a license key for the Linux machines, you can test the license using the Gurobi interactive shell.
Log out and back in to the server and run the ``gurobi.sh`` command::

    telserv3> quit
    Logging out...
    $ ssh telserv3
    telserv3> gurobi.sh
    Using license file /user/($logname)/gurobi/gurobi.lic
    Academic license - for non-commercial use only

    Gurobi Interactive Shell (mac64), Version 9.0.1
    Copyright (c) 2020, Gurobi Optimization, LLC
    Type "help()" for help

    gurobi>

If the above output is produced by the Gurobi shell, the license is functioning correctly.

.. warning::

    If the Gurobi shell doesn't produce the desired output, there's a problem with your license. Visit the
    `Gurobi website <https://www.gurobi.com/documentation/9.0/quickstart_mac/testing_your_license.html#subsection:testlicense>`_
    for troubleshooting tips.




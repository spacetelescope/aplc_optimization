.. _servers:

###########################
Using STScI Linux Servers
###########################

Due to the wide parameter space open to the SCDA study, apodizer optimization jobs are often computationally expensive,
with some high-resolution design cases requiring several days to complete. In order to scale up the feasible volume of
optimization trails, we utilize back-room servers at STScI.
There are three shared Telescope servers that can be used for the purpose of coronagraph design work: telserv1,
telserv2 and telserv3. They are physical machines, not virtual machines, and run the Linux operating system.

The Linux virtual workstations can be accessed from any STScI Macintosh or Windows system, on-site or through VPN, and
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

In order to make files from your personal machine available on the back-room servers, and visa versa, you will need to
have your own Central Storage directory. To verify that your Central Storage directory exists, you can use the
following command in a terminal::

    ls /user/username

replacing ``username`` with your username. If you get a ``not found`` or ``permission denied`` error, you will have to
get in touch with ITSD.

.. _Login:

Logging into a server
----------------------

The command to remotely connect to a given server is ``ssh <SERVERNAME>`` (e.g. ``ssh telserv3``).  It will attempt to
log you in with your current username and prompt you for your password. To enable X Windows forwarding, you can use
the ``-XY`` options (e.g. ``ssh -XY telserv3``). This allows for displays that use X11 to be forwarded to your monitor.
The -X and -Y options do the same thing, with some slight differences; some remote applications may work better with
one than the other, so using both options can prevent unforseen headaches.

Verify that you are able to log in to the TELSERV3 server, as follows::

    $ ssh telserv3

By default, the STScI servers run on ``tsch`` (not ``bash``) when you log in. To switch to a ``bash`` shell use::

    $ bash

Note that any commands run in an ssh session will execute on the server, not your personal machine.

Software Installation
======================
.. _Storage:

Software installed on your local machine will not be accessible to programs running on the servers at STScI. One
important consequence of this is that you will be responsible for installing the same things into both environments.

Installing Conda
-----------------
It is recommended that the ``aplc_optimization`` package be installed using Conda (see :ref:`installing`). While you may
have a version of Conda installed to your local machine, as :ref:`previously <Storage>` mentioned, this version will
not be accessible to programs running on the servers at STScI (see :ref:`above <Storage>`).

All the servers have network-mounted home directories (``/home/<username>``), access to the Central Store
(``/user/<username>``), and locally-mounted disks (e.g. ``/internal/data1``). Conda wants to install in the home
directory by default, however due to disk quotas (each user has a quota of 50GB) it is generally discouraged to install
software to your home directory. Your Central Storage directory has considerably more space (a quota of 500 GB)
and we therefore we will explain how to install Conda to your Central Store directory. This has the advantage of
only needing to be done once, and then it's usable on all of the Linux servers.

To install Conda on Central store, first ssh into one of the servers (see :ref:`Login`)::

    $ ssh telserv3
    telserv3>

By default, these servers run ``tcsh`` when you log in. The ``conda`` command works with ``bash``, not ``tcsh``,
so start a bash shell::

    [username@telserv3 ~]$ bash
    (base) bash$

Change directories to your Central Store directory::

    (base) bash$ cd /user/username

Because the servers at STScI run Linux, you will need to install a *Linux* version of Conda. Download and install
the latest Miniconda for linux using curl::

    (base) bash$ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    (base) bash$ bash Miniconda3-latest-Linux-x86_64.sh -p "/user/$(logname)/miniconda3"

The -p option changes the default installation path. Answer the installer prompts as before. When asked if you want
to run conda init say "yes"; this will automatically update your .bashrc and other necessary files.
When it's finished exit your current bash shell, then start a new one to automatically activate your new Conda installation.

Start a Bash shell to verify that your new Central Store Conda installation is availabl. Then check that you are using
the Conda version of Python::

    [username@telserv3 ~]$
    [username@telserv3 ~]$ bash
    (base) bash-4.1$ which python
    /user/username/miniconda3/bin/python

Now that you have the conda command available, we can install ``aplc_optimization`` on the server.

Installing ``aplc_optimization``
---------------------------------

Because of the storage limitations associated with the ``home`` diretory, we recommend installing ``aplc_optimization``
in your Central Store directory. Storing the toolkit on Central Store will also make your files available to your
personal machine (when connected to the STScI internal network).

If you haven't already, ``ssh`` into the server and open a
``bash`` shell::

    ssh telserv3
    [username@telserv3 ~]$ bash
    (base) bash$

Go to your user directory on Central Store and clone the ``aplc_optimization`` repository::

    (base) bash$ cd /user/$(logname)/
    (base) bash$ git clone -b develop https://github.com/spacetelescope/aplc_optimization.git

Create a new Conda environment on the Server using the ``environment.yml`` file::

    (base) bash$ cd aplc_optimization
    (base) bash$ conda env create --file environment.yml
    (base) bash$ source activate aplc_optimization
    (aplc_optimization) bash$

Installing Gurobi
----------------------
Again, due storage limitations, we will opt to install Gurobi in your Central Store directory. Whilst still in a ssh
session (see :ref:`Login`), go to your Central Store directory and download the latest 64-bit Linux version of the
Gurobi optimizer, using ``curl``::

    (aplc_optimization) bash$ cd /user/bnickson
    (aplc_optimization) bash$ curl -OL https://packages.gurobi.com/9.0/gurobi9.0.2_linux64.tar.gz

Open and extract the ``.tar.gz`` file::

    (aplc_optimization) bash$ tar xvzf gurobi9.0.2_linux64.tar.gz

This command will create a sub-directory ``gurobi912/linux64`` that contains the complete Linux Remote Services
distribution.

The Gurobi Optimizer makes use of several executable files. In order to allow these files to be found when needed,
you will have to modify your search path. Specifically, your PATH environment variable should be extended to include
``gurobi912/linux64/bin``. Add the following lines to your .bashrc file (which is located in the ``/home/username/`` directory)::

    export GUROBI_HOME="/user/username/gurobi912/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${GUROBI_HOME}/lib"

Unless LD_LIBRARY_PATH is already set, in which case add the following line instead::

    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

Obtain a new Gurobi license
'''''''''''''''''''''''''''''
Academic users are able to install and license Gurobi for their own use on more than one machine, however each
individual academic license can only be installed on a single physical machine. Consequently, in order to use the
optimizer on any of the Linux machines, a new individual academic license must be obtained (aside from the license
previously obtained for your personal machine).

In order to request a new license, Log into the Gurobi website and visit the
`Academic License page <https://www.gurobi.com/downloads/end-user-license-agreement-academic/>`_.
Once obtained, this new license will be visible in your
`Current Gurobi Licenses <https://www.gurobi.com/downloads/licenses/>`_ library.

Go to your `Current Gurobi Licenses <https://www.gurobi.com/downloads/licenses/>`_ library and open
the newly created license by clicking on the **Licence ID** number.

Retrieve and set up the Gurobi license on the server
''''''''''''''''''''''''''''''''''''''''''''''''''''''
The next step is to install the new Gurobi license on the Linux machine. After navigating to the **License Detail**
page on the Gurobi website and copy the ``grbgetkey`` command provided. During a ssh session (see :ref:`Login`), open a
``bash`` shell and run the ``grbgetkey`` command to to install the ``gurobi.lic`` file on the server.

Test the license on the server
''''''''''''''''''''''''''''''''
Once you have obtained a license key for a particular server, you can test the license using the Gurobi interactive shell.
Log out and back in to the server and run the ``gurobi.sh`` command::

    [username@telserv3 ~]$ quit
    Logging out...
    $ ssh telserv3
    [username@telserv3 ~]$ gurobi.sh
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


The Screen Utility
===================
When working with the remote servers, it is highly recommended to use the Linux screen utility, which lets you detach
from a running session, go away (shut down your computer, disconnect, etc.), and then reattach to it later.

To start a screen session, connects to a server with ``ssh`` and then start a screen session::

    $ ssh telserv3
    [username@telserv3 ~]$ screen

This will open a screen session, create a new window, and start a shell in that window.

It is also possible to create named sessions when you wish to run multiple screen sessions simultaneously.
To create a named session, run the screen command with the following arguments::

    [username@telserv3 ~]$ screen -S session_name

Where ``session_name`` is a descriptive session name of your choice.

In order to detach from a screen session, use the following key combination: Ctrl-a and then Ctrl-d.
You may then re-attach to the screen at anytime with::

    [username@telserv3 ~]$ screen -r session_name

If you would like to get a list of available screens, use::

    [username@telserv3 ~]$ screen -ls


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
----------------
Verifying network mounts
```````````````````````````
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
``````````````````````
The command to remotely connect to a given server is ``ssh <SERVERNAME>`` (where ``<SERVERNAME>`` is the server hostname,
e.g. ``ssh telserv3``). To enable X forwarding, you can use the ``-XY`` option (e.g. ``ssh -XY science1``). This allows
for displays that use X (Linux's default backend) to be forwarded to your monitor.

Go ahead and verify that you are able to log in to the TELSERV3 server, as follows::

    $ ssh telserv3

By default, the STScI servers run on ``tsch`` (not ``bash``) when you log in. To switch to a ``bash`` shell use::

    $ bash


Software Installation
```````````````````````

.. _Storage:

Storage
.................

When logged into a server (see :ref:`below <Login>`), any commands run in that SSH session will be executed on the server, not on your personal machine.
Thus, files and directories on your local machine will not be accessible when running the servers. Instead, all of the
servers have access to network-mounted home directories (i.e. ``/home/<username>``) and access to the ``/user`` space
on Central Store.

Home directories
================

Each user has a home directory (``/home/<username>``) that is mounted on each of the Linux machines. This home directory
is mounted and shared among all the machines, so that any file in that directory will be available on each machine.
However, space in this directory is limited (each user has a 10 GB quota) and installing software in the home directory
is usually discouraged.

``/user`` space
================
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
.................

It is recommended that the ``aplc-optimization`` package be installed using Conda (see :ref:`installing`).
As mentioned :ref:`above <Storage>`, files and directories on your personal machine cannot be seen by the servers. Therefore,
while you may have Conda installed locally, it won't be accessible to programs run on the one of the servers. Furthermore, the servers
run Linux, and will require a *Linux* compatable version of Miniconda/ Anaconda.

To use Conda for our computing needs, we will need to install Conda onto to one of the Linux servers---fortunately, we need
only go through the process on *one* server to make conda available on *all* servers. Because of the storage limitations
of Linux home diretory (see above), we will recommend installing conda to your Central Store directory instead.

First ssh into one of the servers::

    $ ssh telserv3
    telserv3>

Start a ``bash`` shell:

    telserv3> bash
    bash$

Change directories to your Central Store directory::

    bash$ cd "/user/$(logname)"



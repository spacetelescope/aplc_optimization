########################################################
aplc-optimization: an apodized pupil Lyot coronagraph design survey toolkit
########################################################

.. warning::

   Pardon our dust, this space is under construction.

``aplc-optimization`` is a Python software toolkit for exploring apodized pupil lyot coronagraph (APLC)
solutions for segmented telescope pupils. Its object-orientated approach simplifies
the interface for sampling large parameter spaces, and enables flexibility for implementing various mask architectures
and symmetry cases. The core module, example notebooks, and documentation are privately hosted at
github.com/spacetelescope/aplc-optimization.

.. figure:: ./HiCAT_1944_LS9_6pix_analysis.png
   :align: center
   :alt: Analysis results for HiCAT

   Figure 1: Analysis results from a HiCAT design study.

**Contributors:** ``aplc-optimization`` has been developed with the help of Emeil Por, Mamadou N'Daiye, Remi Flamary, Katherine St Laurent, Remi Soummer, Jamie Noss, Marshall Perrin, Bryony Nickson and Kelsey Glazer.

.. _getting_started:

Getting Started
================

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   introduction.rst
   workflow.rst
   server.rst
   tutorials.rst


-------------------------------

.. admonition:: How to cite ``aplc-optimization``

   In addition to this documentation, the ``aplc-optmization`` is described in the following references.  Users of ``aplc-optimization`` are encouraged to cite one of these.

    * St. Laurent et al. 2018, `"Apodized pupil Lyot coronagraphs designs for future segmented space telescopes", <https://doi.org/10.1117/12.2313902>`_, Proc. SPIE. 109682W
    * Zimmerman et al. 2016, `"Lyot coronagraph design study for large, segmented space telescope apertures” <https://doi.org/10.1117/12.2233205>`_,  Proc. SPIE. 99041Y

   If there is no appropriate place in the body of text to cite these proceedings, please include something along the
   lines of the following in your acknowledgements:


      *"This research made use if aplc-optimization, an object-orientated toolkit for performing Pupil Lyot coronagraph design surveys for segmented telescope apertures."*


----------------------------------------------

Detailed API Documentation
============================

.. toctree::
   :maxdepth: 1
   :caption: APLC-Optimization Package

   modules/inputs_generation


------------------------------------------------

Acknowledgements:
==================
- The `Space Telescope Science Institute collaborators <https://www.stsci.edu/>`_, in particular, the Segmented Coronagraph Design and Analysis (SCDA) team.
- The ``aplc-optimization`` *package was created in support of the Segmented Coronagraph Design and Analysis (SCDA) study,
funded by NASA's Exoplanet Exploration Program (ExEP). The goal of this study is to develop viable coronagraph
instrument concepts for a LUVOIR-type mission. The apodized pupil Lyot coronagraph (APLC) is one of several
coronagraph design families that the SCDA study is assessing.


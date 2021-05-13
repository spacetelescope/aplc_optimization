########################################################
aplc-optimization: an apodized pupil Lyot coronagraph design survey toolkit
########################################################

``aplc-optimization`` is a Python software toolkit for exploring apodized pupil lyot coronagraph (APLC)
solutions for abitrary telescope apertures. It's object-orientated approach simplifies
the interface for sampling large parameter spaces, and enables flexibility for implementing various mask architectures
and symmetry cases.

.. figure:: ./HiCAT_1944_LS9_6pix_analysis.png
   :align: center
   :alt: Analysis results for HiCAT

   Figure 1: Analysis results from a HiCAT design study.

The ``aplc-optimization`` toolkit was developed by the `Segmented Coronagraph Design & Analysis (SCDA) research team <https://www.stsci.edu/stsci-research/research-topics-and-programs/russell-b-makidon-optics-laboratory/meet-the-team>`_ at the Space Telescope Science Institute
(STScI) with the support of the `NASA Exoplanet Exploration Program (ExEP) <https://exoplanets.nasa.gov/exep/technology/SCDA/>`_ and is privately hosted at `github.com/spacetelescope/aplc-optimization <github.com/spacetelescope/aplc-optimization>`_.

.. _getting_started:

======================================
Getting started with aplc-optimization
======================================

.. toctree::
   :maxdepth: 1
   :caption: Contents

   introduction.rst
   installing.rst
   server.rst
   workflow.rst

.. toctree::
   :maxdepth: 1
   :caption: Advanced usage

   API Documentation



-------------------------------

.. admonition:: How to cite ``aplc-optimization``

   In addition to this documentation, the ``aplc-optmization`` toolkit is described in the following references.  Users of ``aplc-optimization`` are encouraged to cite one of these.

    * St. Laurent et al. 2018, `"Apodized pupil Lyot coronagraphs designs for future segmented space telescopes", <https://doi.org/10.1117/12.2313902>`_, Proc. SPIE. 109682W
    * Zimmerman et al. 2016, `"Lyot coronagraph design study for large, segmented space telescope apertures” <https://doi.org/10.1117/12.2233205>`_,  Proc. SPIE. 99041Y

   If there is no appropriate place in the body of text to cite these proceedings, please include something along the
   lines of the following in your acknowledgements:


      *"This research made use if aplc-optimization, an object-orientated toolkit for performing Pupil Lyot coronagraph design surveys for segmented telescope apertures."*



------------------------------------------------

Acknowledgements:
==================
- The `Space Telescope Science Institute collaborators <https://www.stsci.edu/>`_, in particular, the Segmented Coronagraph Design and Analysis (SCDA) team.
- The ``aplc-optimization`` *package was created in support of the Segmented Coronagraph Design and Analysis (SCDA) study, funded by NASA's Exoplanet Exploration Program (ExEP). The goal of this study is to develop viable coronagraph instrument concepts for a LUVOIR-type mission. The apodized pupil Lyot coronagraph (APLC) is one of several coronagraph design families that the SCDA study is assessing.


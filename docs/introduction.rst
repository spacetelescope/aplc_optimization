.. _introduction:

###################################################
Introduction
###################################################

The Apodized Pupil Lyot Coronagraph (APLC) is one of several coronagraph design families that is being assessed as part
of NASA's Exoplanet Exploration Program `Segmented aperture coronagraph design and analysis (SCDA) study <https://exoplanets.nasa.gov/exep/technology/SCDA/>`_, to further our understanding of the ability of coronagraphs to operate with segmented/ obscured apertures.
The APLCupil L is a Lyot-style coronagraph that suppresses starlight through a series of amplitude operations on the on-axis field.
A diagram of the APLC concept is shown in Figure 1. The mask architecture combines an entrance pupil apodization in plane A,
a downstream focal plane mask (FPM) in plane B, and the Lyot stop in the relayed pupil plane C to form the coronagraphic
image on a detector located in the re-imaged focal plane D. A perfect solution is when the transmission profile of the
apodizer in plane A is optimized so that the two scalar field components of the field in the Lyot plane approximately cancel.

.. figure:: ./APLC_schematic.png
   :align: center
   :width: 550
   :alt: schematic diagram for the APLC (Zimmerman 2016)

   **Figure 1**: The schematic layout of the APLC with an example for light propagation through the coronagraph. Label B-
   means just before plane B, and B+ means just after plane B. All images are on logarithmic scales.


Design Survey Strategy
=======================
APLCs are sensitive to several design parameters: telescope aperture shape and segmentation, central obscuration,
lyot stop shape and size, focal plane mask size, dark hole size, and bandwidth. In order to explore
this multi-dimensional parameter space, the ``aplc-optimization`` toolkit has been developed to simplify the organization,
execution and evaluation of large APLC design surveys.

For a given design survey, the toolkit generates a batch of linear programs to be executed on a computing cluster.
These linear programs rely on the third-party solver `Gurobi <Gurobi.com>`_ to determine the apodizer mask solution with maximum
off-axis transmission for a given set of design constraints (namely the raw contrast goal, dark zone extent and spectral bandwidth, telescope pupil, occulting
mask, IWA, OWA, and Lyot stop profile). Details about the apodizer optimization method,
including discrete algebraic models for the on-axis field propagation and definitions of the linear program objectives
and constraints, are given in the appendix of `Zimmerman et al. 2016 <https://www.spiedigitallibrary.org/conference-proceedings-of-spie/9904/1/Lyot-coronagraph-design-study-for-large-segmented-space-telescope-apertures/10.1117/12.2233205.full#c26>`_.


Getting aplc-optimization
=========================
For instructions on how to install the ``aplc-optimization`` toolkit, please see :ref:`installation`.

Running a design survey
========================
For instructions on how to run a design survey, please see :ref:`workflow`.

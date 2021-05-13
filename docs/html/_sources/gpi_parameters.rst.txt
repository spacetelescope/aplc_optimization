##############################
GPI APLC design parameters
##############################

* **Contrast**: 7
 (:ref:`Soummer et al. 2011 <Soummer11>`)


* **bandpass**: 20 *%*
 (:ref:`Soummer et al. 2011 <Soummer11>`)


* **ap_spid**: **False**
 (:ref:`Soummer et al. 2011 <Soummer11>`)
 Including the very thin Gemini Telescope secondary mirror support structures (1cm) would lead to a non-rotationally
 symmetric apodizer with additional system-level complexities that are not worth the potential gain, especially since
 the effect of these features can be mitigated in the Lyot plane (:ref:`Sivaramakrishnan & Lloyd 2005 <Siv05>`).
 Because the Gemini telescopes are optimized for the thermal infrared, the secondary supports are very thin
 compared to most telescopes, and have correspondingly little effect on the PSF. However we mask them out in the Lyot planes regardless.

Apodizer Details


* OWA:
 44 lambda_0/D

contrast is defined as the coronagraphic PSF normalized by the peak of an off-axis source that is not affected by the FPM.
This criterion accounts for the intensity reduction of an off-axis source due to the apodizer and Lyot stop


* Outer diameter undersizing foraction of 2% the diameter, to be compatavke with alignment tolerances.
* Optimal inner diameter oversizing is 98%

References
############

.. _Soummer11:
Soummer et al. 2011, `ApJ <https://iopscience.iop.org/article/10.1088/0004-637X/729/2/144>`_, `729 144 <https://ui.adsabs.harvard.edu/abs/2011ApJ...729..144S/abstract>`_

.. _Siv05:
Sivaramakrishnan, A., & Lloyd, J. P. 2005, `ApJ <https://iopscience.iop.org/article/10.1086/432903>`_, `633, 528 <https://ui.adsabs.harvard.edu/abs/2005ApJ...633..528S/abstract>`_




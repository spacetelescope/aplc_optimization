High-contrast imager for Complex Aperture Telescopes (HiCAT)
##############################################################

The HiCAT testbed is designed to study and develop high-contrast techniques for segmented telescopes
using wavefront control and coronagraphic starlight supression. The testbed design has the flexibility to enable studies
with increasing complexity for telescope aperture geometries starting with off-axis telescopes, then on-axis
telescopes with central obstruction and support structures - e.g. the Wide Field Infrared Survey Telescope (WFIRST) - up to on-axis segmented telescopes, including various concepts
for a Large UV, Optical, IR telescope (LUVOIR).

.. figure:: ./HiCAT.jpg
   :align: center
   :alt: View of the crucial HiCAT components

   View of the crucial HiCAT components, which consists of two continous deformable mirrors, a segmented deformable mirror, and an apodizer.

Hardware
###########

The testbed includes the following elements:

 - A segmented telescope simulator with a central obscuration and support structures using a
37-segment, aperture DM (IRIS AO). The central obscuration and support structures cal be added in the first pupil plane
using a laser-cut mask. This truly-segmented telescope simulater includes real co-phasing wavefront errors and control
ability, with an open-looop calibrated surface error of 9nm rms.
 - A custom Apodized Pupil Lyot Coronagraph (APLC) including a shaped pupil apodizer manufacted by
`Advanced Nanophotonics Inc. <https://www.spiedigitallibrary.org/conference-proceedings-of-spie/7761/77610F/Multiwalled-carbon-nanotubes-for-stray-light-suppression-in-space-flight/10.1117/12.864386.short>`_.
using carbon nanotubes coatings, a circular hard-edge focal plane mask re-using a mirror with an ultra-precise circular hole originally
from the Lyot project, which was first



HiCAT APLC
------------
HiCAT features a custom Apodized Pupil Lyot Coronagraph (APLC) optimized for the HiCAT aperture, which is similar to
one of the possible geometries considered for LUVOIR. The testbed also includes several external metrology
features for rapid replacement of parts, and in particular the ability to test multiple apodizers readily.

The APLC design is optimized for an aperture with geometric features slightly oversized/undersized to account for
alignment tolerancing. The FPM diameter is 8.543 $λ/D_apod$. The APLC dark zone is optimized for an IWA of
3.75 to 15 λ/Dapod. Based on the projected aperture on the deformable mirror, the maxinimum controllable frequency is
14.9 λ/Dapod. Note that the APLC IWA (3.75λ/Dapod) is smaller than the projected radius of the focal plane mask (4.3λ/Dapod),
this configuration provides enhanced robustness to tip/tilt, stellar diameter and low-order aberrations. The Lyot stop
outer diameter is undersized by a factor 0.83 from the apodizer circumscribed aperture, so that the final focal plane IWA and OWA are 3.11 λ/DLyot
and 12.45 λ/Dlyot.

The HiCAT APLC has been optimized to increase robustness to the Lyot stop alignment (0.1% alignment robustness).
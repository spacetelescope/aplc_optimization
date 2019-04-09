from hcipy import *
import matplotlib.pyplot as plt
import numpy as np

contrast = 1e-8
num_pix = 1024 #px
tau = 0.4
q_sci = 3 #px / (lambda_0/D)
iwa = 3.75 # lambda_0/D radius
owa = 15 # lambda_0/D radius
num_pix_foc = 50 # px diameter
foc_inner = 8.543 #lambda_0/D diameter
spectral_bandwidth = 0.1 # fractional
num_wavelengths = 3

pupil_grid = make_uniform_grid(num_pix, [1,1])
focal_grid = make_focal_grid(pupil_grid, 8, owa*1.2)

prop = FraunhoferPropagator(pupil_grid, focal_grid)

#aperture = read_fits('apodizers/HiCAT-N1024_NFOC0050_DZ0375_3000_C080_BW10_NLAM03_SHIFT10_01LS_ADAP4.fits')
aperture = read_fits('apodizers/HiCAT-N1024_NFOC0050_DZ0375_1500_C080_BW10_NLAM03_SHIFT20_05LS_ADAP4.fits')
#aperture2 = read_fits('masks/SYM-HiCAT-Aper_F-N0486_Hex3-Ctr0972-Obs0195-SpX0017-Gap0004.fits')
aperture = Field(aperture.ravel(), pupil_grid)
#aperture2 = Field(aperture2.ravel(), pupil_grid)

q_foc = num_pix_foc / foc_inner
x_foc = (np.arange(num_pix_foc) + 0.5 - num_pix_foc / 2) / q_foc
small_focal_grid = CartesianGrid(RegularCoords(1.0 / q_foc, [num_pix_foc, num_pix_foc], x_foc.min()))

focal_plane_mask = 1 - evaluate_supersampled(circular_aperture(foc_inner), small_focal_grid, 8)
focal_plane_mask2 = 1 - circular_aperture(foc_inner)(focal_grid)

lyot_stop = read_fits('masks/ehpor_lyot_mask_1024_bw.fits')
lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

#a = ((aperture * lyot_stop).sum() / (aperture2 * lyot_stop).sum())**2
#print(a)

coro = LyotCoronagraph(pupil_grid, focal_plane_mask, lyot_stop)
coro_without_lyot = LyotCoronagraph(pupil_grid, focal_plane_mask, None)

wf = Wavefront(aperture)
img = prop(coro(wf))
img_foc = prop(wf)
img_ref = prop(Wavefront(aperture * lyot_stop))
lyot = coro_without_lyot(wf)

plt.subplot(2,3,1)
imshow_field(aperture, cmap='gray')
plt.colorbar()
plt.subplot(2,3,2)
imshow_field(np.log10(img_foc.intensity / img_foc.intensity.max()), vmin=-5, vmax=0, cmap='hot')
plt.colorbar()
plt.subplot(2,3,3)
imshow_field(np.log10((img_foc.intensity * focal_plane_mask2 + 1e-20) / img_foc.intensity.max()), vmin=-5, vmax=0, cmap='hot')
plt.colorbar()
plt.subplot(2,3,4)
imshow_field(np.log10(lyot.intensity / lyot.intensity.max()), vmin=-3, vmax=0, cmap='hot')
plt.colorbar()
plt.subplot(2,3,5)
imshow_field(np.log10((lyot.intensity * lyot_stop + 1e-20) / lyot.intensity.max()), vmin=-3, vmax=0, cmap='hot')
plt.colorbar()
plt.subplot(2,3,6)
imshow_field(np.log10(img.intensity / img_ref.intensity.max()), vmin=-9, vmax=-5, cmap='hot')
plt.colorbar()
plt.show()

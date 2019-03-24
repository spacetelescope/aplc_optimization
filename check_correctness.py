from hcipy import *
import matplotlib.pyplot as plt
import numpy as np

fullres = read_fits('final_solution_aplc_fullres.fits')
refined = read_fits('final_solution_aplc_refined_1.fits')

print('Relative transmission between an interpolated and full optimization:', refined.sum() / fullres.sum())

contrast = 1e-8
num_pix = 486 #px
tau = 0.4
q_sci = 3 #px / (lambda_0/D)
iwa = 3.75 # lambda_0/D radius
owa = 15 # lambda_0/D radius
num_pix_foc = 50 # px diameter
foc_inner = 8.543 #lambda_0/D diameter
spectral_bandwidth = 0.1 # fractional
num_wavelengths = 3
testing = False

#num_pix_foc = 256

if testing:
	num_pix = 256
	owa = 8
	num_wavelengths = 2

pupil_grid = make_pupil_grid(num_pix)
focal_grid = make_focal_grid(pupil_grid, 8, owa*1.5)

prop = FraunhoferPropagator(pupil_grid, focal_grid)
if not testing:
    aperture = read_fits('final_solution_aplc_refined_1_hicat.fits')
    #aperture = read_fits('HiCAT-Apod_F-N0486_nImg0032_Hex3-Ctr0972-Obs0195-SpX0017-Gap0004_GreyFPM8543-M050_LS-Ann-gy-ID0345-OD0807-SpX0036_DZ-C030-080-Sep037-150_Bw10-Lam3_shiftXY050.fits')
    aperture = Field(aperture.ravel(), pupil_grid)
else:
    aperture = read_fits('')

small_focal_grid = make_pupil_grid(num_pix_foc, foc_inner)
focal_plane_mask = 1 - circular_aperture(foc_inner)(small_focal_grid)
focal_plane_mask2 = 1 - circular_aperture(foc_inner)(focal_grid)

if testing:
	lyot_stop = make_obstructed_circular_aperture(0.95, 0.3, 3, 0.02)
	lyot_stop = evaluate_supersampled(lyot_stop, pupil_grid, 4)
else:
	lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
	lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

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
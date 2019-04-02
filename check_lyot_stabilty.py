from hcipy import *
import matplotlib.pyplot as plt
import numpy as np

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

pupil_grid = make_uniform_grid(num_pix, [1,1])
focal_grid = make_focal_grid(pupil_grid, 8, owa*1.2)

prop = FraunhoferPropagator(pupil_grid, focal_grid)

aperture = read_fits('full_resolution_lyot_robust_apod.fits')
aperture = read_fits('HiCAT-Apod_F-N0486_nImg0032_Hex3-Ctr0972-Obs0195-SpX0017-Gap0004_GreyFPM8543-M050_LS-Ann-gy-ID0345-OD0807-SpX0036_DZ-C030-080-Sep037-150_Bw10-Lam3_shiftXY050.fits')

#aperture = np.repeat(np.repeat(aperture, 2, 0), 2, 1)
aperture = Field(aperture.ravel(), pupil_grid)

q_foc = num_pix_foc / foc_inner
x_foc = (np.arange(num_pix_foc) + 0.5 - num_pix_foc / 2.0) / q_foc
small_focal_grid = CartesianGrid(SeparatedCoords((x_foc, x_foc)))
focal_plane_mask = 1 - circular_aperture(foc_inner)(small_focal_grid)
focal_plane_mask2 = 1 - circular_aperture(foc_inner)(focal_grid)

lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-gy-ID0345-OD0807-SpX0036_shiftX+000.fits')
#lyot_stop = np.repeat(np.repeat(lyot_stop, 2, 0), 2, 1)
lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

dark_zone = (circular_aperture(owa * 2)(focal_grid) - circular_aperture(iwa * 2)(focal_grid)).astype('bool')

coro = LyotCoronagraph(pupil_grid, focal_plane_mask, lyot_stop)
coro_without_lyot = LyotCoronagraph(pupil_grid, focal_plane_mask, None)

##########################################################################
# Centered Lyot stop

def show_summary(aperture, coro, coro_without_lyot, focal_plane_mask2):
    wf = Wavefront(aperture)
    img = prop(coro(wf))
    img_foc = prop(wf)
    img_ref = prop(Wavefront(aperture * lyot_stop))
    lyot = coro_without_lyot(wf)

    plt.figure()
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

show_summary(aperture, coro, coro_without_lyot, focal_plane_mask2)

#############################################################################
# Shifted Lyot stop

dxs = np.array([-2, -1, 0, 1, 2])

dither_grid = CartesianGrid(SeparatedCoords((dxs, dxs)))

E = []
mean_intensity = []
plt.figure(figsize=(8,8))
for i, (dx, dy) in enumerate(dither_grid.points):
    shifted_lyot_stop = np.roll(np.roll(lyot_stop.shaped, dx, 1), dy, 0).ravel()
    coro = LyotCoronagraph(pupil_grid, focal_plane_mask, shifted_lyot_stop)

    img = prop(coro(Wavefront(aperture)))
    img_ref = prop(Wavefront(aperture * lyot_stop))

    x = i % len(dxs)
    y = i // len(dxs)
    print(x + (len(dxs) - y - 1) * len(dxs) + 1)

    plt.subplot(len(dxs), len(dxs), x + (len(dxs) - y - 1) * len(dxs) + 1)
    imshow_field(np.log10(img.intensity / img_ref.intensity.max()), vmin=-8.5, vmax=-5, cmap='hot')
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])

    E.append(img.electric_field.at((4.5,0)))
    mean_intensity.append(img.intensity[dark_zone].mean() / img_ref.intensity.max())

plt.subplots_adjust(top=0.98, bottom=0.02, left=0.02, right=0.98, wspace=0, hspace=0)

plt.figure()
imshow_field(np.log10(mean_intensity), dither_grid, cmap='inferno')
plt.colorbar()

E = np.array(E).reshape((len(dxs), len(dxs)))

plt.figure()
plt.plot(E.real)
plt.plot(E.imag, ls='--')
plt.show()
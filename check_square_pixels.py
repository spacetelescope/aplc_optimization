from hcipy import *
import matplotlib.pyplot as plt
import numpy as np

class LyotCoronagraphSquarePixels(LyotCoronagraph):
	def __init__(self, input_grid, focal_plane_mask, lyot_stop=None):
		LyotCoronagraph.__init__(self, input_grid, focal_plane_mask, lyot_stop)
	
	def forward(self, wavefront):
		wf_foc = self.prop.forward(wavefront)
		TF = np.sinc(wf_foc.electric_field.grid.x / wavefront.wavelength * wavefront.electric_field.grid.delta[0])
		TF *= np.sinc(wf_foc.electric_field.grid.y / wavefront.wavelength * wavefront.electric_field.grid.delta[1])
		print(TF.min())
		wf_foc.electric_field -= self.focal_plane_mask.forward(wf_foc).electric_field
		wf_foc.electric_field *= TF

		lyot = self.prop.backward(wf_foc)
		lyot.electric_field[:] = wavefront.electric_field - lyot.electric_field

		if self.lyot_stop is not None:
			lyot = self.lyot_stop.forward(lyot)
		
		return lyot
	
	def backward(self, wavefront):
		raise NotImplementedError()

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
pupil_oversamplings = np.arange(20)+1
dx = 0
dy = 0

imgs = []
for pupil_oversampling in pupil_oversamplings:
	print(pupil_oversampling)
	pupil_grid = make_uniform_grid(num_pix * pupil_oversampling, [1,1])
	focal_grid = make_focal_grid(pupil_grid, 8, owa*1.2)

	prop = FraunhoferPropagator(pupil_grid, focal_grid)

	aperture = read_fits('full_resolution_lyot_robust_apod.fits')
	print(aperture.shape)
	aperture = np.repeat(np.repeat(aperture, pupil_oversampling, 0), pupil_oversampling, 1)
	aperture = Field(aperture.ravel(), pupil_grid)

	small_focal_grid = make_pupil_grid(num_pix_foc, foc_inner)
	focal_plane_mask = 1 - circular_aperture(foc_inner)(small_focal_grid)
	focal_plane_mask2 = 1 - circular_aperture(foc_inner)(focal_grid)

	lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
	lyot_stop = np.roll(np.roll(lyot_stop, dx, 1), dy, 0)
	lyot_stop = np.repeat(np.repeat(lyot_stop, pupil_oversampling, 0), pupil_oversampling, 1)
	lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

	coro = LyotCoronagraph(pupil_grid, focal_plane_mask, lyot_stop)

	wf = Wavefront(aperture)
	img = prop(coro(wf)).intensity
	img_ref = prop(Wavefront(aperture * lyot_stop)).intensity
	imgs.append(img / img_ref.max())

# With square pixels Lyot coronagraph
pupil_grid = make_uniform_grid(num_pix, [1,1])
focal_grid = make_focal_grid(pupil_grid, 8, owa*1.2)

prop = FraunhoferPropagator(pupil_grid, focal_grid)

aperture = read_fits('full_resolution_lyot_robust_apod.fits')
aperture = Field(aperture.ravel(), pupil_grid)

small_focal_grid = make_pupil_grid(num_pix_foc, foc_inner)
focal_plane_mask = 1 - circular_aperture(foc_inner)(small_focal_grid)
focal_plane_mask2 = 1 - circular_aperture(foc_inner)(focal_grid)

lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
lyot_stop = np.roll(np.roll(lyot_stop, dx, 1), dy, 0)
lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

coro = LyotCoronagraphSquarePixels(pupil_grid, focal_plane_mask, lyot_stop)

wf = Wavefront(aperture)
img = prop(coro(wf)).intensity
img_ref = prop(Wavefront(aperture * lyot_stop)).intensity
imgs.append(img / img_ref.max())

# Analyze
dark_zone = (circular_aperture(owa * 2)(focal_grid) - circular_aperture(iwa * 2)(focal_grid)).astype('bool')

img_ref = imgs[-2]
max_diff = [np.max((img - img_ref)[dark_zone]) for img in imgs]
mean_diff = [np.mean((img - img_ref)[dark_zone]) for img in imgs]

plt.figure(figsize=(5,3))
plt.plot(pupil_oversamplings[:-1], mean_diff[:-2], '-o', c=colors.blue, label='Mean difference in dark zone')
plt.plot(pupil_oversamplings[:-1], max_diff[:-2], '-v', c=colors.blue, label='Maximum difference in dark zone')
plt.plot([1], [mean_diff[-1]], 'o', c=colors.red, label='Sinc approximate mean difference')
plt.plot([1], [max_diff[-1]], 's', c=colors.red, label='Sinc approximate max difference')
plt.xlabel('Number of pixels per superpixel')
plt.ylabel('Difference in raw contrast')
plt.legend()
plt.yscale('log')
plt.show()

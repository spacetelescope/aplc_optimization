from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from scipy.ndimage.morphology import grey_erosion, grey_dilation

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

if testing:
	num_pix = 128
	owa = 8
	num_wavelengths = 2

pupil_grid = make_pupil_grid(num_pix)
focal_grid = make_focal_grid(pupil_grid, q_sci, owa * (1 + spectral_bandwidth))
focal_mask = (circular_aperture(owa*2)(focal_grid) - circular_aperture(iwa*2)(focal_grid)).astype('bool')

fraun_prop = FraunhoferPropagator(pupil_grid, focal_grid)

if testing:
	aperture = circular_aperture(1)(pupil_grid)
else:
	aperture = read_fits('masks/SYM-HiCAT-Aper_F-N0486_Hex3-Ctr0972-Obs0195-SpX0017-Gap0004.fits')
	aperture = Field(aperture.ravel(), pupil_grid)

small_focal_grid = make_pupil_grid(num_pix_foc, foc_inner)
focal_plane_mask = 1 - circular_aperture(foc_inner)(small_focal_grid)

if testing:
	lyot_stop = circular_aperture(0.95)(pupil_grid)
else:
	lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
	lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

coro = LyotCoronagraph(pupil_grid, focal_plane_mask, lyot_stop)
coro_prop = OpticalSystem([coro, fraun_prop])

def optimize_at_scale(pupil_grid, focal_grid, focal_mask, prop, aperture, lyot_stop, last_optim, subsampling, contrast, wavelengths):
	'''
	last_optim=None means no previous knowledge is known; all pixels will be optimized
	'''
	inds = np.arange(pupil_grid.size).reshape((pupil_grid.shape[1]//subsampling, subsampling, pupil_grid.shape[0]//subsampling, subsampling))
	inds = np.swapaxes(inds, 1, 2).reshape((pupil_grid.size//(subsampling**2), -1))

	aperture_subsampled = subsample_field(aperture, subsampling)
	lyot_stop_subsampled = subsample_field(lyot_stop, subsampling)

	blind = last_optim is None

	if blind:
		last_optim = Field(np.ones(aperture_subsampled.grid.size) * 0.5, aperture_subsampled.grid)
	else:
		last_optim = np.repeat(np.repeat(last_optim.shaped, 2, 1), 2, 0).ravel()
		last_optim = Field(last_optim, aperture_subsampled.grid)

	structure = np.array([[0,1,0],[1,1,1],[0,1,0]])
	structure = np.array([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]])
	structure = np.array([[0,1,1,0],[1,1,1,1],[1,1,1,1],[0,1,1,0]])

	pixels_to_optimize = (grey_dilation(last_optim.shaped, structure=structure) - grey_erosion(last_optim.shaped, structure=structure)).ravel() - 2
	pixels_to_optimize = Field(pixels_to_optimize, last_optim.grid)

	#imshow_field(pixels_to_optimize)
	#plt.show()

	pixels_to_optimize = np.abs(pixels_to_optimize) > 1e-3
	pixels_to_optimize2 = np.logical_and(last_optim < (1 - 1e-3), last_optim > 1e-3)

	pixels_to_optimize = np.logical_or(pixels_to_optimize, pixels_to_optimize2)

	if not blind:
		imshow_field(pixels_to_optimize, cmap='Reds', mask=aperture_subsampled, vmin=0, vmax=1.5)
		plt.savefig('pixels_to_optimize%d.pdf' % subsampling)
		plt.clf()

	x0 = Field(np.zeros(pupil_grid.size), pupil_grid)

	print('Calculating mask for %d/%d variables.' % (np.sum(pixels_to_optimize * (aperture_subsampled > 0)), np.sum(aperture_subsampled > 0)))
	print('Number of constraints: %d' % (np.sum(focal_mask) * 2))
	print('Number of wavelengths: %d' % len(wavelengths))

	# Create model
	model = gp.Model('lp')
	model.Params.Threads = 0
	model.Params.Crossover = 0
	model.Params.Method = 2

	optimize_mask = np.logical_and(pixels_to_optimize, aperture_subsampled > 0)
	n = np.sum(optimize_mask)
	n = int(n)
	x_vars = model.addVars(n, lb=0, ub=1)

	M_max = (aperture_subsampled * lyot_stop_subsampled)[optimize_mask]
	obj = gp.quicksum((x_vars[i] * M_max[i] for i in range(n)))
	model.setObjective(obj, gp.GRB.MAXIMIZE)

	for i, wavelength in enumerate(wavelengths):
		j = 0
		mat = []
		base_electric_field = np.zeros(focal_grid.size, dtype='complex')

		for ind, amp, to_optimize in zip(inds, last_optim, pixels_to_optimize):
			if np.sum(aperture[ind]) < 1e-3:
				continue
			else:
				x = Field(np.zeros(pupil_grid.size), pupil_grid)
				x[ind] = aperture[ind]
				
				if not to_optimize:
					x0 += x * amp
				else:
					y = prop(Wavefront(x, wavelength)).electric_field[focal_mask]
					mat.append(y)
					j += 1
					if j % 1000 == 0:
						print('Wavelength: %d/%d; Variables: %d/%d' % (i+1, len(wavelengths), j, n))

		base_electric_field = prop(Wavefront(x0, wavelength)).electric_field[focal_mask]

		M = np.array(mat).T
		M = np.vstack([M.real, M.imag])

		E_0 = np.concatenate((base_electric_field.real, base_electric_field.imag))

		m, n = M.shape
		print(m, n)

		contrast_requirement = np.ones(m) * np.sqrt(contrast)

		for i, ee in enumerate(M):
			e = gp.quicksum((x_vars[i] * ee[i] for i in range(n)))
			model.addConstr(e <= -E_0[i] + contrast_requirement[i])
			model.addConstr(e >= -E_0[i] - contrast_requirement[i])
		
		del M
	
	model.optimize()

	solution = np.array([x_vars[i].x for i in range(n)])

	sol = last_optim
	sol[optimize_mask] = solution
	sol = Field(sol, make_subsampled_grid(pupil_grid, subsampling))
	imshow_field(sol, cmap='gray', mask=aperture_subsampled)
	plt.savefig('aplc_%d.pdf' % subsampling)
	plt.clf()

	return sol

last_optim = None
if num_wavelengths > 1:
	wavelengths = np.linspace(-spectral_bandwidth / 2, spectral_bandwidth / 2, num_wavelengths) + 1
else:
	wavelengths = [1]

for subsampling in [2,1]:
	last_optim = optimize_at_scale(pupil_grid, focal_grid, focal_mask, coro_prop, aperture, lyot_stop, last_optim, subsampling, contrast*tau, wavelengths)
	aperture_subsampled = subsample_field(aperture, subsampling)
	write_fits(last_optim * aperture_subsampled, 'final_solution_aplc_refined_%d.fits' % subsampling)

optim = optimize_at_scale(pupil_grid, focal_grid, focal_mask, coro_prop, aperture, lyot_stop, None, subsampling, contrast*tau, wavelengths)
write_fits(optim * aperture, 'final_solution_aplc_fullres.fits')
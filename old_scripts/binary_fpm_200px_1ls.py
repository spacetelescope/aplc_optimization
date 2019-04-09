run_without_gurobipy = False

import matplotlib as mpl
mpl.use('Agg')
from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
if not run_without_gurobipy:
	import gurobipy as gp
from scipy.ndimage.morphology import grey_erosion, grey_dilation

def calculate_pixels_to_optimize(last_optim, pupil_subsampled):
	"""Calculate the pixels to be used for the optimization for the adaptive algorithm.

	The exact selection of the pixels to be taken into account depends on future research.

	Parameters
	----------
	last_optim : Field
		The previously-optimized apodizer.
	pupil_subsampled : Field
		The telescope pupil subsampled to the same resolution as `last_optim`.

	Returns
	-------
	Field
		A boolean mask to indicate which pixels to take into account.
	"""
	if last_optim is None:
		return pupil_subsampled > 0

	structure = np.array([[0,1,0],[1,1,1],[0,1,0]])
	#structure = np.array([[0,1,1,0],[1,1,1,1],[1,1,1,1],[0,1,1,0]])
	#structure = np.array([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]])

	a = (grey_dilation(last_optim.shaped, structure=structure) - grey_erosion(last_optim.shaped, structure=structure)).ravel() - 2
	a = np.abs(a) > 1e-3

	b = np.logical_and(last_optim < (1 - 1e-3), last_optim > 1e-3)

	c = np.logical_or(a, b)

	return np.logical_and(c, pupil_subsampled > 0)

def optimize_aplc(pupil, focal_plane_mask, lyot_stops, dark_zone_mask, wavelengths, contrast, num_scalings=1,
	force_no_x_symmetry=False, force_no_y_symmetry=False, do_apodizer_throughput_maximization=False):
	"""Optimize an APLC with a (optional) iterative algorithm taking into account x and/or y mirror symemtries.

	Parameters
	----------
	pupil : Field
		The telescope pupil
	focal_plane_mask : Field
		The focal plane mask of the APLC. The grid is assumed to be in lambda_0/D.
	lyot_stops : list of Fields
		A list of the Lyot stops used in the optimization problem.
	dark_zone_mask : Field
		A binary field indicating which pixels belong to the dark zone.
	wavelengths : array_like
		An array of wavelengths as fractions of lambda_0.
	contrast : scalar
		The contrast that needs to be achieved in the dark zone.
	num_scalings : int
		The number of scalings that need to be performed in the adaptive algorithm. A value of 1 effectively turns
		the adaptive algorithm off.
	force_no_x_symmetry : boolean
		Force the algorithm to ignore any x mirror symmetry that might exist in the problem.
	force_no_y_symmetry : boolean
		Force the algorithm to ignore any y mirror symmetry that might exist in the problem.

	Returns
	-------
	Field
		The optimized apodizer for the APLC. This has not been multiplied by the telescope pupil.
	"""
	pupil_grid = pupil.grid
	focal_grid = dark_zone_mask.grid

	aplc = LyotCoronagraph(pupil_grid, focal_plane_mask)
	prop = FraunhoferPropagator(pupil_grid, focal_grid)

	focal_grid_0 = CartesianGrid(UnstructuredCoords([np.array([0.0]), np.array([0.0])]), np.array([1.0]))
	prop_0 = FraunhoferPropagator(pupil_grid, focal_grid_0)

	prior = None
	subsamplings = 2**np.arange(num_scalings)[::-1]

	# Determine pupil symmetries
	x_symm_pupil = np.allclose(pupil.shaped[:,::-1], pupil.shaped)
	y_symm_pupil = np.allclose(pupil.shaped[::-1,:], pupil.shaped)

	print('Telescope pupil:')
	print('   Mirror symmetry in x: %s' % ('yes' if x_symm_pupil else 'no'))
	print('   Mirror symmetry in y: %s' % ('yes' if y_symm_pupil else 'no'))
	print('')

	# Determine lyot stop symmetries
	x_symm_lyot_stops = [np.allclose(lyot_stop.shaped[:,::-1], lyot_stop.shaped) for lyot_stop in lyot_stops]
	y_symm_lyot_stops = [np.allclose(lyot_stop.shaped[::-1,:], lyot_stop.shaped) for lyot_stop in lyot_stops]

	lyot_stop_duplication = np.zeros(len(lyot_stops), dtype='bool')
	lyot_stop_duplication_reason = [[] for i in range(len(lyot_stops))]

	x_symm = x_symm_pupil
	y_symm = y_symm_pupil

	lyot_stop_status = []
	for i, a in enumerate(lyot_stops):
		print('Lyot stop #%d:' % i)

		if lyot_stop_duplication[i]:
			print('   Will be ignored due to symmetries with Lyot stops #' + str(lyot_stop_duplication_reason[i]))
			continue

		print('   Mirror symmetry in x: %s' % ('yes' if x_symm_lyot_stops[i] else 'no'))
		print('   Mirror symmetry in y: %s' % ('yes' if y_symm_lyot_stops[i] else 'no'))

		if x_symm_pupil and not x_symm_lyot_stops[i]:
			print('   Searching for mirror symmetric Lyot stops in x...')
			for j, b in enumerate(lyot_stops):
				if j <= i:
					continue

				if np.allclose(a.shaped[:,::-1], b.shaped):
					print('      Found Lyot stop #%d to fit.' % j)
					lyot_stop_duplication[j] = True
					lyot_stop_duplication_reason[j].append(i)
					x_symm = True
					break
			else:
				print('      No Lyot stop found with this symmetry. This breaks the mirror symmetry in x of the optimization.')
				x_symm = False

		if y_symm_pupil and not y_symm_lyot_stops[i]:
			print('   Searching for mirror symmetric Lyot stops in y...')
			for j, b in enumerate(lyot_stops):
				if j <= i:
					continue

				if np.allclose(a.shaped[::-1,:], b.shaped):
					print('      Found Lyot stop #%d to fit.' % j)
					lyot_stop_duplication[j] = True
					lyot_stop_duplication_reason[j].append(i)
					y_symm = True
					break
			else:
				print('      No Lyot stop found with this symmetry. This breaks the mirror symmetry in y of the optimization.')
				y_symm = False
	print('')

	print('Complete APLC:')
	print('   Mirror symmetry in x: %s' % ('yes' if x_symm else 'no'))
	print('   Mirror symmetry in y: %s' % ('yes' if y_symm else 'no'))

	if force_no_x_symmetry:
		print('   The user forced me to ignore mirror-symmetries in x.')
		x_symm = False
	if force_no_y_symmetry:
		print('   The user forced me to ignore mirror-symmetries in y.')
		y_symm = False
	print('')

	# Find number of constraints per focal-plane point
	num_constraints_per_focal_point = []
	for i, lyot_stop in enumerate(lyot_stops):
		if lyot_stop_duplication[i]:
			# Lyot stop already taken into account due to symmetries
			num_constraints_per_focal_point.append(0)
		elif (x_symm and x_symm_lyot_stops[i]) and (y_symm and y_symm_lyot_stops[i]):
			# Only constrain real part
			num_constraints_per_focal_point.append(1)
		else:
			# Constrain both real and imag part
			num_constraints_per_focal_point.append(2)

	# We are optimizing amplitude (=real part) only
	dark_zone_mask *= focal_grid.x > 0

	if x_symm or y_symm:
		dark_zone_mask *= focal_grid.y > 0

	dark_zone_mask = dark_zone_mask.astype('bool')
	mm = int(np.sum(dark_zone_mask))

	# Calculate number of constraints per wavelength
	m = int(np.sum(num_constraints_per_focal_point)) * mm

	# Iterate from lowest to highest resolution
	for subsampling in subsamplings:
		# Calculated subsampled pupil and lyot stops
		pupil_subsampled = subsample_field(pupil, subsampling)
		lyot_stops_subsampled = [subsample_field(lyot_stop, subsampling) for lyot_stop in lyot_stops]
		pupil_grid_subsampled = pupil_subsampled.grid

		# Calculate which pixels belong to which superpixel
		inds = np.arange(pupil_grid.size).reshape((pupil_grid.shape[1]//subsampling, subsampling, pupil_grid.shape[0]//subsampling, subsampling))
		inds = np.swapaxes(inds, 1, 2).reshape((pupil_grid_subsampled.shape[0], pupil_grid_subsampled.shape[1], -1))#.reshape((pupil_grid.size//(subsampling**2), -1))

		# Apply x,y-mirror-symmetries
		mask = Field(np.ones(pupil_grid_subsampled.size), pupil_grid_subsampled)
		if x_symm:
			inds = np.concatenate((inds, inds[:,::-1,:]), axis=2)
			mask *= pupil_grid_subsampled.x < 0
		if y_symm:
			inds = np.concatenate((inds, inds[::-1,:,:]), axis=2)
			mask *= pupil_grid_subsampled.y < 0

		mask = mask.astype('bool')
		inds = inds.reshape((pupil_grid_subsampled.size, -1))
		#inds = [inds[i,:] for i in range(pupil_grid_subsampled.size) if mask[i]]

		# Upscale last optim to current resolution
		blind = prior is None
		if blind:
			# No prior information; assume totally dark apodizer
			last_optim = Field(np.zeros(pupil_grid_subsampled.size), pupil_grid_subsampled)
			prior = Field(np.zeros(pupil_grid.size), pupil_grid)
		else:
			# Upscale prior information by factor 2
			#last_optim = np.repeat(np.repeat(last_optim.shaped, 2, 1), 2, 0).ravel()
			last_optim = subsample_field(prior, subsampling)
			#last_optim = Field(last_optim, pupil_grid_subsampled)

		# Get pixels to optimize
		optimize_mask = np.logical_and(calculate_pixels_to_optimize(last_optim, pupil_subsampled), mask)
		if blind:
			optimize_mask[:] = np.logical_and(pupil_subsampled > 0, mask)
		n = int(np.sum(optimize_mask*mask))

		print('Starting optimization at scale %d with %d variables and %d constraints.' % (subsampling, n, m*len(wavelengths)))
		print('Creating model...')

		# Create Gurobi model
		if not run_without_gurobipy:
			model = gp.Model('lp')
			model.Params.Threads = 0
			model.Params.Crossover = 0
			model.Params.Method = 2
			x_vars = model.addVars(n, lb=0, ub=1)

		print('Calculating and adding constraints...')

		# Create problem matrix for one wavelength but for all Lyot stops
		M = np.empty((m, n))

		# Add constraints for each wavelength
		for wl_i, wavelength in enumerate(wavelengths):
			j = 0
			x0 = Field(np.zeros(pupil_grid.size), pupil_grid)
			x = Field(np.zeros(pupil_grid.size), pupil_grid)

			# Calculate norm electric field for each Lyot stop
			norms = []
			for lyot_stop in lyot_stops:
				norms.append(prop_0(Wavefront(pupil * lyot_stop, wavelength)).electric_field[0])

			for ind, amp, to_optimize, masked in zip(inds, last_optim, optimize_mask, mask):
				x[:] = 0
				x[ind] = pupil[ind]

				if not to_optimize:
					# Do not optimize this pixel
					# Add to accumulator pupil-plane wavefront
					x0 += x * amp
				else:
					# Calculate field before the Lyot stop
					lyot = aplc(Wavefront(x, wavelength))

					k = 0
					for i, lyot_stop in enumerate(lyot_stops):
						if num_constraints_per_focal_point[i] == 0:
							continue

						# Apply the Lyot stop and get focal-plane electric field
						lyot_copy = lyot.copy()
						lyot_copy.electric_field *= lyot_stop
						E = prop(lyot_copy).electric_field[dark_zone_mask]
						E /= norms[i]

						# Add only imaginary or both imaginary and real constraints depending on symmetry
						M[k*mm:(k+1)*mm,j] = E.real
						if num_constraints_per_focal_point[i] == 2:
							M[(k+1)*mm:(k+2)*mm,j] = E.imag

						k += num_constraints_per_focal_point[i]
					j += 1

					if j % 1000 == 0:
						print('Wavelength %d/%d; Variable %d/%d' % (wl_i + 1, len(wavelengths), j, n))

			# Calculate base electric field
			base_electric_field = []
			for i, lyot_stop in enumerate(lyot_stops):
				if num_constraints_per_focal_point[i] == 0:
					continue

				wf = aplc(Wavefront(x0, wavelength))
				wf.electric_field *= lyot_stop
				img = prop(wf)

				base_electric_field.append(img.electric_field.real)
				if num_constraints_per_focal_point[i] == 2:
					base_electric_field.append(img.electric_field.imag)
			base_electric_field = np.concatenate(base_electric_field)

			# Calculate contrast requirement
			contrast_requirement = np.tile((np.ones(dark_zone_mask.size) * np.sqrt(contrast))[dark_zone_mask], int(np.sum(num_constraints_per_focal_point)))

			# Add constraints
			for ee, e0, c0 in zip(M, base_electric_field, contrast_requirement):
				e = gp.LinExpr(ee, x_vars.values())
				model.addConstr(e <= (-e0 + c0))
				model.addConstr(e >= (-e0 - c0))

		del M

		# Use central Lyot stop for throughput metric (assume that this is the unshifted Lyot stop)
		if do_apodizer_throughput_maximization:
			M_max = pupil_subsampled[optimize_mask]
		else:
			M_max = (pupil_subsampled * lyot_stops_subsampled[0])[optimize_mask]
		obj = gp.LinExpr(M_max, x_vars.values())
		model.setObjective(obj, gp.GRB.MAXIMIZE)

		# Optimize model
		print('Start optimization...')
		model.optimize()
		print('Optimization finished!')

		# Extract solution from Gurobi
		solution = np.array([x_vars[i].x for i in range(n)])

		# Integrate solution into upsampled old solution
		sol = prior
		j = 0
		for ind, to_optimize in zip(inds, optimize_mask):
			if to_optimize:
				sol[ind] = solution[j]
				j += 1
		sol = Field(sol, pupil_grid)
		prior = sol

	return prior

if __name__ == '__main__':
	contrast = 1e-8
	num_pix = 486
	q_sci = 2 # px / (lambda_0/D)
	iwa = 3.75 # lambda_0/D
	owa = 15 # lambda_0/D
	n_foc = 200 # px diameter
	foc_inner = 8.543 # lambda_0/D diameter
	spectral_bandwidth = 0.1 # fractional
	num_wavelengths = 3
	lyot_stop_robustness = False
	lyot_stop_shift = 1 # px
	tau = 0.55 # expected planet peak intensity (relative to without focal plane mask)

	fname_pupil = 'masks/ehpor_apodizer_mask_486_gy.fits'
	fname_lyot_stop = 'masks/ehpor_lyot_mask_486_gy.fits'

	fname = 'apodizers/HiCAT-N%04d_NFOC%04d_DZ%04d_%04d_C%03d_BW%02d_NLAM%02d_SHIFT%02d_BWFPM' % (num_pix, n_foc, iwa*100, owa*100, -10*np.log10(contrast), spectral_bandwidth*100, num_wavelengths, lyot_stop_shift*10)
	print('Apodizer will be saved to:')
	print('   ' + fname + '.fits')
	print('')

	pupil_grid = make_uniform_grid((num_pix, num_pix), 1)

	pupil = read_fits(fname_pupil)
	pupil = Field(pupil.ravel(), pupil_grid)

	lyot_stop = read_fits(fname_lyot_stop)
	lyot_stop = Field(lyot_stop.ravel(), pupil_grid)
	lyot_stops = [lyot_stop]

	if lyot_stop_robustness:
		lyot_stop_pos_x = np.roll(lyot_stop.shaped, lyot_stop_shift, 1).ravel()
		lyot_stop_neg_x = np.roll(lyot_stop.shaped, -lyot_stop_shift, 1).ravel()
		lyot_stop_pos_y = np.roll(lyot_stop.shaped, lyot_stop_shift, 0).ravel()
		lyot_stop_neg_y = np.roll(lyot_stop.shaped, -lyot_stop_shift, 0).ravel()

		lyot_stops.extend([lyot_stop_pos_x, lyot_stop_neg_x, lyot_stop_pos_y, lyot_stop_neg_y])

	n_sci = int((np.ceil(owa) + 1) * q_sci) * 2
	x_sci = (np.arange(n_sci) + 0.5 - n_sci / 2) / q_sci
	focal_grid = CartesianGrid(SeparatedCoords((x_sci, x_sci)))

	dark_zone_mask = circular_aperture(owa * 2)(focal_grid) - circular_aperture(iwa * 2)(focal_grid)

	q_foc = n_foc / foc_inner
	x_foc = (np.arange(n_foc) + 0.5 - n_foc / 2) / q_foc
	focal_mask_grid = CartesianGrid(SeparatedCoords((x_foc, x_foc)))

	focal_plane_mask = 1 - circular_aperture(foc_inner)(focal_mask_grid)

	if num_wavelengths == 1:
		wavelengths = [1]
	else:
		wavelengths = np.linspace(-spectral_bandwidth / 2, spectral_bandwidth / 2, num_wavelengths) + 1
	
	apodizer = optimize_aplc(pupil, focal_plane_mask, lyot_stops, dark_zone_mask, wavelengths, contrast * tau, num_scalings=1, force_no_x_symmetry=False, force_no_y_symmetry=False, do_apodizer_throughput_maximization=True)
	write_fits(apodizer, fname + '.fits')

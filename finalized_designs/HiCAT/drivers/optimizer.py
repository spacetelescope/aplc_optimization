import os

import gurobipy as gp
import numpy as np
from hcipy import *
from scipy.ndimage.morphology import grey_erosion, grey_dilation


def calculate_pixels_to_optimize(last_optim, pupil_subsampled, edge_width_for_prior):
	"""Calculate the pixels to be used for the optimization for the adaptive algorithm.
    The exact selection of the pixels to be taken into account depends on future research.

    Parameters
    ----------
    last_optim : Field
        The previously-optimized apodizer.
    pupil_subsampled : Field
        The telescope pupil subsampled to the same resolution as `last_optim`.
    edge_pixels_for_prior : int
        The number of pixels along an edge.

    Returns
    -------
    Field
        A boolean mask to indicate which pixels to take into account.
    """
	if last_optim is None:
		return pupil_subsampled > 0

	if edge_width_for_prior == 2:
		structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
	elif edge_width_for_prior == 4:
		structure = np.array([[0, 1, 1, 0], [1, 1, 1, 1], [1, 1, 1, 1], [0, 1, 1, 0]])
	elif edge_width_for_prior == 6:
		structure = np.array([[0, 1, 1, 1, 0], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [0, 1, 1, 1, 0]])
	else:
		print('The number of edge pixels for prior is not allowed. It needs to be in [2, 4, 6].')
		exit()

	a = (grey_dilation(last_optim.shaped, structure=structure) - grey_erosion(last_optim.shaped, structure=structure)).ravel() - 2
	a = np.abs(a) > 1e-3
	b = np.logical_and(last_optim < (1 - 1e-3), last_optim > 1e-3)
	c = np.logical_or(a, b)

	return np.logical_and(c, pupil_subsampled > 0)


def optimize_aplc(pupil, focal_plane_mask, lyot_stops, dark_zone_mask, wavelengths, contrast,
				  starting_scale=1, ending_scale=1, force_no_x_symmetry=False, force_no_y_symmetry=False,
				  force_no_hermitian_symmetry=False, maximize_planet_throughput=True, num_throughput_iterations=2,
				  initial_throughput_estimate=1, edge_width_for_prior=2, solver_num_threads=0, solver_crossover=0,
				  solver_method=2, debug=False):
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
    starting_scale : int
        The number of pixels per unit cell for the initial solution. This is used for the
        adaptive algorithm. It must be a power of 2 times the `ending_scale`.
    ending_scale : int
        The number of pixels per unit cell for the final solution. If this is the same as `starting_scale`,
        the adaptive algorithm is essentially turned off.
    force_no_x_symmetry : boolean
        Force the algorithm to ignore any x mirror symmetry that might exist in the problem.
    force_no_y_symmetry : boolean
        Force the algorithm to ignore any y mirror symmetry that might exist in the problem.
    force_no_hermitian_symmetry : boolean
        Force the algorithm to ignore the (always present) Hermitian symmetry in the problem.
    maximize_planet_throughput : boolean
        If this is True, the throughput through the apodizer times Lyot stop is maximized, a metric for
        planet throughput at large angular separations. Otherwise the throughput of the apodizer is
        maximized.
    num_throughput_iterations : integer
        The number of iterations to let the throughput factor converge. Too few iterations and the contrast might be
        a fraction too high, too high and you are wasting computation time. Setting this to 1 turns iterations off.
    initial_throughput_estimate : float
        The expected relative throughput of the coronagraph compared to without apodizer or focal-plane
        mask (but including Lyot stop). A good estimate (for example from a lower resolution optimization)
        can often eliminate the need for throughput iterations altogether.
    edge_width_for_prior : [2,4,6]
        The width of the optimized regions along the edges. This is only used for the adaptive algorithm.
    solver_num_threads : integer
        The number of threads that Gurobi is allowed to use. A value of 0 indicates that Gurobi is free to
        choose the number.
    solver_crossover : integer
        The crossover strategy that Gurobi needs to use. See Gurobi documentation for more information.
    solver_method : integer
        The algorithm that Gurobi needs to use. See Gurobi documentation for more information.
    debug : boolean
        Output specific debugging information to a debug folder. This is mostly used to monitor the adaptive algorithm.

    Returns
    -------
    Field
        The optimized apodizer for the APLC. This has not been multiplied by the telescope pupil.
    """
	if debug:
		debug_dir = './debug/'
		os.makedirs(debug_dir)

		import matplotlib as mpl
		mpl.use('Agg')
		import matplotlib.pyplot as plt

	pupil_grid = pupil.grid
	focal_grid = dark_zone_mask.grid

	aplc = LyotCoronagraph(pupil_grid, focal_plane_mask)

	focal_grid_0 = CartesianGrid(UnstructuredCoords([np.array([0.0]), np.array([0.0])]), np.array([1.0]))
	prop_0 = FraunhoferPropagator(pupil_grid, focal_grid_0)

	prior = None

	# Calculate subsamplings
	num_subsamplings = int(round(np.log2(starting_scale / ending_scale))) + 1
	subsamplings = ending_scale * 2 ** np.arange(num_subsamplings)[::-1]

	# Determine pupil symmetries
	x_symm_pupil = np.allclose(pupil.shaped[:, ::-1], pupil.shaped)
	y_symm_pupil = np.allclose(pupil.shaped[::-1, :], pupil.shaped)

	print('Telescope pupil:')
	print('   Mirror symmetry in x: %s' % ('yes' if x_symm_pupil else 'no'))
	print('   Mirror symmetry in y: %s' % ('yes' if y_symm_pupil else 'no'))
	print('')

	# Determine lyot stop symmetries
	x_symm_lyot_stops = [np.allclose(lyot_stop.shaped[:, ::-1], lyot_stop.shaped) for lyot_stop in lyot_stops]
	y_symm_lyot_stops = [np.allclose(lyot_stop.shaped[::-1, :], lyot_stop.shaped) for lyot_stop in lyot_stops]

	lyot_stop_duplication = np.zeros(len(lyot_stops), dtype='bool')
	lyot_stop_duplication_reason = [[] for i in range(len(lyot_stops))]

	x_symm = x_symm_pupil
	y_symm = y_symm_pupil

	# X mirror symmetry for Lyot stops
	for i, a in enumerate(lyot_stops):
		print('Lyot stop #%d:' % i)

		if lyot_stop_duplication[i]:
			print('   Will be ignored due to symmetries with Lyot stops #' + str(lyot_stop_duplication_reason[i]))
			continue

		print('   Mirror symmetry in x: %s' % ('yes' if x_symm_lyot_stops[i] else 'no'))

		if x_symm_pupil and not x_symm_lyot_stops[i]:
			print('   Searching for mirror symmetric Lyot stops in x...')
			for j, b in enumerate(lyot_stops):
				if j <= i:
					continue

				if np.allclose(a.shaped[:, ::-1], b.shaped):
					print('      Found Lyot stop #%d to fit.' % j)
					lyot_stop_duplication[j] = True
					lyot_stop_duplication_reason[j].append(i)
					x_symm = True
					break
			else:
				print(
					'      No Lyot stop found with this symmetry. This breaks the mirror symmetry in x of the optimization.')
				x_symm = False

	# Y mirror symmetry for Lyot stops
	for i, a in enumerate(lyot_stops):
		print('Lyot stop #%d:' % i)

		if lyot_stop_duplication[i]:
			print('   Will be ignored due to symmetries with Lyot stops #' + str(lyot_stop_duplication_reason[i]))
			continue

		print('   Mirror symmetry in y: %s' % ('yes' if y_symm_lyot_stops[i] else 'no'))

		if y_symm_pupil and not y_symm_lyot_stops[i]:
			print('   Searching for mirror symmetric Lyot stops in y...')
			for j, b in enumerate(lyot_stops):
				if j <= i:
					continue

				if np.allclose(a.shaped[::-1, :], b.shaped):
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

	num_constraints_per_focal_point = np.array(num_constraints_per_focal_point)

	# Calculate dark zone masks, and propagators for each Lyot stops.
	# These can be different due to the symmetry properties of each Lyot stop,
	# and thus need to be calculated separately.
	dark_zone_masks = []
	dark_zone_masks_full = []
	propagators = []
	num_focal_points = []
	for i, lyot_stop in enumerate(lyot_stops):
		if lyot_stop_duplication[i]:
			# Lyot stop was removed due to symmetries
			dark_zone_masks.append(None)
			dark_zone_masks_full.append(np.zeros(focal_grid.size, dtype='bool'))
			propagators.append(None)
			num_focal_points.append(0)
			continue

		# Separated coords for sub focal grid
		x, y = focal_grid.separated_coords

		# Hermitian symmetry (all planes are real)
		if not force_no_hermitian_symmetry:
			m = focal_grid.x > 0
			x = x[x > 0]

		# Only optimize quarter of roi is mirror symmetric in one or two axes
		if (x_symm and x_symm_lyot_stops[i]) or (y_symm and y_symm_lyot_stops[i]):
			m *= focal_grid.y > 0
			y = y[y > 0]

		# Make grid with subset of focal grid
		focal_grid_sub = CartesianGrid(SeparatedCoords((x, y)))
		# focal_grid_sub = focal_grid

		# Make propagator for this Lyot stop
		propagators.append(FraunhoferPropagator(pupil_grid, focal_grid_sub))

		# Recalculate dark zone mask for the sub focal grid
		dark_zone_masks.append(Field(dark_zone_mask[m > 0], focal_grid_sub).astype('bool'))
		dark_zone_masks_full.append((dark_zone_mask * m).astype('bool'))
		# dark_zone_masks.append((dark_zone_mask * m).astype('bool'))

		# Calculating number of focal points
		num_focal_points.append(int(np.sum(dark_zone_masks[-1])))

	num_focal_points = np.array(num_focal_points)

	# Calculate number of constraints per wavelength
	m = int(np.sum(num_constraints_per_focal_point * num_focal_points))
	cum_mms = np.concatenate(([0], np.cumsum(num_constraints_per_focal_point * num_focal_points))).astype('int')

	throughput_estimate = initial_throughput_estimate

	# Iterate from lowest to highest resolution (=highest to lowest subsampling)
	for subsampling_index, subsampling in enumerate(subsamplings):
		# Calculated subsampled pupil and lyot stops
		pupil_subsampled = subsample_field(pupil, subsampling)
		lyot_stops_subsampled = [subsample_field(lyot_stop, subsampling) for lyot_stop in lyot_stops]
		pupil_grid_subsampled = pupil_subsampled.grid

		# Calculate which pixels belong to which superpixel
		inds = np.arange(pupil_grid.size).reshape(
			(pupil_grid.shape[1] // subsampling, subsampling, pupil_grid.shape[0] // subsampling, subsampling))
		inds = np.swapaxes(inds, 1, 2).reshape((pupil_grid_subsampled.shape[0], pupil_grid_subsampled.shape[1],
												-1))  # .reshape((pupil_grid.size//(subsampling**2), -1))

		# Apply x,y-mirror-symmetries
		symmetry_mask = Field(np.ones(pupil_grid_subsampled.size), pupil_grid_subsampled)
		if x_symm:
			inds = np.concatenate((inds, inds[:, ::-1, :]), axis=2)
			symmetry_mask *= pupil_grid_subsampled.x < 0
		if y_symm:
			inds = np.concatenate((inds, inds[::-1, :, :]), axis=2)
			symmetry_mask *= pupil_grid_subsampled.y < 0

		symmetry_mask = symmetry_mask.astype('bool')
		inds = inds.reshape((pupil_grid_subsampled.size, -1))

		# Upscale last optim to current resolution
		blind = prior is None
		if blind:
			# No prior information; assume totally dark apodizer
			last_optim = Field(np.zeros(pupil_grid_subsampled.size), pupil_grid_subsampled)
			prior = Field(np.zeros(pupil_grid.size), pupil_grid)
		else:
			# Upscale prior information by factor 2
			last_optim = subsample_field(prior, subsampling)

		# Write prior to file
		if debug:
			write_fits(prior * pupil, debug_dir + 'prior%d.fits' % subsampling)
			write_fits(last_optim * pupil_subsampled, debug_dir + 'last_optim%d.fits' % subsampling)

		# Get pixels to optimize
		optimize_mask = np.logical_and(calculate_pixels_to_optimize(last_optim, pupil_subsampled, edge_width_for_prior),
									   symmetry_mask)
		if blind:
			optimize_mask[:] = np.logical_and(pupil_subsampled > 0, symmetry_mask)
		n = int(np.sum(optimize_mask * symmetry_mask))

		# Write optimize mask to file
		if debug:
			write_fits(optimize_mask.astype('int'), debug_dir + 'optimize_mask%d.fits' % subsampling)

		print('Starting optimization at scale %d with %d variables and %d constraints.' % (
			subsampling, n, m * len(wavelengths)))
		print('Creating model...')

		# Create Gurobi model
		model = gp.Model('lp')
		model.Params.Threads = solver_num_threads
		model.Params.Crossover = solver_crossover
		model.Params.Method = solver_method
		x_vars = model.addVars(n, lb=0, ub=1)

		print('Calculating and adding constraints...')

		# Create problem matrix for one wavelength but for all Lyot stops
		M = np.empty((m, n))
		base_electric_fields = []
		contrast_requirements = []
		constraints = []

		# Add constraints for each wavelength
		for wl_i, wavelength in enumerate(wavelengths):
			j = 0
			x0 = Field(np.zeros(pupil_grid.size, dtype='complex'), pupil_grid)
			x = Field(np.zeros(pupil_grid.size, dtype='complex'), pupil_grid)

			# Calculate norm electric field for each Lyot stop
			norms = []
			for lyot_stop in lyot_stops:
				norms.append(prop_0(Wavefront(pupil * lyot_stop, wavelength)).electric_field[0])

			for ind, amp, to_optimize, masked_by_symmetry in zip(inds, last_optim, optimize_mask, symmetry_mask):
				if not to_optimize:
					# Do not optimize this pixel
					# Add to accumulator pupil-plane wavefront
					if masked_by_symmetry:
						x0[ind] += pupil[ind] * amp
				else:
					x[:] = 0
					x[ind] = pupil[ind]

					# Calculate field before the Lyot stop
					lyot = aplc(Wavefront(x, wavelength))

					k = 0
					for i, lyot_stop in enumerate(lyot_stops):
						if num_constraints_per_focal_point[i] == 0:
							continue

						# Apply the Lyot stop and get focal-plane electric field
						lyot_copy = lyot.copy()
						lyot_copy.electric_field *= lyot_stop
						E = propagators[i](lyot_copy).electric_field[dark_zone_masks[i]]
						E /= norms[i]

						# Add only imaginary or both imaginary and real constraints depending on symmetry
						if num_constraints_per_focal_point[i] == 1:
							M[cum_mms[i]:cum_mms[i + 1], j] = E.real
						if num_constraints_per_focal_point[i] == 2:
							M[cum_mms[i]:cum_mms[i] + num_focal_points[i], j] = E.real + E.imag
							M[cum_mms[i] + num_focal_points[i]:cum_mms[i + 1], j] = E.real - E.imag

						k += num_constraints_per_focal_point[i]
					j += 1

					if j % 1000 == 0:
						print('Wavelength %d/%d; Variable %d/%d' % (wl_i + 1, len(wavelengths), j, n))

			# Display x0
			if debug:
				imshow_field(x0.real)
				plt.colorbar()
				plt.savefig(debug_dir + 'base_aperture%d.pdf' % subsampling)
				plt.clf()

			# Calculate base electric field
			base_electric_field = []
			for i, lyot_stop in enumerate(lyot_stops):
				if num_constraints_per_focal_point[i] == 0:
					continue

				wf = aplc(Wavefront(x0, wavelength))
				wf.electric_field *= lyot_stop
				E = propagators[i](wf).electric_field[dark_zone_masks[i]]
				E /= norms[i]

				if num_constraints_per_focal_point[i] == 1:
					base_electric_field.append(E.real)
				if num_constraints_per_focal_point[i] == 2:
					base_electric_field.append(E.real + E.imag)
					base_electric_field.append(E.real - E.imag)
			base_electric_field = np.concatenate(base_electric_field)
			base_electric_fields.append(base_electric_field)

			# Calculate contrast requirement
			contrast_requirement = []
			for p, q, r in zip(num_focal_points, num_constraints_per_focal_point, dark_zone_masks_full):
				temp = np.repeat((np.ones(focal_grid.size) * np.sqrt(contrast))[r], q)
				contrast_requirement.append(temp)
			contrast_requirement = np.concatenate(contrast_requirement)
			contrast_requirements.append(contrast_requirement.copy())

			contrast_requirement *= np.sqrt(throughput_estimate)

			# Add constraints
			for ee, e0, c0 in zip(M, base_electric_field, contrast_requirement):
				e = gp.LinExpr(ee, x_vars.values())
				c1 = model.addConstr(e <= (c0 - e0))
				c2 = model.addConstr(e >= (-c0 - e0))

				constraints.extend([c1, c2])

		base_electric_fields = np.concatenate(base_electric_fields)
		contrast_requirements = np.concatenate(contrast_requirements)

		del M

		# Use central Lyot stop for throughput metric (assume that this is the unshifted Lyot stop)
		if maximize_planet_throughput:
			M_max = (pupil_subsampled * lyot_stops_subsampled[0])[optimize_mask]
		else:
			M_max = pupil_subsampled[optimize_mask]
		obj = gp.LinExpr(M_max, x_vars.values())
		model.setObjective(obj, gp.GRB.MAXIMIZE)

		# Optimize model
		print('Start optimization...')
		model.optimize()
		print('Optimization finished!')

		for throughput_iter in range(num_throughput_iterations):
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

			throughput_estimate = ((pupil * prior * lyot_stop).sum() / (pupil * lyot_stop).sum()) ** 2

			# Only do throughput iterations on the first iteration for the adaptive algorithm.
			# Doing it on later iterations might make the model infeasible.
			# TODO: This might be overkill, but this is to be investigated.
			# Also, only modify and reoptimize if this is not the last iteration (mimicing do-while loop).
			if not (throughput_iter == num_throughput_iterations - 1) and subsampling_index == 0:
				# Update throughput estimate and rerun
				rhs = np.empty(base_electric_fields.size * 2)
				rhs[::2] = contrast_requirements * np.sqrt(throughput_estimate) - base_electric_fields
				rhs[1::2] = -contrast_requirements * np.sqrt(throughput_estimate) - base_electric_fields
				for constr, r in zip(constraints, rhs):
					constr.RHS = r

				# Rerun optimizer
				print('Restarting optimization for throughput iteration %d...' % (throughput_iter + 1))
				model.optimize()
				print('Optimization finished!')

		# Write result to display
		if debug:
			write_fits(prior * (pupil > 0), debug_dir + 'result%d.fits' % subsampling)

	return prior * (pupil > 0)

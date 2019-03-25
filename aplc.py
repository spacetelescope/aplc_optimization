import matplotlib as mpl
mpl.use('Agg')
from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from scipy.ndimage.morphology import grey_erosion, grey_dilation

def calculate_pixels_to_optimize(last_optim, pupil_subsampled):
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

def optimize_aplc(pupil, focal_plane_mask, lyot_stops, dark_zone_mask, wavelengths, contrast, num_scalings=1):
	pupil_grid = pupil.grid
	focal_grid = dark_zone_mask.grid

	aplc = LyotCoronagraph(pupil_grid, focal_plane_mask)
	prop = FraunhoferPropagator(pupil_grid, focal_grid)

	last_optim = None
	subsamplings = 2**np.arange(num_scalings)[::-1]

	# Determine pupil symmetries
	x_symm_pupil = np.allclose(pupil.shaped[:,::-1], pupil.shaped)
	y_symm_pupil = np.allclose(pupil.shaped[::-1,:], pupil.shaped)

	if x_symm_pupil:
		if y_symm_pupil:
			print('Telescope pupil is both x- and y-symmetric.')
		else:
			print('Telescope pupil is x-symmetric.')
	elif y_symm_pupil:
		print('Telescope pupil is y-symmetric.')
	else:
		print('Telescope pupil contains no mirror symmetries.')
	
	# Determine lyot stop symmetries
	x_symm_lyot_stops = [np.allclose(lyot_stop.shaped[:,::-1], lyot_stop.shaped) for lyot_stop in lyot_stops]
	y_symm_lyot_stops = [np.allclose(lyot_stop.shaped[::-1,:], lyot_stop.shaped) for lyot_stop in lyot_stops]

	# Determine problem symmetries of the entire optimization problem
	if x_symm_pupil:
		if np.all(x_symm_lyot_stops):
			x_symm = True
		else:
			x_symm = True
			for a in lyot_stops:
				if not np.any([np.allclose(a[:,::-1], b) for b in lyot_stops]):
					x_symm = False
	else:
		x_symm = False
	
	if y_symm_pupil:
		if np.all(y_symm_lyot_stops):
			y_symm = True
		else:
			y_symm = True
			for a in lyot_stops:
				if not np.any([np.allclose(a[::-1,:], b) for b in lyot_stops]):
					y_symm = False
	else:
		y_symm = False

	# Remove Lyot stops from the list if they are covered by symmetry from others
	lyot_stop_duplication = np.zeros(len(lyot_stops), dtype='bool')
	for i, a in enumerate(lyot_stops):
		if x_symm:
			if not x_symm_lyot_stops[i]:
				for j, b in enumerate(lyot_stops):
					if j <= i:
						continue
					if np.allclose(a[:,::-1], b):
						lyot_stop_duplication[i] = True
		if y_symm:
			if not y_symm_lyot_stops[i]:
				for j, b in enumerate(lyot_stops):
					if j <= i:
						continue
					if np.allclose(a[::-1,:], b):
						lyot_stop_duplication[i] = True

	# Find number of constraints per focal-plane point
	num_constraints_per_focal_point = []
	for i, lyot_stop in enumerate(lyot_stops):
		if lyot_stop_duplication[i]:
			# Lyot stop already taken into account due to symmetries
			num_constraints_per_focal_point.append(0)
		elif (x_symm_pupil and x_symm_lyot_stops[i]) and (y_symm_pupil and y_symm_lyot_stops[i]):
			# Only constrain real part
			num_constraints_per_focal_point.append(1)
		else:
			# Constrain both real and imag part
			num_constraints_per_focal_point.append(2)
	mm = int(np.sum(dark_zone_mask))
	
	# We are optimizing amplitude (=real part) only
	dark_zone_mask *= focal_grid.x > 0

	if x_symm or y_symm:
		dark_zone_mask *= focal_grid.y > 0
	
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
		inds = np.swapaxes(inds, 1, 2)#.reshape((pupil_grid.size//(subsampling**2), -1))

		# Apply x,y-mirror-symmetries
		mask = Field(np.ones(pupil_grid_subsampled.size), pupil_grid_subsampled)
		if x_symm:
			inds = np.concatenate((inds, inds[:,::-1,:]))
			mask *= pupil_grid_subsampled.x > 0
		if y_symm:
			inds = np.concatenate((inds, inds[::-1,:,:]))
			mask *= pupil_grid_subsampled.y > 0
		
		mask = mask.astype('bool')
		inds = inds.reshape((pupil_grid_subsampled.size, -1))
		inds = [inds[i,:] for i in range(pupil_grid_subsampled.size) if mask[i]]

		# Upscale last optim to current resolution
		blind = last_optim is None
		if blind:
			# No prior information; assume totally dark apodizer
			last_optim = Field(np.zeros(pupil_grid_subsampled.size), pupil_grid_subsampled)
		else:
			# Upscale prior information by factor 2
			last_optim = np.repeat(np.repeat(last_optim.shaped, 2, 1), 2, 0).ravel()
			last_optim = Field(last_optim, pupil_grid_subsampled)
		
		# Get pixels to optimize
		optimize_mask = calculate_pixels_to_optimize(last_optim, pupil_subsampled)[mask]
		if blind:
			optimize_mask[:] = pupil_subsampled > 0
		n = int(np.sum(optimize_mask))

		# Create Gurobi model
		model = gp.Model('lp')
		model.Params.Threads = 0
		model.Params.Crossover = 0
		model.Params.Method = 2
		x_vars = model.addVars(n, lb=0, ub=1)

		# Create problem matrix for one wavelength but for all Lyot stops
		M = np.empty((m, n))

		# Add constraints for each wavelength
		for wavelength in wavelengths:
			j = 0
			x0 = Field(np.zeros(pupil_grid.size), pupil_grid)
			x = Field(np.zeros(pupil_grid.size), pupil_grid)
			
			for ind, amp, to_optimize in zip(inds, last_optim, optimize_mask):
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

						# Add only imaginary or both imaginary and real constraints depending on symmetry
						M[k*mm:(k+1)*mm,j] = E.imag
						if num_constraints_per_focal_point[i] == 2:
							M[(k+1)*mm:(k+2)*mm,j] = E.real
						
						k += num_constraints_per_focal_point[i]
					j += 1

			# Calculate base electric field
			base_electric_field = []
			for i, lyot_stop in enumerate(lyot_stops):
				if num_constraints_per_focal_point[i] == 0:
					continue
				
				wf = aplc(Wavefront(x0, wavelength))
				wf.electric_field *= lyot_stop
				img = prop(wf)
				
				base_electric_field.append(img.electric_field.imag)
				if num_constraints_per_focal_point[i] == 2:
					base_electric_field.append(img.electric_field.real)
			base_electric_field = np.concatenate(base_electric_field)

			# Calculate contrast requirement
			contrast_requirement = np.tile([np.ones(m) * np.sqrt(contrast[dark_zone_mask])], int(np.sum(num_constraints_per_focal_point)))

			# Add constraints
			for ee, e0, c0 in zip(M, base_electric_field, contrast_requirement):
				e = gp.LinExpr(ee, x_vars.values())
				model.addConstr(e <= (-e0 + c0))
				model.addConstr(e >= (-e0 - c0))
		
		del M

		# Use central Lyot stop for throughput metric (assume that this is the unshifted Lyot stop)
		M_max = (pupil_subsampled * lyot_stops_subsampled[0])
		obj = gp.LinExpr(M_max, x_vars.values())
		model.setObjective(obj, gp.GRB.MAXIMIZE)

		# Optimize model
		model.optimize()

		# Extract solution from Gurobi
		solution = np.array([x_vars[i].x for i in range(n)])

		# Integrate solution into upsampled old solution
		sol = last_optim
		sol[optimize_mask] = solution
		sol = Field(sol, make_subsampled_grid(pupil_grid, subsampling))
		
		last_optim = sol
	
	return last_optim
	
if __name__ == '__main__':
	contrast = 1e-8
	num_pix = 486
	q_sci = 2 # px / (lambda_0/D)
	iwa = 3.75 # lambda_0/D
	owa = 15 # lambda_0/D
	n_foc = 50 # px diameter
	foc_inner = 8.543 # lambda_0/D diameter
	spectral_bandwidth = 0.1 # fractional
	num_wavelengths = 3
	
	testing = True

	if testing:
		num_pix = 256
		owa = 8
		num_wavelengths = 2

	pupil_grid = make_pupil_grid(num_pix)

	if testing:
		pupil = evaluate_supersampled(make_hicat_aperture(True), pupil_grid, 4)
	else:
		pupil = read_fits('masks/SYM-HiCAT-Aper_F-N0486_Hex3-Ctr0972-Obs0195-SpX0017-Gap0004.fits')
		pupil = Field(pupil.ravel(), pupil_grid)
	
	if testing:
		lyot_stop = evaluate_supersampled(make_hicat_lyot_stop(True), pupil_grid, 4)
	else:
		lyot_stop = read_fits('masks/HiCAT-Lyot_F-N0486_LS-Ann-bw-ID0345-OD0807-SpX0036.fits')
		lyot_stop = Field(lyot_stop.ravel(), pupil_grid)

	lyot_stops = [lyot_stop]
	
	n_sci = int((np.ceil(owa) + 1) * q_sci) * 2
	x_sci = (np.arange(n_sci) + 0.5 - n_sci / 2) / q_sci
	focal_grid = CartesianGrid(SeparatedCoords((x_sci, x_sci)))

	dark_zone_mask = circular_aperture(owa * 2)(focal_grid) - circular_aperture(iwa * 2)(focal_grid)

	q_foc = n_foc / foc_inner
	x_foc = (np.arange(n_foc) + 0.5 - n_foc / 2) / q_foc
	focal_mask_grid = CartesianGrid(SeparatedCoords((x_foc, x_foc)))

	focal_plane_mask = circular_aperture(foc_inner)(focal_mask_grid)

	if num_wavelengths == 1:
		wavelengths = [1]
	else:
		wavelengths = np.linspace(-spectral_bandwidth / 2, spectral_bandwidth / 2, num_wavelengths) + 1

	res = optimize_aplc(pupil, focal_plane_mask, lyot_stops, dark_zone_mask, wavelengths, contrast, num_scalings=1)

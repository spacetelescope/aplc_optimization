from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from scipy.ndimage.morphology import grey_erosion, grey_dilation

pupil_grid = make_pupil_grid(512)
focal_grid = CartesianGrid(SeparatedCoords((np.arange(3.5, 14, 0.5), np.arange(-14, 14.1, 0.5))))

prop = FraunhoferPropagator(pupil_grid, focal_grid).get_monochromatic_propagator(1)

aperture = make_obstructed_circular_aperture(1, 0.2, 3, 0.01)
aperture = evaluate_supersampled(aperture, pupil_grid, 4)

#imshow_field(aperture)
#plt.show()

def f(i, prop):
	x = Field(np.zeros(prop.input_grid.size), prop.input_grid)
	x[i] = 1
	return prop(Wavefront(x)).electric_field

def get_propagation_matrix(prop, aperture, subsampling):
	pupil_grid = prop.input_grid
	focal_grid = prop.output_grid

	inds = np.arange(pupil_grid.size).reshape((pupil_grid.shape[1]//subsampling, subsampling, pupil_grid.shape[0]//subsampling, subsampling))
	inds = np.swapaxes(inds, 1, 2).reshape((pupil_grid.size//(subsampling**2), -1))

	mat = []
	for ind in inds:
		#print(ind)
		x = Field(np.zeros(prop.input_grid.size), prop.input_grid)
		x[ind] = aperture[ind]
		y = prop(Wavefront(x)).electric_field
		mat.append(y)
	
	return np.array(mat)

def optimize_at_scale(prop, aperture, last_optim, subsampling, contrast):
	pupil_grid = prop.input_grid
	focal_grid = prop.output_grid

	inds = np.arange(pupil_grid.size).reshape((pupil_grid.shape[1]//subsampling, subsampling, pupil_grid.shape[0]//subsampling, subsampling))
	inds = np.swapaxes(inds, 1, 2).reshape((pupil_grid.size//(subsampling**2), -1))

	last_optim = np.repeat(np.repeat(last_optim.shaped, 2, 1), 2, 0).ravel()
	aperture_subsampled = subsample_field(aperture, subsampling)
	last_optim = Field(last_optim, aperture_subsampled.grid)

	mat = []
	base_electric_field = np.zeros(focal_grid.size, dtype='complex')
	optimize_mask = []

	pixels_to_optimize = (grey_dilation(last_optim.shaped, 4) - grey_erosion(last_optim.shaped, 4)).ravel()
	pixels_to_optimize = Field(pixels_to_optimize, last_optim.grid)
	pixels_to_optimize = np.abs(pixels_to_optimize) > 1e-3

	pixels_to_optimize2 = np.logical_and(last_optim < (1 - 1e-3), last_optim > 1e-3)

	pixels_to_optimize = np.logical_or(pixels_to_optimize, pixels_to_optimize2)

	imshow_field(pixels_to_optimize, cmap='Reds', mask=aperture_subsampled, vmin=0, vmax=1.5)
	plt.savefig('pixels_to_optimize%d.pdf' % subsampling)
	plt.clf()

	for ind, amp, to_optimize in zip(inds, last_optim, pixels_to_optimize):
		if np.sum(aperture[ind]) < 1e-3:
			optimize_mask.append(False)
		else:
			x = Field(np.zeros(pupil_grid.size), pupil_grid)
			x[ind] = aperture[ind]
			y = prop(Wavefront(x)).electric_field
			
			if not to_optimize:
				base_electric_field += amp * y
				optimize_mask.append(False)
			else:
				mat.append(y)
				optimize_mask.append(True)

	optimize_mask = np.array(optimize_mask)
	base_electric_field = np.array(base_electric_field)

	M = np.array(mat).T
	M = np.vstack([M.real, M.imag])

	M_max = aperture_subsampled[optimize_mask]

	E_0 = np.concatenate((base_electric_field.real, base_electric_field.imag))

	m, n = M.shape
	print(m, n)

	contrast_requirement = np.ones(m) * np.sqrt(contrast)

	model = gp.Model('lp')
	model.Params.Threads = 4

	x = model.addVars(n, lb=0, ub=1)

	obj = gp.quicksum((x[i] * M_max[i] for i in range(n)))
	model.setObjective(obj, gp.GRB.MAXIMIZE)

	for i, ee in enumerate(M):
		e = gp.quicksum((x[i] * ee[i] for i in range(n)))
		model.addConstr(e <= -E_0[i] + contrast_requirement[i])
		model.addConstr(e >= -E_0[i] - contrast_requirement[i])
	
	model.optimize()

	solution = np.array([x[i].x for i in range(n)])

	sol = last_optim
	sol[optimize_mask] = solution
	sol = Field(sol, make_subsampled_grid(pupil_grid, subsampling))
	imshow_field(sol, cmap='gray', mask=aperture_subsampled)
	plt.savefig('shaped_pupil_%d.pdf' % subsampling)
	plt.clf()

	return sol

	print(M.sum())
	print(M_max.shape)
	return
	m, n = M.shape
	print(m, n)

g = make_subsampled_grid(pupil_grid, 32)
last_optim = Field(np.ones(g.size) * 0.5, g)
contrast = 1e-5

for subsampling in [16,8,4,2,1]:
	last_optim = optimize_at_scale(prop, aperture, last_optim, subsampling, contrast)
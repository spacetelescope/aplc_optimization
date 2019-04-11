from hcipy import *
import numpy as np
import asdf
from astropy.io import fits
import os
import pprint
from por_aplc_optimizer import optimize_aplc

print('Optimizer called with the following parameters:')
pprint.pprint(parameters, width=1)
print('')
print('File organization:')
pprint.pprint(file_organization, width=1)
print('')

# Getting values from the template
# Driver is automatically written with "parameters" and "file_organization"
pup_fname = parameters['pupil']['filename']
fpm_radius = parameters['focal_plane_mask']['radius']
fpm_num_pix = parameters['focal_plane_mask']['num_pix']
fpm_grayscale = parameters['focal_plane_mask']['grayscale']
ls_fname = parameters['lyot_stop']['filename']
ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
ls_num_stops = parameters['lyot_stop']['num_lyot_stops']
img_contrast = 10**(-parameters['image']['contrast'])
img_iwa = parameters['image']['iwa']
img_owa = parameters['image']['owa']
img_num_wavelengths = parameters['image']['num_wavelengths']
img_bandwidth = parameters['image']['bandwidth']
img_resolution = parameters['image']['resolution']
method_force_no_x_mirror_symmetry = parameters['method']['force_no_x_mirror_symmetry']
method_force_no_y_mirror_symmetry = parameters['method']['force_no_y_mirror_symmetry']
method_force_no_hermitian_symmetry = parameters['method']['force_no_hermitian_symmetry']
method_starting_scale = parameters['method']['starting_scale']
method_ending_scale = parameters['method']['ending_scale']
method_edge_width_for_prior = parameters['method']['edge_width_for_prior']
method_num_throughput_iterations = parameters['method']['num_throughput_iterations']
method_initial_throughput_estimate = parameters['method']['initial_throughput_estimate']
method_maximize_planet_throughput = parameters['method']['maximize_planet_throughput']
solver_num_threads = parameters['solver']['num_threads']
solver_crossover = parameters['solver']['crossover']
solver_method = parameters['solver']['method']

if not os.path.isabs(pup_fname):
	pup_fname = os.path.join(file_organization['input_files_dir'], pup_fname)
if not os.path.isabs(ls_fname):
	ls_fname = os.path.join(file_organization['input_files_dir'], ls_fname)

pupil = read_fits(pup_fname)

num_pix = pupil.shape[0]
pupil_grid = make_uniform_grid(num_pix, [1, 1])
pupil = Field(pupil.ravel(), pupil_grid)

try:
	lyot_stops = [Field(read_fits(ls_fname.format(i)).ravel(), pupil_grid) for i in range(ls_num_stops)]
except:
	lyot_stop = Field(read_fits(ls_fname).ravel(), pupil_grid)

	# Building Lyot stops according to alignment tolerance
	if ls_num_stops in [1, 5, 9]:
		lyot_stops = [lyot_stop]
	else:
		lyot_stops = []

	if ls_num_stops in [4, 5, 8, 9]:
		lyot_stop_pos_x = np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1).ravel()
		lyot_stop_neg_x = np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1).ravel()
		lyot_stop_pos_y = np.roll(lyot_stop.shaped, ls_alignment_tolerance, 0).ravel()
		lyot_stop_neg_y = np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 0).ravel()

		lyot_stops.extend([lyot_stop_pos_x, lyot_stop_neg_x, lyot_stop_pos_y, lyot_stop_neg_y])

	if ls_num_stops in [8, 9]:
		lyot_stop_pos_x_pos_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1), ls_alignment_tolerance, 0).ravel()
		lyot_stop_pos_x_neg_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1), -ls_alignment_tolerance, 0).ravel()
		lyot_stop_neg_x_pos_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1), ls_alignment_tolerance, 0).ravel()
		lyot_stop_neg_x_neg_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1), -ls_alignment_tolerance, 0).ravel()

		lyot_stops.extend([lyot_stop_pos_x_pos_y, lyot_stop_pos_x_neg_y, lyot_stop_neg_x_pos_y, lyot_stop_neg_x_neg_y])

# Build science focal grid
n_sci = int((np.ceil(img_owa) + 1) * img_resolution) * 2
x_sci = (np.arange(n_sci) + 0.5 - n_sci / 2) / img_resolution
focal_grid = CartesianGrid(SeparatedCoords((x_sci, x_sci)))

dark_zone_mask = circular_aperture(img_owa * 2)(focal_grid) - circular_aperture(img_iwa * 2)(focal_grid)

# Build focal plane mask
q_foc = fpm_num_pix / (fpm_radius * 2)
x_foc = (np.arange(fpm_num_pix) + 0.5 - fpm_num_pix / 2) / q_foc
focal_mask_grid = CartesianGrid(RegularCoords(1.0 / q_foc, [fpm_num_pix, fpm_num_pix], x_foc.min()))

if fpm_grayscale:
	focal_plane_mask = 1 - evaluate_supersampled(circular_aperture(2 * fpm_radius), focal_mask_grid, 8)
else:
	focal_plane_mask = 1 - circular_aperture(2 * fpm_radius)(focal_mask_grid)

if img_num_wavelengths == 1:
	wavelengths = [1]
else:
	wavelengths = np.linspace(-img_bandwidth / 2, img_bandwidth / 2, img_num_wavelengths) + 1

# Optimize
apodizer = optimize_aplc(
	pupil=pupil, 
	focal_plane_mask=focal_plane_mask, 
	lyot_stops=lyot_stops, 
	dark_zone_mask=dark_zone_mask, 
	wavelengths=wavelengths, 
	contrast=img_contrast, 
	starting_scale=method_starting_scale, 
	ending_scale=method_ending_scale, 
	force_no_x_symmetry=method_force_no_x_mirror_symmetry,
	force_no_y_symmetry=method_force_no_y_mirror_symmetry,
	force_no_hermitian_symmetry=method_force_no_hermitian_symmetry,
	maximize_planet_throughput=method_maximize_planet_throughput,
	num_throughput_iterations=method_num_throughput_iterations,
	initial_throughput_estimate=method_initial_throughput_estimate,
	edge_width_for_prior=method_edge_width_for_prior,
	solver_num_threads=solver_num_threads,
	solver_crossover=solver_crossover,
	solver_method=solver_method)

# TODO: Field stop

# Write to file
hdu_list = fits.HDUList()
hdu_list.append(fits.ImageHDU((apodizer * (pupil > 0)).shaped, name='APOD'))

tree = {'parameters': parameters, 'file_organization': file_organization, 'apodizer': hdu_list['APOD'].data}

ff = asdf.fits_embed.AsdfInFits(hdu_list, tree)
ff.write_to(solution_fname, overwrite=True)

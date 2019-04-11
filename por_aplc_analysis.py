from hcipy import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import asdf
import os

def create_coronagraph(solution_filename):
	solution = asdf.open(solution_filename)
	parameters = solution.tree['parameters']
	file_organization = solution.tree['file_organization']
	apodizer = solution.tree['apodizer']

	pup_fname = parameters['pupil']['filename']
	fpm_radius = parameters['focal_plane_mask']['radius']
	fpm_num_pix = parameters['focal_plane_mask']['num_pix']
	fpm_grayscale = parameters['focal_plane_mask']['grayscale']
	ls_fname = parameters['lyot_stop']['filename']
	ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
	ls_num_stops = parameters['lyot_stop']['num_lyot_stops']

	if not os.path.isabs(pup_fname):
		pup_fname = os.path.join(file_organization['input_files_dir'], pup_fname)
	if not os.path.isabs(ls_fname):
		ls_fname = os.path.join(file_organization['input_files_dir'], ls_fname)

	pupil = read_fits(pup_fname)

	num_pix = pupil.shape[0]
	pupil_grid = make_uniform_grid(num_pix, [1, 1])
	pupil = Field(pupil.ravel(), pupil_grid)
	apodizer = Field(apodizer.ravel(), pupil_grid)

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
		
	# Build focal plane mask
	q_foc = fpm_num_pix / (fpm_radius * 2)
	x_foc = (np.arange(fpm_num_pix) + 0.5 - fpm_num_pix / 2) / q_foc
	focal_mask_grid = CartesianGrid(RegularCoords(1.0 / q_foc, [fpm_num_pix, fpm_num_pix], x_foc.min()))

	if fpm_grayscale:
		focal_plane_mask = 1 - evaluate_supersampled(circular_aperture(2 * fpm_radius), focal_mask_grid, 8)
	else:
		focal_plane_mask = 1 - circular_aperture(2 * fpm_radius)(focal_mask_grid)
	
	return pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization

def analyze_contrast_monochromatic(solution_filename, pdf=None):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	owa = parameters['image']['owa']
	contrast = parameters['image']['contrast']

	lyot_stop = lyot_stops[0]

	coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
	focal_grid = make_focal_grid(pupil.grid, 8, owa * 1.2)
	prop = FraunhoferPropagator(pupil.grid, focal_grid)

	wf = Wavefront(apodizer)
	img = prop(coro(wf)).intensity
	img_ref = prop(Wavefront(apodizer * lyot_stop)).intensity
	
	plt.clf()
	imshow_field(np.log10(img / img_ref.max()), vmin=-contrast-1, vmax=-contrast+4, cmap='inferno')
	plt.colorbar()
	if pdf is not None:
		pdf.savefig()
	else:
		plt.show()
	
	return {'normalized_irradiance_image': img / img_ref.max()}

def analyze_max_throughput(solution_filename, pdf=None):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	lyot_stop = lyot_stops[0]

	maximum_integrated_throughput = ((pupil * lyot_stop * apodizer).sum() / (pupil * lyot_stop).sum())**2

	return {'maximum_integrated_throughput': maximum_integrated_throughput}

def analyze_offaxis_throughput(solution_filename, pdf=None):
	pass

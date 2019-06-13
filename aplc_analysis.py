import matplotlib as mpl
mpl.use('Agg')

from hcipy import *
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import asdf
import os

#stand alone analysis script

def create_coronagraph(solution_filename):
		
	parameters = {'pupil': {'filename': 'hicat_test/inputs/hicat_apodizer_mask_1944_bw.fits', 'N': 1944}, 
	'focal_plane_mask': {'radius': 8.543/2, 'num_pix': 80, 'grayscale': True, 'field_stop_radius': -1.0}, 
	'lyot_stop': {'filename': 'hicat_test/inputs/hicat_lyot_mask_1944_gy_0.fits', 'alignment_tolerance': 6, 'num_lyot_stops': 9}, 
	'image': {'contrast': 8, 'iwa': 3.75, 'owa': 15.0, 'num_wavelengths': 4, 'bandwidth': 0.1, 'resolution': 2}, 
	'method': {'force_no_x_mirror_symmetry': False, 'force_no_y_mirror_symmetry': False, 'force_no_hermitian_symmetry': False, 'starting_scale': 1, 'ending_scale': 1, 'edge_width_for_prior': 2, 'num_throughput_iterations': 2, 'initial_throughput_estimate': 1, 'maximize_planet_throughput': True}, 
	'solver': {'num_threads': 0, 'crossover': 0, 'method': 2}}


                     
	file_organization = {'survey_dir': 'hicat_test', 
					   'solution_dir': 'hicat_test/solutions', 
					   'analysis_dir': 'hicat_test/analysis', 
					   'drivers_dir': 'hicat_test/scripts', 
					   'log_dir': 'hicat_test/logs', 
					   'input_files_dir': 'hicat_test/inputs'}

	pup_fname = parameters['pupil']['filename']
	fpm_radius = parameters['focal_plane_mask']['radius']
	fpm_num_pix = parameters['focal_plane_mask']['num_pix']
	fpm_grayscale = parameters['focal_plane_mask']['grayscale']
	ls_fname = parameters['lyot_stop']['filename']
	ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
	ls_num_stops = parameters['lyot_stop']['num_lyot_stops']
	
	
	pupil = read_fits(pup_fname)
	apodizer = read_fits(solution_filename)

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
			#lyot_stops.update({'lyot_stop_pos_x':lyot_stop_pos_x,'lyot_stop_neg_x':lyot_stop_neg_x,'lyot_stop_pos_y':lyot_stop_pos_y,'lyot_stop_neg_y':lyot_stop_neg_y})

		if ls_num_stops in [8, 9]:
			lyot_stop_pos_x_pos_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1), ls_alignment_tolerance, 0).ravel()
			lyot_stop_pos_x_neg_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1), -ls_alignment_tolerance, 0).ravel()
			lyot_stop_neg_x_pos_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1), ls_alignment_tolerance, 0).ravel()
			lyot_stop_neg_x_neg_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1), -ls_alignment_tolerance, 0).ravel()

			lyot_stops.extend([lyot_stop_pos_x_pos_y, lyot_stop_pos_x_neg_y, lyot_stop_neg_x_pos_y, lyot_stop_neg_x_neg_y])
			#lyot_stops.update({'lyot_stop_pos_x_pos_y':lyot_stop_pos_x_pos_y,'lyot_stop_pos_x_neg_y':lyot_stop_pos_x_neg_y,''})

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
	iwa = parameters['image']['iwa']
	owa = parameters['image']['owa']
	contrast = parameters['image']['contrast']
	radius_fpm = parameters['focal_plane_mask']['radius']

	lyot_stop = lyot_stops[0]

	coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
	focal_grid = make_focal_grid(pupil.grid, 8, owa * 1.2)
	prop = FraunhoferPropagator(pupil.grid, focal_grid)

	wf = Wavefront(pupil * apodizer)
	img = prop(coro(wf)).intensity
	img_ref = prop(Wavefront(apodizer * lyot_stop)).intensity

	plt.title('Monochromatic normalized irradiance')
	imshow_field(np.log10(img / img_ref.max()), vmin=-contrast-1, vmax=-contrast+4, cmap='inferno')
	plt.colorbar()
	plt.axis('off')
	if pdf is not None:
		pdf.savefig()
		plt.close()
	else:
		plt.show()

	r, profile, std_profile, n_profile = radial_profile(img / img_ref.max(), 0.2)

	plt.title('Monochromatic normalized irradiance (radial average)')
	plt.plot(r, profile)
	plt.axvline(iwa, color=colors.red)
	plt.axvline(owa, color=colors.red)
	plt.axvline(radius_fpm, color='k')
	plt.axhline(10**(-contrast), xmin=0, xmax=owa*1.2, linewidth=1, color='k', linestyle='--')
	
	plt.yscale('log')
	plt.xlim(0, owa*1.2)
	plt.ylim(5e-12, 2e-5)
	plt.ylabel('Normalized irradiance')
	plt.xlabel(r'Angular separation ($\lambda_0/D$)')
	
	if pdf is not None:
		pdf.savefig()
		plt.close()
	else:
		plt.show()

	return {'normalized_irradiance_image': img / img_ref.max(), 'normalized_irradiance_radial': (r, profile, std_profile, n_profile)}

def analyze_max_throughput(solution_filename, pdf=None):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	lyot_stop = lyot_stops[0]

	maximum_integrated_throughput = ((pupil * lyot_stop * apodizer).sum() / (pupil * lyot_stop).sum())**2
		
	return {'maximum_integrated_throughput': maximum_integrated_throughput}

def analyze_offaxis_throughput(solution_filename, pdf=None):
	pass

def analyze_summary(solution_filename, pdf=None):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	lyot_stop = lyot_stops[0]
	owa = parameters['image']['owa']
	bandwidth = parameters['image']['bandwidth']
	radius_fpm = parameters['focal_plane_mask']['radius']
	contrast = parameters['image']['contrast']

	coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
	coro_without_lyot = LyotCoronagraph(pupil.grid, focal_plane_mask)
	focal_grid = make_focal_grid(pupil.grid, 8, owa * 1.2)
	prop = FraunhoferPropagator(pupil.grid, focal_grid)
	wavelengths = np.linspace(-bandwidth / 2, bandwidth / 2, 11) + 1

	focal_plane_mask_large = 1 - circular_aperture(2 * radius_fpm)(focal_grid)

	img = 0
	img_ref = 0
	img_foc = 0
	lyot = 0

	for wl in wavelengths:
		wf = Wavefront(pupil * apodizer, wl)
		img += prop(coro(wf)).intensity
		img_foc += prop(wf).intensity
		img_ref += prop(Wavefront(pupil * lyot_stop)).intensity
		lyot += coro_without_lyot(wf).intensity
	
	font = {'family' : 'DejaVu Sans', 'weight' : 'medium', 'size'   : 10}
	matplotlib.rc('font', **font)
	
	#apodizer and telescope aperture
	plt.subplot(2,3,1)
	imshow_field(pupil * apodizer, cmap='Greys_r')
	plt.title('Apodizer and Telescope Aperture')
	plt.axis('off')
	
	#image plane
	plt.subplot(2,3,2)
	im = imshow_field(np.log10(img_foc / img_foc.max()), vmin=-5, vmax=0, cmap='inferno')
	plt.title('Image plane')
	plt.colorbar(im,fraction=0.046, pad=0.04)
	plt.axis('off')
	
	#image plane masked by focal plane mask
	plt.subplot(2,3,3)
	im = imshow_field(np.log10(img_foc / img_foc.max() * (1e-20 + focal_plane_mask_large)), vmin=-5, vmax=0, cmap='inferno')
	plt.title('Image plane w/FPM')
	plt.colorbar(im,fraction=0.046, pad=0.04)
	plt.axis('off')
	
	#lyot plane
	plt.subplot(2,3,4)
	im = imshow_field(np.log10(lyot / lyot.max()), vmin=-3, vmax=0, cmap='inferno')
	plt.title('Lyot plane')
	plt.colorbar(im,fraction=0.046, pad=0.04)
	plt.axis('off')
	
	#lyot plane masked by lyot stop
	plt.subplot(2,3,5)
	im = imshow_field(np.log10(lyot / lyot.max() * (1e-20 + lyot_stop)), vmin=-3, vmax=0, cmap='inferno')
	plt.title('Lyot plane w/lyot stop')
	plt.colorbar(im,fraction=0.046, pad=0.04)
	plt.axis('off')
	
	#final image plane
	plt.subplot(2,3,6)
	im = imshow_field(np.log10(img / img_ref.max()), vmin=-contrast-1, vmax=-contrast+4, cmap='inferno')
	plt.title('Final image plane')
	plt.colorbar(im,fraction=0.046, pad=0.04)
	plt.axis('off')
	
	if pdf is not None:
		pdf.savefig()
		plt.close()
	else:
		plt.show()

	return {}

def analyze_lyot_robustness(solution_filename, pdf=None):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	
	num_pix = parameters['pupil']['N'] #px
	iwa = parameters['image']['iwa']
	owa = parameters['image']['owa']
	bandwidth = parameters['image']['bandwidth']
	radius_fpm = parameters['focal_plane_mask']['radius']
	contrast = parameters['image']['contrast']
	
	ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
	
	wavelengths = np.linspace(-bandwidth / 2, bandwidth / 2, 11) + 1
	
	focal_grid = make_focal_grid(pupil.grid, 8, owa * 1.2)
	prop = FraunhoferPropagator(pupil.grid, focal_grid)
	
	lyot_stop = lyot_stops[0]
	
	dxs = np.array([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])
	#dxs = np.array(range(-ls_alignment_tolerance-1,+ls_alignment_tolerance+1,2))
	
		
	dither_grid = CartesianGrid(SeparatedCoords((dxs, dxs)))

	E = []
	mean_intensity = []
	plt.figure(figsize=(8,8))
	
	for i, (dx, dy) in enumerate(dither_grid.points):
		shifted_lyot_stop = np.roll(np.roll(lyot_stop.shaped, dx, 1), dy, 0).ravel()
		coro = LyotCoronagraph(pupil.grid, focal_plane_mask, shifted_lyot_stop)
		
		img = 0
		img_ref = 0
		
		for wl in wavelengths:
			wf = Wavefront(pupil * apodizer, wl)
			img += prop(coro(wf)).intensity
			img_ref += prop(Wavefront(pupil * lyot_stop)).intensity

		x = i % len(dxs)
		y = i // len(dxs)
	
		plt.subplot(len(dxs), len(dxs), x + (len(dxs) - y - 1) * len(dxs) + 1)
		
		imshow_field(np.log10(img / img_ref.max()), vmin=-contrast-1, vmax=-contrast+4, cmap='inferno')
		
		frame1 = plt.gca()
		frame1.axes.xaxis.set_ticklabels([])
		frame1.axes.yaxis.set_ticklabels([])

	plt.subplots_adjust(top=0.98, bottom=0.02, left=0.02, right=0.98, wspace=0, hspace=0)
	
	
	if pdf is not None:
		pdf.savefig()
		plt.close()
	else:
		plt.show()
	
	return {}

if __name__ == '__main__':
		
		identifier = 'HICAT_1944_LS9_6pix_telserv3'
		
		solution_filename = '/user/kstlaurent/git/progressive_refinement_coronagraphy/hicat_test/solutions/'+identifier+'.fits'
		
		pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
		
		analysis_summary_filename = os.path.join(file_organization['analysis_dir'], identifier + '.pdf')
		analysis_metrics_filename = os.path.join(file_organization['analysis_dir'], identifier + '.asdf')

		metrics = {}

		with PdfPages(analysis_summary_filename) as pdf:
			
			res = analyze_contrast_monochromatic(solution_filename, pdf)
			metrics.update(res)
			
			res = analyze_summary(solution_filename, pdf)
			metrics.update(res)
			
			res = analyze_lyot_robustness(solution_filename, pdf)
			metrics.update(res)
			
			res = analyze_max_throughput(solution_filename, pdf)
			metrics.update(res)


		# Write out metrics to a file
		f = asdf.AsdfFile(metrics)
		f.write_to(analysis_metrics_filename)

'''

def analyze_calculate_throughput(solution_filename):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	




	
'''
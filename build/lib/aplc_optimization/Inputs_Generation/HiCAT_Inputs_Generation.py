
from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
import os

def HiCAT_inputs_gen(input_files_dict):
	"""Generate HiCAT aperture and Lyot stop(s) input mask files for the APLC coronagraph.

    Parameters
    ----------
    input_files_dict: dict
        A dictionary of input mask parameters.

    Returns
    -------
    pup_filename: str
    	The name of the generated telescope aperture FITS file.
    ls_filename: list
    	The name(s) of the generated Lyot stop FITS file(s).
	"""

	# Contour
	p3_apodizer_size = 19.725e-3 # m
	gamma_31 = 1.008
	gamma_51 = 0.979
	
	filepath    = input_files_dict['directory']
	N           = input_files_dict['N']
	ap_spid     = input_files_dict['aperture']['ap_spid']
	ap_gap      = input_files_dict['aperture']['ap_gap']
	ap_grey		= input_files_dict['aperture']['ap_grey'] 
	pup_diam	= input_files_dict['aperture']['pup_diam']   	
	
	ls_spid		= input_files_dict['lyot_stop']['ls_spid'] 
	ls_grey		= input_files_dict['lyot_stop']['ls_grey'] 
	LS_OD       = input_files_dict['lyot_stop']['LS_OD']
	LS_ID       = input_files_dict['lyot_stop']['LS_ID']	
	
	if ap_grey:
		oversamp    = 4
	else:
		oversamp    = 1
	
	pup_filename = filepath+'ApodMask_HiCAT_ovsamp{0:02d}_N{1:04d}_Spider{2:}_Gap{3:}.fits'.format(oversamp, N,ap_spid,ap_gap)
		
	grid                        = make_uniform_grid(N, [p3_apodizer_size, p3_apodizer_size])

	print("test statement")
	print(Path().parent.resolve())
	config = Path('masks/'+pup_filename)
	if config.is_file():
		print('{0:s} exists'.format('masks/'+pup_filename))
	else:
		HiCAT_ap, header           = make_a_hicat_aperture(with_spiders=ap_spid, with_segment_gaps=ap_gap, return_header=True)
		
		pupil                      = evaluate_supersampled(HiCAT_ap,grid,oversamp)
	
		
		hdr = fits.Header()
		hdr.set('TELESCOP', header['TELESCOP'])
		hdr.set('P3_APOD', header['P3_APOD'],'m: p3 apodizer size')
		hdr.set('P3_CENT_SEG',header['P3_CENT_SEG'],'m: p3 apodizer mask central segment')
		
		if ap_gap:
			hdr.set('P3_GAP',header['P3_GAP'],'m: actual gap size between segments in P3')
			hdr.set('P3_GAP_OVER',header['P3_GAP_OVER'],'apod gap oversize factor wrt irisao ')
			
		if ap_spid:
			hdr.set('P3_STRUT',header['P3_STRUT'],'m: P3 apodizer mask spiders thickness')
		
		hdr.set('PROV',header['PROV'])
	
		fits.writeto('masks/'+pup_filename, pupil.shaped, hdr,overwrite=True)
		print('{0:s} has been written to file'.format('masks/'+pup_filename))
		
	ls_filenames = []
	
	if ls_grey:
		oversamp    = 4
	else:
		oversamp    = 1
	
	for ls_od in LS_OD:
		for ls_id in LS_ID:

			ls_filename  = filepath+'LS_HiCAT_ID{0:04d}_OD{1:04d}_Spider{2:}_ovsamp{3:02d}_N{4:04d}.fits'.format(int((ls_id)*1000),\
                                                                                    								int((ls_od)*1000), \
                                                                                        							ls_spid,oversamp, N)

			config = Path('masks/'+ls_filename)
			if config.is_file():
				print('{0:s} exists'.format('masks/'+ls_filename))
				
			else:
							
				HICAT_ls, ls_header = make_a_hicat_lyot_stop(lyot_inner=ls_id, lyot_outer=ls_od, with_spiders=ls_spid, return_header=True)

				magnification_tolerance = 0.000  # not investigated yet, in case we want to add robustness to magnification later
				nlyotstops = 1  # for the magnification robustness


				lyot_stop = evaluate_supersampled(HICAT_ls, grid.scaled(gamma_51 / gamma_31), oversamp)
				
				hdr = fits.Header()
				ls_header['OVERSAMP'] = oversamp
				  		
				hdr.set('TELESCOP', ls_header['TELESCOP'])
				hdr.set('P3_APOD', ls_header['P3_APOD'],'m: p3 apodizer size')
				hdr.set('P3_CENT_SEG',ls_header['P3_CENT_SEG'],'m: p3 apodizer mask central segment')
				hdr.set('LS_CENT', ls_header['LS_CENT'], 'm p5 ls mask central segment size')
				hdr.set('LS_SIZE', ls_header['LS_SIZE'], 'm p5 ls size')
				
				if ls_spid:
					hdr.set('P5_STRUT',ls_header['P5_STRUT'],'m: P5 ls mask spiders thickness')
				
				fits.writeto('masks/'+ls_filename, lyot_stop.shaped, hdr,overwrite=True)
				print('{0:s} has been written to file'.format('masks/'+ls_filename))
				
				
			ls_filenames.append(ls_filename)
			#remove duplicates
			s = set(ls_filenames)
			ls_filenames = list(s)
				
	
		return pup_filename, ls_filenames


def make_a_hicat_lyot_stop(normalized=False, with_spiders=True, lyot_inner=0.2, lyot_outer=0.9,
						 return_header=False):
	'''Make a HiCAT Lyot stop.

	Parameters
	----------
	normalized : boolean
		If this is True, the outer diameter will be scaled to 1. Otherwise, the
		diameter of the pupil will be 15.0 meters.
	with_spiders : boolean
		Include the secondary mirror support structure in the aperture.
	inner_diameter_fraction : scalar
		The fractional size of the circular central obstruction as fraction of the pupil diameter.
	outer_diameter_fraction : scalar
		The fractional size of the circular outer edge as fraction of the pupil diameter.
	return_header : boolean
		If this is True, a header will be returned giving all important values for the
		created aperture for reference.

	Returns
	-------
	lyot_stop : Field generator
		A field generator for the Lyot stop.
	header : dict
		A dictionary containing all important values for the created aperture. Only returned
		if `return_header` is True.
	'''
	gamma_21 = 0.423
	gamma_31 = 1.008
	gamma_51 = 0.979

	p2_irisao_segment_size = 1.4e-3  # m (note: point to point)
	p2_irisao_segment_side_length = p2_irisao_segment_size / 2
	p2_irisao_inscribed_circle_size = 10 * p2_irisao_segment_side_length
	p2_irisao_flat_to_flat_size = 14 * np.sqrt(3) / 2 * p2_irisao_segment_side_length
	p2_irisao_circumscribed_circle_size = np.sqrt(p2_irisao_flat_to_flat_size ** 2 + p2_irisao_segment_side_length ** 2)

	p3_apodizer_mask_central_segment_size = 3.950e-3  # m
	p3_apodizer_size = 19.725e-3  # m
	p5_apodizer_size = p3_apodizer_size * gamma_51 / gamma_31

	p5_lyot_stop_size = lyot_outer # m
	p5_irisao_inscribed_circle_size = p2_irisao_inscribed_circle_size * gamma_51 / gamma_21
	lyot_stop_mask_undersize_contour_wrt_inscribed_circle = p5_lyot_stop_size / p5_irisao_inscribed_circle_size

	p5_irisao_flat_to_flat_size = p2_irisao_flat_to_flat_size * gamma_51 / gamma_21
	p5_irisao_circumscribed_circle_size = p2_irisao_circumscribed_circle_size * gamma_51 / gamma_21

	# Central segment
	p5_lyot_stop_mask_central_segment_size = lyot_inner  # m
	p5_apodizer_mask_central_segment_size = p3_apodizer_mask_central_segment_size * gamma_51 / gamma_31

	p5_irisao_segment_size = p2_irisao_segment_size * gamma_51 / gamma_21
	lyot_stop_mask_central_segment_oversize_factor_wrt_apodizer_mask = p5_lyot_stop_mask_central_segment_size / p5_apodizer_mask_central_segment_size

	# Spiders
	p5_lyot_stop_mask_spiders_thickness = 0.700e-3  # m
	lyot_stop_mask_spiders_thickness_ratio = p5_lyot_stop_mask_spiders_thickness / p5_irisao_circumscribed_circle_size

	if normalized:
		p5_lyot_stop_size /= p5_apodizer_size
		p5_lyot_stop_mask_central_segment_size /= p5_apodizer_size
		p5_lyot_stop_mask_spiders_thickness /= p5_apodizer_size

	central_obscuration = circular_aperture(p5_lyot_stop_mask_central_segment_size)
	outer_diameter = circular_aperture(p5_lyot_stop_size)

	header = {'TELESCOP': 'HiCAT', 'P3_APOD': p3_apodizer_size, 'P3_CENT_SEG': p3_apodizer_mask_central_segment_size,
			  'LS_CENT': p5_lyot_stop_mask_central_segment_size, 'LS_SIZE': p5_lyot_stop_size,
			  'P5_STRUT': p5_lyot_stop_mask_spiders_thickness}

	if with_spiders:
		spider1 = make_spider_infinite([0, 0], 60, p5_lyot_stop_mask_spiders_thickness)
		spider2 = make_spider_infinite([0, 0], 120, p5_lyot_stop_mask_spiders_thickness)
		spider3 = make_spider_infinite([0, 0], -60, p5_lyot_stop_mask_spiders_thickness)
		spider4 = make_spider_infinite([0, 0], -120, p5_lyot_stop_mask_spiders_thickness)

	def func(grid):
		res = outer_diameter(grid) - central_obscuration(grid)

		if with_spiders:
			res *= spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid)

		return Field(res, grid)

	if return_header:
		return func, header
	else:
		return func


def make_a_hicat_aperture(normalized=False, with_spiders=True, with_segment_gaps=True, return_header=False,
						return_segments=False):
	'''Make the HiCAT P3 apodizer mask

	Parameters
	----------
	normalized : boolean
		If this is True, the outer diameter will be scaled to 1. Otherwise, the
		diameter of the pupil will be 15.0 meters.
	with_spiders : boolean
		Include the secondary mirror support structure in the aperture.
	with_segment_gaps : boolean
		Include the gaps between individual segments in the aperture.
	return_header : boolean
		If this is True , a header will be returned giving all important values for the
		created aperture for reference.
	return_segments : boolean
		If this is True, the segments will also be returned as a list of Field generators.

	Returns
	-------
	aperture : Field generator
		The HiCAT aperture.
	header : dict
		A dictionary containing all important values for the created aperture. Only returned
		if `return_header` is True.
	segments : list of Field generators
		The segments. Only returned when `return_segments` is True.
	'''
	gamma_21 = 0.423
	gamma_31 = 1.008

	# P2 - Iris AO
	p2_irisao_segment_size = 1.4e-3  # m (note: point to point)
	p2_irisao_segment_side_length = p2_irisao_segment_size / 2
	p2_irisao_segment_gap_size = 12e-6  # m

	p2_irisao_distance_between_segments = p2_irisao_segment_side_length * np.sqrt(3)
	p2_irisao_segment_circumdiameter = (2 * p2_irisao_segment_side_length) - (
				2 / np.sqrt(3)) * p2_irisao_segment_gap_size

	# P1 - Pupil Mask
	# Central segment
	p1_pupil_mask_central_segment_size = 3.600e-3  # m

	# P3 - Apodizer
	# Contour
	p3_apodizer_size = 19.725e-3  # m
	p2_apodizer_size = p3_apodizer_size * gamma_21 / gamma_31

	# Gap
	p3_apodizer_mask_gap_size = 0.090e-3  # m
	p3_irisao_segment_gap_size = p2_irisao_segment_gap_size * gamma_31 / gamma_21
	apodizer_mask_gap_oversize_factor_wrt_irisao = p3_apodizer_mask_gap_size / p3_irisao_segment_gap_size

	# Central segment
	p3_apodizer_mask_central_segment_size = 3.950e-3  # m
	p3_pupil_mask_central_segment_size = p1_pupil_mask_central_segment_size * gamma_31
	apodizer_mask_central_segment_oversize_factor_wrt_pupil_mask = p3_apodizer_mask_central_segment_size / p3_pupil_mask_central_segment_size
	p3_irisao_segment_size = p2_irisao_segment_size * gamma_31 / gamma_21

	# Spiders
	p3_apodizer_mask_spiders_thickness = 0.350e-3  # m

	header = {'TELESCOP': 'HiCAT', 'P3_APOD': p3_apodizer_size, 'P3_CENT_SEG': p3_apodizer_mask_central_segment_size,
			  'P3_GAP': p3_apodizer_mask_gap_size, 'P3_GAP_OVER': apodizer_mask_gap_oversize_factor_wrt_irisao,
			  'P3_STRUT': p3_apodizer_mask_spiders_thickness, 'PROV': 'HiCAT spreadsheet'}

	p3_irisao_segment_circumdiameter = p2_irisao_segment_circumdiameter * gamma_31 / gamma_21
	p3_irisao_distance_between_segments = p2_irisao_distance_between_segments * gamma_31 / gamma_21
	p3_apodizer_segment_circumdiameter = p3_irisao_segment_circumdiameter + (
				-p3_apodizer_mask_gap_size + p3_irisao_segment_gap_size) * (2 / np.sqrt(3))

	if normalized:
		p3_apodizer_segment_circumdiameter /= p3_apodizer_size
		p3_irisao_distance_between_segments /= p3_apodizer_size
		p3_apodizer_mask_central_segment_size /= p3_apodizer_size
		p3_apodizer_mask_spiders_thickness /= p3_apodizer_size
		p3_apodizer_size = 1

	segment = hexagonal_aperture(p3_apodizer_segment_circumdiameter, np.pi / 2)
	segment_positions = make_hexagonal_grid(p3_irisao_distance_between_segments, 3, False)
	segmentation = make_segmented_aperture(segment, segment_positions, return_segments=return_segments)

	if return_segments:
		segmentation, segments = segmentation

	segment = hexagonal_aperture(p3_apodizer_size / 7 / np.sqrt(3) * 2, np.pi / 2)
	distance_between_segments = p3_apodizer_size / 7
	segment_positions = make_hexagonal_grid(distance_between_segments, 3)
	contour = make_segmented_aperture(segment, segment_positions)

	central_segment = hexagonal_aperture(p3_apodizer_mask_central_segment_size, np.pi / 2)

	if return_segments:
		# Use function to return the lambda, to avoid incorrect binding of variables
		def segment_obstructed(segment):
			return lambda grid: segment(grid) * (contour(grid) - central_segment(grid))

		segments = [segment_obstructed(s) for s in segments]

	if with_spiders:
		spider1 = make_spider_infinite([0, 0], 60, p3_apodizer_mask_spiders_thickness)
		spider2 = make_spider_infinite([0, 0], 120, p3_apodizer_mask_spiders_thickness)
		spider3 = make_spider_infinite([0, 0], -60, p3_apodizer_mask_spiders_thickness)
		spider4 = make_spider_infinite([0, 0], -120, p3_apodizer_mask_spiders_thickness)

	if with_spiders and return_segments:
		# Use function to return the lambda, to avoid incorrect binding of variables
		def segment_with_spider(segment):
			return lambda grid: segment(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid)

		segments = [segment_with_spider(s) for s in segments]

	def func(grid):
		res = contour(grid) - central_segment(grid)

		if with_segment_gaps:
			res *= segmentation(grid)

		if with_spiders:
			res *= spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid)

		return Field(res, grid)

	if return_header:
		if return_segments:
			return func, header, segments
		else:
			return func, header
	elif return_segments:
		return func, segments
	else:
		return func

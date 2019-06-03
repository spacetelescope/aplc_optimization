from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

def HiCAT_inputs_gen(input_files_dict):
	
	# Contour
	p3_apodizer_size = 19.725e-3 # m
	
	filepath    = input_files_dict['directory']
	N           = input_files_dict['N']
	ap_spid     = input_files_dict['aperture']['ap_spid']
	ap_gap      = input_files_dict['aperture']['ap_gap']
	ap_grey		= input_files_dict['aperture']['ap_grey']   	
	lyot_ref_diam = input_files_dict['lyot_stop']['lyot_ref_diam']
	LS_OD         = input_files_dict['lyot_stop']['LS_OD']
	LS_ID         = input_files_dict['lyot_stop']['LS_ID']	
	
	if ap_grey:
		oversamp    = 4
	else:
		oversamp    = 4
	
	pup_filename = filepath+'ApodMask_HiCAT_ovsamp{0:02d}_N{1:04d}_Spider{2:}_Gap{3:}.fits'.format(oversamp, N,ap_spid,ap_gap)
	
	grid                        = make_uniform_grid(N, [p3_apodizer_size, p3_apodizer_size])
	
	config = Path('masks/'+pup_filename)
	if config.is_file():
		print('{0:s} exists'.format('masks/'+pup_filename))
	else:
		HiCAT_ap, header           = make_hicat_aperture(normalized=True, with_spiders=ap_spid, with_segment_gaps=ap_gap)
		pupil                      = evaluate_supersampled(HiCAT_ap,grid,oversamp)
	
		
		hdr = fits.Header()
		hdr.set('TELESCOP', header['TELESCOP'])
		hdr.set('D_CIRC', header['D_CIRC'],'m: circumscribed diameter')
		hdr.set('D_SEG_CIRC',header['D_SEG_CIRC'],'m: segment circumscribed diameter')
		
		if ap_gap:
			hdr.set('SEG_GAP',header['SEG_GAP'],'m: actual gap size between segments')
		if ap_spid:
			hdr.set('STRUT_W',header['STRUT_W'],'m: actual support strut width')
		
		hdr.set('NORM',header['NORM'],'normalization keyword, OD scaled to 1 by Dcirc')
		hdr.set('SEG_TRAN',header['SEG_TRAN'],'The transmission for each of the segments')
		hdr.set('PROV',header['PROV'])
	
		fits.writeto('masks/'+pup_filename, pupil.shaped, hdr,overwrite=True)
		print('{0:s} has been written to file'.format('masks/'+pup_filename))
		
	ls_filenames = []
	
	for ls_od in LS_OD:
		for ls_id in LS_ID:

			if ap_spid:
				ls_filename  = filepath+'LS_HiCAT_ID{0:04d}_OD{1:04d}_Spider{2:}_ovsamp{3:02d}_N{4:04d}.fits'.format(int((ls_id/lyot_ref_diam)*1000),\
                                                                                    								int((ls_od/lyot_ref_diam)*1000), \
                                                                                        							ap_spid,oversamp, N)
            
			ls_filename  = filepath+'LS_HiCAT_ID{0:04d}_OD{1:04d}_ovsamp{2:02d}_N{2:04d}.fits'.format(int((ls_id/lyot_ref_diam)*1000),\
                                                                                        int((ls_od/lyot_ref_diam)*1000), \
                                                                                        oversamp, N)

			config = Path('masks/'+ls_filename)
			if config.is_file():
				print('{0:s} exists'.format('masks/'+ls_filename))
				
			else:
							
				LUVOIR_ls, ls_header = make_hicat_lyot_stop(lyot_ref_diam, ls_id, ls_od,normalized=True, with_spiders=ap_spid)
				lyot_stop = evaluate_supersampled(LUVOIR_ls, grid, oversamp)
				
				hdr = fits.Header()
				ls_header['OVERSAMP'] = oversamp
				  		
				hdr.set('TELESCOP', ls_header['TELESCOP'])
				hdr.set('D_CIRC', ls_header['D_CIRC'],'m: circumscribed diameter')
				hdr.set('LS_REF_D',ls_header['LS_REF_D'],'m: used to reference given LS id and od')
				hdr.set('LS_ID', ls_header['LS_ID'], 'LS inner d, fraction of LS_REF_D')
				hdr.set('LS_OD', ls_header['LS_OD'], 'LS outer d, fraction of LS_REF_D')
				
				if ap_spid:
					hdr.set('STRUT_W',ls_header['LS_STRUT_W'],'m: actual support strut width')
	
				hdr.set('NORM',ls_header['NORM'],'normalization keyword, OD scaled to 1 by Dcirc')
				hdr.set('OVERSAMP',ls_header['OVERSAMP'],'oversampling factor, # grey levels')

				fits.writeto('masks/'+ls_filename, lyot_stop.shaped, hdr,overwrite=True)
				print('{0:s} has been written to file'.format('masks/'+ls_filename))
				
				
			ls_filenames.append(ls_filename)
			#remove duplicates
			s = set(ls_filenames)
			ls_filenames = list(s)
				
	
		return pup_filename, ls_filenames
	
		

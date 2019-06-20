from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

def HiCAT_inputs_gen(input_files_dict):
	
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
	
	config = Path('masks/'+pup_filename)
	if config.is_file():
		print('{0:s} exists'.format('masks/'+pup_filename))
	else:
		HiCAT_ap, header           = make_hicat_aperture(with_spiders=ap_spid, with_segment_gaps=ap_gap)
		
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

			ls_filename  = filepath+'LS_HiCAT_ID{0:04d}_OD{1:04d}_Spider{2:}_ovsamp{3:02d}_N{4:04d}.fits'.format(int((ls_id/pup_diam)*1000),\
                                                                                    								int((ls_od/pup_diam)*1000), \
                                                                                        							ls_spid,oversamp, N)
            
			config = Path('masks/'+ls_filename)
			if config.is_file():
				print('{0:s} exists'.format('masks/'+ls_filename))
				
			else:
							
				LUVOIR_ls, ls_header = make_hicat_lyot_stop(ls_id, ls_od, with_spiders=ls_spid)
				lyot_stop = evaluate_supersampled(LUVOIR_ls, grid.scaled(gamma_51 / gamma_31), oversamp)
				
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
	
		

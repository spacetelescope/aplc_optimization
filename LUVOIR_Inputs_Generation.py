from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

def LUVOIR_inputs_gen(input_files_dict):
	
	filepath    = input_files_dict['directory']
	N           = input_files_dict['N']
	oversamp    = input_files_dict['oversamp']
	gap_padding = input_files_dict['aperture']['seg_gap_pad']
	
	lyot_ref_diam = input_files_dict['lyot_stop']['lyot_ref_diam']
	LS_SPID       = input_files_dict['lyot_stop']['ls_spid']
	ls_spid_ov    = input_files_dict['lyot_stop']['ls_spid_ov']
	LS_OD         = input_files_dict['lyot_stop']['LS_OD']
	LS_ID         = input_files_dict['lyot_stop']['LS_ID']	
	
	pup_filename = filepath+'TelAp_LUVOIR_gap_pad{0:02d}_bw_ovsamp{1:02d}_N{2:04d}.fits'.format(gap_padding, oversamp, N)
	
	grid                        = make_pupil_grid(N)
	
	config = Path('masks/'+pup_filename)
	if config.is_file():
		print('{0:s} exists'.format('masks/'+pup_filename))
	else:
		LUVOIR_ap, header           = make_luvoir_a_aperture(gap_padding, header = True)
		pupil                       = evaluate_supersampled(LUVOIR_ap,grid,oversamp)
	
		header['EFF_GAP']  = header['GAP_PAD']*header['SEG_GAP']
	
		hdr = fits.Header()
		hdr.set('TELESCOP', header['TELESCOP'])
		hdr.set('D_CIRC', header['D_CIRC'],'m: circumscribed diameter')
		hdr.set('D_INSC', header['D_INSC'],'m: inscribed diameter')
		hdr.set('SEG_F2F',header['SEG_F2F_D'],'m: actual segment flat-to-flat diameter')
		hdr.set('SEG_GAP',header['SEG_GAP'],'m: actual gap size between segments')
		hdr.set('GAP_PAD',header['GAP_PAD'],'arbitrary multiplicative padding of gaps')
		hdr.set('EFF_GAP',header['EFF_GAP'],'m: effective gap size after padding')
		hdr.set('STRUT_W',header['STRUT_W'],'m: actual support strut width')
		hdr.set('STRUT_ST',header['STRUT_ST'],'m: lower spider starting point d from center')
		hdr.set('STRUT_AN',header['STRUT_AN'],'deg: angle lower spiders offset from vertical')
		hdr.set('NORM',header['NORM'],'normalization keyword, OD scaled to 1 by Dcirc')
		hdr.set('SEG_TRAN',header['SEG_TRAN'],'The transmission for each of the segments')
		hdr.set('EDGE','bw','black and white, or grey pixels')
		hdr.set('PROV',header['PROV'])
	
		fits.writeto('masks/'+pup_filename, pupil.shaped, hdr,overwrite=True)
		print('{0:s} has been written to file'.format('masks/'+pup_filename))
		
	ls_filenames = []
	
	
	
	for ls_od in LS_OD:
		for ls_id in LS_ID:
			#ls_id = (0.19/0.937)*ls_od
			
			#filename key for struts or no struts
			#black and white or grey pixels
			if ls_id == 0:
				LS_SPID = False
				
			if oversamp == 1:
				edge = 'bw'
			elif oversamp > 1:
				edge = 'gy'
			
			if LS_SPID:
				strut_key = 'struts'
			else:
				strut_key = 'no_struts'
				
				
				
			if LS_SPID:
				ls_filename  = filepath+'LS_LUVOIR_ID{0:04d}_OD{1:04d}_{2:s}_pad{3:02d}_{4:s}_ovsamp{5:d}_N{6:04d}.fits'.format(int(ls_id*1000),\
                                                                                                int(ls_od*1000), \
                                                                                                strut_key, ls_spid_ov,\
                                                                                                edge, oversamp, N)
			else:
				ls_filename  = filepath+'LS_LUVOIR_ID{0:04d}_OD{1:04d}_{2:s}_{3:s}_ovsamp{4:d}_N{5:04d}.fits'.format(int(ls_id*1000), \
                                                                                                int(ls_od*1000), \
                                                                                                strut_key, edge, \
                                                                                                oversamp, N)
			
			config = Path('masks/'+ls_filename)
			if config.is_file():
				print('{0:s} exists'.format('masks/'+ls_filename))
				
			else:
				
			
				LUVOIR_ls, ls_header = make_luvoir_a_lyot_stop(ls_id, ls_od, lyot_ref_diam, spid_oversize=ls_spid_ov, spiders=LS_SPID, header = True)
				lyot_stop = evaluate_supersampled(LUVOIR_ls, grid, oversamp)

				header.update(ls_header)
				hdr = fits.Header()
				header['OVERSAMP'] = oversamp
				header['EDGE']     = edge
  		
				hdr.set('TELESCOP', header['TELESCOP'])
				hdr.set('D_CIRC', header['D_CIRC'],'m: circumscribed diameter')
				hdr.set('D_INSC', header['D_INSC'],'m: inscribed diameter')
				hdr.set('LS_REF_D',header['LS_REF_D'],'m: used to reference given LS id and od')
				hdr.set('LS_ID', header['LS_ID'], 'LS inner d, fraction of LS_REF_D')
				hdr.set('LS_OD', header['LS_OD'], 'LS outer d, fraction of LS_REF_D')
			
				if LS_SPID:
					hdr.set('STRUT_W',header['STRUT_W'],'m: actual support strut width')
					hdr.set('STRUT_ST',header['STRUT_ST'],'m: lower spider starting point d from center')
					hdr.set('STRUT_AN',header['STRUT_AN'],'deg: angle lower spiders offset from vertical')
					hdr.set('STRUT_P',header['STRUT_P'], 'spider padding factor')
   
				hdr.set('NORM',header['NORM'],'normalization keyword, OD scaled to 1 by Dcirc')
				hdr.set('EDGE',header['EDGE'],'black and white, or grey pixels')
				hdr.set('OVERSAMP',header['OVERSAMP'],'oversampling factor, # grey levels')

				fits.writeto('masks/'+ls_filename, lyot_stop.shaped, hdr,overwrite=True)
				print('{0:s} has been written to file'.format('masks/'+ls_filename))
				
				
			ls_filenames.append(ls_filename)
			#remove duplicates
			s = set(ls_filenames)
			ls_filenames = list(s)
				
	
		return pup_filename, ls_filenames
	
		

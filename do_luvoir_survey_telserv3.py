from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

'''
n = 250


input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':1,\
					'aperture': {'seg_gap_pad':10}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.19], 'LS_OD':[0.95]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)


for ls_od in input_files_dict['lyot_stop']['LS_OD']:
	
	ls_id = (0.19/0.937)*ls_od

	ls_old = fits.getdata("/user/kstlaurent/LUVOIR_APLCs_SENT_TO_KRIST_09_2018/LS_full_luvoir_ann19D94_clear_N0250.fits")
	ls_new = fits.getdata("masks/LUVOIR/LS_LUVOIR_ID{0:04d}_OD{1:04d}_no_struts_bw_ovsamp1_N0250.fits".format(int(ls_id*1000),int(ls_od*1000)))


	diff = ls_old - ls_new

	fits.writeto("LS_diff_OD{0:04d}.fits".format(int(ls_od*1000)), diff, overwrite=True)


pup_old = fits.getdata("/user/kstlaurent/LUVOIR_APLCs_SENT_TO_KRIST_09_2018/TelAp_full_luvoirss100cobs1gap2_N0250.fits")
pup_new = fits.getdata("masks/"+pup_filename)
diff = pup_old - pup_new

fits.writeto("PUP_diff.fits", diff, overwrite=True)


exit()
survey_parameters_small = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.85},
                     'image': {'contrast':10,'iwa':3.75,'owa':12.0,'bandwidth':0.10,'num_wavelengths':4}}


luvoir_small = DesignParameterSurvey(PorAPLC, survey_parameters_small, 'surveys/luvoir_baseline_N{:04d}_small_telserv3/'.format(n), 'masks/')
luvoir_small.describe()

luvoir_small.write_drivers(True)
luvoir_small.run_optimizations(True)
luvoir_small.run_analyses(True)




input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':2}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0, 0.19], 'LS_OD':[0.937]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_middle = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':6.82},
                     'image': {'contrast':10,'iwa':6.72,'owa':21.50,'bandwidth':0.10,'num_wavelengths':4}}


luvoir_mid = DesignParameterSurvey(PorAPLC, survey_parameters_middle, 'surveys/luvoir_baseline_N{:04d}_middle_telserv3/'.format(n), 'masks/')
luvoir_mid.describe()

luvoir_mid.write_drivers(True)
luvoir_mid.run_optimizations(True)
luvoir_mid.run_analyses(True)
'''

n = 500

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':2}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0, 0.19], 'LS_OD':[0.937]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_large = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':12.14},
                     'image': {'contrast':10,'iwa':12.04,'owa':38.53,'bandwidth':0.10,'num_wavelengths':4}}


luvoir_lg = DesignParameterSurvey(PorAPLC, survey_parameters_large, 'surveys/luvoir_baseline_N{:04d}_large_telserv3/'.format(n), 'masks/')
luvoir_lg.describe()

luvoir_lg.write_drivers(True)
luvoir_lg.run_optimizations(True)
luvoir_lg.run_analyses(True)


#old LUVOIR baseline -> "10%" -> double check with Neil, fpm 0.10 margin?
#small:  fpm: , drkhl:  3.75  - 12.00
#medium: fpm: , drkhl:  6.72  - 21.50
#large:  fpm: , drkhl: 12.04 - 38.53

#new, from Kevin MCMC -> "18%"
#3.82, 3.00 - 12
#7.54, 6.72 - 21.5
#15.87, 15.05 - 60.20

#H = 0.56
#X = 3.2

#IWA_i = OWA_(i-1) * H
#OWA_i = IWA_i * X

#n=1000, oversamp = 3, seg_gap_pad=1
#n=500,  oversamp = 3, seg_gap_pad=2
#n=200,  oversamp = 4, seg_gap_pad=4

#print out size of darkzone array



					
from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen


n = 1000

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':1}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.195], 'LS_OD':[0.965]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_large = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':7.54},
                     'image': {'contrast':10,'iwa':6.72,'owa':21.50,'bandwidth':0.18,'num_wavelengths':8}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters_large, 'surveys/luvoir_baseline_N{:04d}_medium_AWS/'.format(n), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)


#old LUVOIR baseline -> "10%" -> double check with Neil, fpm 0.10 margin?
#small: 3.85  fpm: , drkhl:  3.75  - 12.00
#medium: 6.82 fpm: , drkhl:  6.72  - 21.50
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



					
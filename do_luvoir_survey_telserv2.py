from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

n = 200

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':4}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':[False,True],'ls_spid_ov':2,'LS_ID':[0, 0.19], 'LS_OD':[0.937]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_small = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.88},
                     'image': {'contrast':10,'iwa':3.00,'owa':12.0,'bandwidth':0.18,'num_wavelengths':8}}

luvoir_small = DesignParameterSurvey(PorAPLC, survey_parameters_small, 'surveys/luvoir_N{:04d}_small_telserv3/'.format(n), 'masks/')
luvoir_small.describe()

luvoir_small.write_drivers(True)
luvoir_small.run_optimizations(True)
luvoir_small.run_analyses(True)

survey_parameters_middle = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':7.60},
                     'image': {'contrast':10,'iwa':6.72,'owa':26.88,'bandwidth':0.18,'num_wavelengths':8}}


luvoir_mid = DesignParameterSurvey(PorAPLC, survey_parameters_middle, 'surveys/luvoir_N{:04d}_middle_telserv3/'.format(n), 'masks/')
luvoir_mid.describe()

luvoir_mid.write_drivers(True)
luvoir_mid.run_optimizations(True)
luvoir_mid.run_analyses(True)

survey_parameters_large = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':14.17},
                     'image': {'contrast':10,'iwa':15.05,'owa':60.2,'bandwidth':0.18,'num_wavelengths':8}}


luvoir_lg = DesignParameterSurvey(PorAPLC, survey_parameters_large, 'surveys/luvoir_N{:04d}_large_telserv3/'.format(n), 'masks/')
luvoir_lg.describe()

luvoir_lg.write_drivers(True)
luvoir_lg.run_optimizations(True)
luvoir_lg.run_analyses(True)

#small:  fpm: 3.88, drkhl:  3.00  - 12.00
#medium: fpm: 7.60, drkhl:  6.72  - 26.88
#large:  fpm: 14.17, drkhl: 15.05 - 60.20

#H = 0.56
#X = 4

#IWA_i = OWA_(i-1) * H
#OWA_i = IWA_i * X

#n=1000, oversamp = 3, seg_gap_pad=1
#n=500,  oversamp = 3, seg_gap_pad=2
#n=200, 



					
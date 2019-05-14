from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen


n = 100

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':4}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.12], 'LS_OD':[0.982]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_large = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.5, 'num_pix': 50},
                     'image': {'contrast':8,'iwa':3.4,'owa':8.0,'bandwidth':0.2,'num_wavelengths':2}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters_large, 'surveys/luvoir_N{:04d}_AWS_test/'.format(n), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)



					
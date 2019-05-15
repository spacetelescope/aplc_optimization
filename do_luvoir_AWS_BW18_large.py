from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen


n = 1000

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':1}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.195], 'LS_OD':[0.965]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':15.87, 'num_pix':400 , 'grayscale': True,},
                     'image': {'contrast':10,'iwa':15.05,'owa':60.20,'bandwidth':0.18,'num_wavelengths':8}, \
                     'method':{'starting_scale': 4}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/luvoir_BW18_large_N{:04d}_AWS/'.format(n), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)



					
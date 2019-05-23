from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from HiCAT_Inputs_Generation import HiCAT_inputs_gen


n = 486

FPM1 = 1
FPM2 = 2

pupil_diameter = 19.9e-3 #m

lyot_inner = 6.8e-3
lyot_outer = 15.9e-3


input_files_dict = {'directory':'HiCAT/', 'N':n, 'oversamp':4,\
					'aperture': {'ap_spid':True,'ap_gap':True}, \
					'lyot_stop':{'lyot_ref_diam':pupil_diameter,'LS_ID':[lyot_inner], 'LS_OD':[lyot_outer]}}

pup_filename, ls_filenames = HiCAT_inputs_gen(input_files_dict)

exit()

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':FPM1, 'num_pix':150 , 'grayscale': True,},
                     'image': {'contrast':8,'iwa':3.75,'owa':15.0,'bandwidth':0.10,'num_wavelengths':5}, \
                     'method':{'starting_scale': 1}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/hicat_test/'.format(n), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)



					
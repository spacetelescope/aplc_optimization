from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from HiCAT_Inputs_Generation import HiCAT_inputs_gen


n = 486

FPM1 = 8.543/2 # lambda_0/D radius
FPM2 = 6.271/2 # lambda_0/D radius

pupil_diameter = 19.9e-3 #m P3 apodizer size
lyot_inner = 6.8e-3 #m p5 lyot stop mask central segment size
lyot_outer = 15.9e-3 #m p5 lyot stop size


input_files_dict = {'directory':'HiCAT/', 'N':n,\
					'aperture': {'ap_spid':True,'ap_gap':True,'ap_grey':False,'pup_diam':pupil_diameter}, \
					'lyot_stop':{'ls_spid':True,'ls_grey':True,'LS_ID':[lyot_inner], 'LS_OD':[lyot_outer]}}

pup_filename, ls_filenames = HiCAT_inputs_gen(input_files_dict)

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames,'alignment_tolerance': 4, 'num_lyot_stops': 9}, \
                     'focal_plane_mask': {'radius':FPM1, 'num_pix': 80, 'grayscale': True,},
                     'image': {'contrast':8,'iwa':3.75,'owa':15.0,'bandwidth':0.10,'num_wavelengths':4}, \
                     'method':{'starting_scale': 1}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/hicat_test/'.format(n), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)



					
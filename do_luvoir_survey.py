from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen


input_files_dict = {'directory':'LUVOIR/', 'N':500, 'oversamp':2,\
					'aperture': {'seg_gap_pad':5}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.18,0.19,0.20], 'LS_OD':[0.92,0.93,0.937,0.94,0.95]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

n = input_files_dict['N']
id = input_files_dict['lyot_stop']['LS_ID']
od = input_files_dict['lyot_stop']['LS_OD']

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'LS_ID': id, 'LS_OD':od, 'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.82},
                     'image': {'contrast':10,'iwa':3,'owa':12,'bandwidth':0.18,'num_wavelengths':8}}

survey = DesignParameterSurvey(PorAPLC, survey_parameters, 'survey/', 'masks/')
survey.describe()

survey.write_drivers(True)
survey.run_optimizations(True)
survey.run_analyses(True)

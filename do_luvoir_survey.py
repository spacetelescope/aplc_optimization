from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen


input_files_dict = {'directory':'LUVOIR/', 'N':256, 'oversamp':2,\
					'aperture': {'seg_gap_pad':5}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.19], 'LS_OD':[0.937]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

print(pup_filename)

print(ls_filenames)

survey_parameters = {'pupil': {'filename': pup_filename}, 'lyot_stop': {'filename': ls_filenames}, 'image': {'owa': 12}}

survey = DesignParameterSurvey(PorAPLC, survey_parameters, 'survey/', 'masks/')
survey.describe()

survey.write_drivers(True)
survey.run_optimizations(True)
survey.run_analyses(True)

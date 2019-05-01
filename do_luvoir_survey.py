from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

n = 1000

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':1}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.19], 'LS_OD':[0.937]}}

input_files_dict_ls_spid = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':1}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':True,'ls_spid_ov':2,'LS_ID':[0.19], 'LS_OD':[0.937]}}
					

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)
pup_filename, ls_filenames_spid = LUVOIR_inputs_gen(input_files_dict_ls_spid)

ls_filenames.append(ls_filenames_spid[0])


survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.82},
                     'image': {'contrast':10,'iwa':3,'owa':12,'bandwidth':0.18,'num_wavelengths':8}}

luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'survey/luvoir/', 'masks/')
luvoir.describe()

#luvoir.write_drivers(True)
#luvoir.run_optimizations(True)
luvoir.run_analyses(True)

#3-12, 6.72 - 26.88, 
#H = 0.56
#X = ~3-4

#IWA_i = OWA_(i-1) * H
#OWA_i = IWA_i * X

#OWA_1 = 12
#IWA_2 = 12 * 0.56 = 6.72
#OWA_2 = 6.72 * 4 = 26.88

#IWA_3 = 26.88 * .56 = 15.05
#OWA_3 = 60.2



					
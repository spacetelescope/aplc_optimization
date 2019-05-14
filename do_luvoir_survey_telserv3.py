import time
from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

t0_total = time.time()

n = 1000


input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':3,\
					'aperture': {'seg_gap_pad':1}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.195], 'LS_OD':[0.965]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_small = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.82, 'num_pix': 250}, \
                     'image': {'contrast':10,'iwa':3.00,'owa':12.0,'bandwidth':0.18,'num_wavelengths':8}}


luvoir_small = DesignParameterSurvey(PorAPLC, survey_parameters_small, 'surveys/luvoir_18_N{:04d}_small_telserv3/'.format(n), 'masks/')
luvoir_small.describe()

luvoir_small.write_drivers(True)
t0 = time.time()
print('small optimization start time: {0:.2f}s'.format(t0-t0_total))
luvoir_small.run_optimizations(True)
t1 = time.time()
print('small optimization fin time: {0:.2f}s'.format(t1-t0_total))
print('small optimization run time: {0:.2f}s'.format(t1-t0))

luvoir_small.run_analyses(True)

n = 300

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':4}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.195], 'LS_OD':[0.965]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_middle = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':6.82, 'num_pix': 250}, \
                     'image': {'contrast':10,'iwa':6.72,'owa':21.50,'bandwidth':0.18,'num_wavelengths':8}}


luvoir_mid = DesignParameterSurvey(PorAPLC, survey_parameters_middle, 'surveys/luvoir_18_N{:04d}_middle_telserv3/'.format(n), 'masks/')
luvoir_mid.describe()

luvoir_mid.write_drivers(True)
t0 = time.time()
print('medium optimization start time: {0:.2f}s'.format(t0-t0_total))
luvoir_mid.run_optimizations(True)
t1 = time.time()
print('medium optimization fin time: {0:.2f}s'.format(t1-t0_total))
print('medium optimization run time: {0:.2f}s'.format(t1-t0))

luvoir_mid.run_analyses(True)

survey_parameters_large = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':15.87, 'num_pix': 250}, \
                     'image': {'contrast':10,'iwa':15.05,'owa':60.20,'bandwidth':0.18,'num_wavelengths':8}}


luvoir_lg = DesignParameterSurvey(PorAPLC, survey_parameters_large, 'surveys/luvoir_18_N{:04d}_large_telserv3/'.format(n), 'masks/')
luvoir_lg.describe()

luvoir_lg.write_drivers(True)
t0 = time.time()
print('large optimization start time: {0:.2f}s'.format(t0-t0_total))
luvoir_lg.run_optimizations(True)
t1 = time.time()
print('large optimization fin time: {0:.2f}s'.format(t1-t0_total))
print('large optimization time: {0:.2f}s'.format(t1-t0))
luvoir_lg.run_analyses(True)

t1_total = time.time()
print('total survey run time: {0:.2f}s'.format(t1_total-t0_total))

#H = 0.56
#X = 3.2

#IWA_i = OWA_(i-1) * H
#OWA_i = IWA_i * X

#n=1000, oversamp = 3, seg_gap_pad=1
#n=500,  oversamp = 3, seg_gap_pad=2
#n=200,  oversamp = 4, seg_gap_pad=4
#n=300,  oversamp = , seg_gap_pad=


					
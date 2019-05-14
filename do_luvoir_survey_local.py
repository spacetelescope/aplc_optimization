import time
from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

t0_total = time.time()

n = 256

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':4}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'LS_oversamp':1,'ls_spid':False,'ls_spid_ov':2,\
					'LS_ID':[0.12],'LS_OD':[0.982]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)


pup_filename = "TelAp.fits"
ls_filenames = "LS.fits"

survey_parameters_small = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.5,'num_pix': 80},
                     'image': {'contrast':10,'iwa':3.4,'owa':12.0,'bandwidth':0.10,'num_wavelengths':4}}

luvoir_small = DesignParameterSurvey(PorAPLC, survey_parameters_small, 'surveys/luvoir_baseline_Neil_inputs_small_local/'.format(n), 'masks/')
luvoir_small.describe()

luvoir_small.write_drivers(True)

t0 = time.time()
print('small optimization start time: {0:.2f}s'.format(t0-t0_total))
luvoir_small.run_optimizations(True)
t1 = time.time()
print('small optimization fin time: {0:.2f}s'.format(t1-t0_total))
print('optimization time: {0:.2f}s'.format(t1-t0))

luvoir_small.run_analyses(True)

t1_total = time.time()
print('total survey run time: {0:.2f}s'.format(t1_total-t0_total))

					
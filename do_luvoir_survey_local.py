import time
from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

n = 200
a = 0.937
b = 0.19
c = b/a

t0_total = time.time()

input_files_dict = {'directory':'LUVOIR/', 'N':n, 'oversamp':4,\
					'aperture': {'seg_gap_pad':4}, \
					'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,\
					'LS_OD':[a,a*1.015,a*1.03]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

survey_parameters_small = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':[3.72,3.82,3.92]},
                     'image': {'contrast':10,'iwa':[3,3.25,3.5],'owa':12.0,'bandwidth':[0.18,0.16],'num_wavelengths':8}}

luvoir_small = DesignParameterSurvey(PorAPLC, survey_parameters_small, 'surveys/luvoir_N{:04d}_small_local/'.format(n), 'masks/')
luvoir_small.describe()

luvoir_small.write_drivers(True)

t0 = time.time()
luvoir_small.run_optimizations(True)
t1 = time.time()
print('optimization time: {0:.2f}s'.format(t1-t0))

luvoir_small.run_analyses(True)

t1_total = time.time()
print('total survey run time: {0:.2f}s'.format(t1_total-t0_total))

#single, baseline design, 3-12, 18% nlam=8, run on local machine
#optimization time: 1889.58s
#total run time: 1892.74s

#old LUVOIR baseline -> "10%" -> double check with Neil, fpm 0.10 margin?
#small:  fpm: , drkhl:  3.75  - 12.00
#medium: fpm: , drkhl:  6.72  - 21.50
#large:  fpm: , drkhl: 12.04 - 38.53

#new, from Kevin MCMC -> "18%"
#3.82, 3.00 - 12
#7.54, 6.72 - 21.5
#15.87, 15.05 - 60.20

#H = 0.56
#X = 3.2

#IWA_i = OWA_(i-1) * H
#OWA_i = IWA_i * X

#n=1000, oversamp = 3, seg_gap_pad=1
#n=500,  oversamp = 3, seg_gap_pad=2
#n=200, 

#print out size of darkzone array

#,a*1.01,a*1.02,a*1.03,a*1.04,a*1.05,a*1.06,a*1.07,a*1.08,a*1.09,a*1.10

					
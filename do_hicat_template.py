from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from HiCAT_Inputs_Generation import HiCAT_inputs_gen

'''
DO NOT RUN THIS SCRIPT - It is the launcher template. 
Default parameters for 1000 pix BW10 small design

Workflow: 
- Make a copy
- rename the copy to designate survey you are running and the machine name you will be 
  running it on (do_luvoir_BW10_small_telserv3.py , for example, designates the BW10 small 
  design run on telserv3)
  
- set survey_name the same as survey name you use in the launcher file name
- Set the machine to be the name of the machine you're running on (need to write this to 
  grab machine name automatically)
  
- run the launcher
- after survey is done, move the launcher script into the survey directory
'''

n = 486 #number of pixels in input and final arrays
survey_name = "test" #survey name
machine = "telserv3" #machine the survey is run on. TODO: set this automatically

#Physical dimenstions of focal plane masks, pupil aperture, and lyot stops
FPM1 = 8.543/2 # lambda_0/D radius
FPM2 = 6.271/2 # lambda_0/D radius

pupil_diameter = 19.9e-3 #m P3 apodizer size
lyot_inner = 6.8e-3 #m p5 lyot stop mask central segment size
lyot_outer = 15.9e-3 #m p5 lyot stop size



#parameters for input files
input_files_dict = {'directory':'HiCAT/', 'N':n,\
					'aperture': {'ap_spid':True,'ap_gap':True,'ap_grey':False,'pup_diam':pupil_diameter}, \
					'lyot_stop':{'ls_spid':True,'ls_grey':True,'LS_ID':[lyot_inner], 'LS_OD':[lyot_outer]}}

pup_filename, ls_filenames = HiCAT_inputs_gen(input_files_dict)

#design parameters
#for multiple design parameters as a grid, input as list
#for multiple design parameters NOT as a grid, make multiple entries of below (as shown in commented out block)
#resulting apodizer are packaged with their inputs in a fits cube
survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames,'alignment_tolerance': 4, 'num_lyot_stops': 9}, \
                     'focal_plane_mask': {'radius':FPM1, 'num_pix': 80, 'grayscale': True,},
                     'image': {'contrast':8,'iwa':3.75,'owa':15.0,'bandwidth':0.10,'num_wavelengths':4}, \
                     'method':{'starting_scale': 1}}


hicat = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/hicat_{}_N{:04d}_{}/'.format(survey_name,n,machine), 'masks/')
hicat.describe()

hicat.write_drivers(True)
hicat.run_optimizations(True)
hicat.run_analyses(True)

'''
survey_parameters2 = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames,'alignment_tolerance': 4, 'num_lyot_stops': 9}, \
                     'focal_plane_mask': {'radius':FPM1, 'num_pix': 80, 'grayscale': True,},
                     'image': {'contrast':8,'iwa':3.75,'owa':15.0,'bandwidth':0.10,'num_wavelengths':4}, \
                     'method':{'starting_scale': 1}}


hicat = DesignParameterSurvey(PorAPLC, survey_parameters2, 'surveys/hicat_{}_N{:04d}_{}/'.format(survey_name,n,machine), 'masks/')
hicat.describe()

hicat.write_drivers(True)
hicat.run_optimizations(True)
hicat.run_analyses(True)
'''

					
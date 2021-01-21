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
instrument = 'hicat'  # instrument name
survey_name = "headertest"  # survey name
machine = "telserv3"  # machine the survey is run on.

nArray = 486  # number of pixels in input (TelAP, LS) and final (apodizer) arrays

# Aperture parameters
pupil_mask_size = 19.85e-3  # m: p1 pupil mask size
pupil_diameter = 19.725e-3  # m: p3 apodizer size
ap_spiders = True  # include secondary support mirror structure in the aperture
ap_gaps = True  # include the gaps between individual segments in the aperture
ap_grey = False  # grey pupil, else b&w

# Lyot stop parameters
lyot_inner = 6.8e-3  # m: P5 lyot stop mask central segment size
lyot_outer = 15.9e-3  # m: p5 lyot stop size
ls_spiders = True  # include secondary support mirror structure in the aperture
ls_grey = True  # grey lyot stop, else b&w

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {'directory': 'HiCAT/', 'N': nArray, \
                    'aperture': {'ap_spid': ap_spiders, 'ap_gap': ap_gaps, 'ap_grey': ap_grey,
                                 'pup_diam': pupil_diameter}, \
                    'lyot_stop': {'ls_spid': ls_spiders, 'ls_grey': ls_grey, 'LS_ID': [lyot_inner],
                                  'LS_OD': [lyot_outer]}}

# INPUT FILE GENERATION
pup_filename, ls_filenames = HiCAT_inputs_gen(input_files_dict)

'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple entries of below 
  (as shown in the commented block, at bottom of this script)
'''

# Lyot stop parameters
ls_alignment_tolerance = 1  # Lyot stop pixel shifts: tolerance to Lyot stop misalignments
ls_num = 1

# Focal plane mask parameters
FPM1 = 8.543 / 2  # lambda_0/D: focal plane mask 1 radius
FPM2 = 6.271 / 2  # lambda_0/D: focal plane mask 2 radius
nFPM = 80  # pixels
FPM_grey = True  # Focal plane mask is grey, otherwise bw

# Optimization parameters
contrast = 8  # 10-<value>: contrast goal in the dark zone of the coronagraphic image
IWA = 3.75  # lamda_0/D: inner edge of dark zone region in coronagraphic image
OWA = 15  # lambda_0/D: outer edge of dark zone region in coronagraphic image
spectral_bandwidth = 0.0  # fractional, dark zone bandpass
nLams = 1  # number of wavelengths spanning the design bandpass

# Optimization method
ending_scale = 1

# SURVEY PARAMETER DICTIONARY
survey_parameters = {'instrument': {'inst_name': instrument.upper()}, \
                     'pupil': {'N': nArray, 'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames, 'alignment_tolerance': ls_alignment_tolerance,
                                   'num_lyot_stops': ls_num}, \
                     'focal_plane_mask': {'radius': FPM1, 'num_pix': nFPM, 'grayscale': FPM_grey},
                     'image': {'contrast': contrast, 'iwa': IWA, 'owa': OWA, 'bandwidth': spectral_bandwidth,
                               'num_wavelengths': nLams}, \
                     'method': {'ending_scale': ending_scale}}

# RUN DESIGN SURVEY
hicat = DesignParameterSurvey(PorAPLC, survey_parameters,
                              'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, nArray, machine), 'masks/')
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



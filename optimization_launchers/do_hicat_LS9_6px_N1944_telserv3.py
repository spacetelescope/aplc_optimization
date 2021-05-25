from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.HiCAT_Inputs_Generation import HiCAT_inputs_gen

instrument = 'hicat'  # instrument name
survey_name = 'hicat_LS9_6px'  # survey name
machine = "telserv3"  # machine the survey is run on.

nArray = 1944  # number of pixels in input (TelAP, LS) and final (apodizer) arrays

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
input_files_dict = {'directory': 'HiCAT/', 'N': nArray,
                    'aperture': {'ap_spid': ap_spiders, 'ap_gap': ap_gaps, 'ap_grey': ap_grey, 'pup_diam': pupil_diameter},
                    'lyot_stop': {'ls_spid': ls_spiders, 'ls_grey': ls_grey, 'LS_ID': [lyot_inner], 'LS_OD': [lyot_outer]}}

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
ls_alignment_tolerance = 6  # Lyot stop pixel shifts: tolerance to Lyot stop misalignments
ls_num = 9 # Number of Lyot stops

# Focal plane mask parameters
FPM1 = 8.543 / 2  # lambda_0/D: focal plane mask 1 radius
FPM2 = 6.271 / 2  # lambda_0/D: focal plane mask 2 radius
nFPM = 80  # pixels
FPM_grey = True  # Focal plane mask is grey, otherwise bw

# Optimization parameters
contrast = 8  # 10-<value>: contrast goal in the dark zone of the coronagraphic image
IWA = 3.75  # lamda_0/D: inner edge of dark zone region in coronagraphic image
OWA = 15  # lambda_0/D: outer edge of dark zone region in coronagraphic image
spectral_bandwidth = 0.1  # fractional, dark zone bandpass
nLams = 1  # number of wavelengths spanning the design bandpass
resolution = 2.5 # px: spatial resolution of final coronagraphic image

# Optimization method
maximize_planet_throughput = False # Maximize throughput of apodizer
ending_scale = 4 # number of pixels per unit cell for the final solution.

# SURVEY PARAMETER DICTIONARY
survey_parameters = {'instrument': {'inst_name': instrument.upper()},
                     'pupil': {'N': nArray, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filenames, 'alignment_tolerance': ls_alignment_tolerance, 'num_lyot_stops': ls_num},
                     'focal_plane_mask': {'radius': FPM1, 'num_pix': nFPM, 'grayscale': FPM_grey},
                     'image': {'contrast': contrast, 'iwa': IWA, 'owa': OWA, 'bandwidth': spectral_bandwidth, 'num_wavelengths': nLams, 'resolution': resolution},
                     'method': {'maximize_planet_throughput': maximize_planet_throughput, 'ending_scale': ending_scale}}

# RUN DESIGN SURVEY
hicat = DesignParameterSurvey(APLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, nArray, machine),
                              'masks/')
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



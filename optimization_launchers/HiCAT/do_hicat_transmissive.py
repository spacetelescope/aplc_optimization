'''
DO NOT RUN THIS SCRIPT - It is the launcher template for HiCAT design surveys.

Workflow:
- Make a copy of this script.
- Rename the copy with the following naming schema: 'do_hicat_<survey_name>_<machine>.py', designating the name of the
  survey you are running and the machine you will be running it on.
- Run the launcher on the designated machine.
'''

import os

os.chdir('../..')
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.HiCAT_Inputs_Generation import HiCAT_inputs_gen

# Survey information
instrument = 'hicat'  # instrument name
survey_name = 'launcher_template'  # survey name
machine = 'local'  # machine the survey is run on.

N = 488  # number of pixels in input (TelAP, LS) and final (apodizer) arrays

'''
Input (aperture and Lyot stop) Array Parameters
-----------------------------------------------
N: int
    The number of pixels in input (TelAP, LS) and final (apodizer) arrays.
ap_spiders: bool
    Whether to include secondary support mirror structure in the aperture.
ap_gaps: bool
    Wether to include the gaps between individual segments in the aperture
ap_grey: bool
    Whether to model a grey pupil, else black and white.
pup_diam: float
    The HiCAT P1 pupil mask size in meters.
lS_ID: float
    The P5 lyot stop mask central segment size in meters.
lS_OD: float
    The P5 lyot stop size in meters.
ls_spid: bool
    Whether to include secondary support mirror structure in the aperture.
ls_grey: bool
    Whether to model a grey lyot stop, else black and white.
'''

# Aperture parameters
pup_diam = 19.52e-3  # m: p3 apodizer size
ap_spid = False
ap_gap = False
ap_grey = False

# Lyot stop parameters
LS_ID = 0
LS_OD = 15.0e-3
ls_spid = False
ls_grey = False

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {'directory': 'HiCAT/', 'N': N,
                    'aperture': {'ap_spid': ap_spid, 'ap_gap': ap_gap, 'ap_grey': ap_grey, 'pup_diam': pup_diam},
                    'lyot_stop': {'ls_spid': ls_spid, 'ls_grey': ls_grey, 'LS_ID': [LS_ID],
                                  'LS_OD': [LS_OD]}}

# INPUT FILE GENERATION
#pup_filename, ls_filenames = HiCAT_inputs_gen(input_files_dict)
pup_filename = 'HiCAT/hicat_apodizer_mask_488_bw.fits'
ls_filenames = 'HiCAT/hicat_lyot_mask_488_gy_0.fits'

'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple entries of below 
  (as shown in the commented block, at bottom of this script)

N: int
    The number of pixels in input (TelAP, LS) and final (apodizer) arrays.
alignment_tolerance: int
    The Lyot stop alignment tolerance in pixels. 
num_lyot_stops: int
    The number of translated Lyot stops to optimize for simultaneously. If 1, design will be non-robust.
radius: float
    The radius of the focal plane mask in lambda_0/D.
num_pix: float
    The number of pixels in the focal plane mask array.
grayscale: bool
    Whether the focal plane mask is gray, else black and white.
contrast: int
    The contrast goal in the dark zone of the coronagraphic image (negative exponent of 10).
iwa: float
    The inner edge of dark zone region in coronagraphic image in lambda_0/D.
owa: float
    The outer edge of dark zone region in coronagraphic image in lambda0/D.
bandwidth: float
    The spectral bandwidth over which to optimize (fractional).
num_wavelengths: int
    The number of wavelengths spanning the design bandpass.
resolution: 
    The image resolution.
starting_scale: int
    The number of pixels per unit cell for the initial solution. This is used for the adaptive algorithm. 
    It must be a power of 2 times the `ending_scale`.
ending_scale: int
    The number of pixels per unit cell for the final solution. If this is the same as `starting_scale`,
    the adaptive algorithm is essentially turned off.
'''

# Lyot stop parameters
alignment_tolerance = 1  # pixels
num_lyot_stops = 1

# Focal plane mask parameters
FPM1 = 8.525 / 2  # lambda_0/D: focal plane mask 1 radius
radius = FPM1
num_pix = 80
grayscale = True

# Optimization parameters
contrast = 8  # 10-<value>
iwa = 3.75  # lamda_0/D
owa = 16  # lambda_0/D
bandwidth = 0.15
num_wavelengths = 1
resolution = 2.5

# Optimization method
starting_scale = 4
ending_scale = 1

# SURVEY PARAMETER DICTIONARY
survey_parameters = {'instrument': {'inst_name': instrument.upper()},
                     'pupil': {'N': N, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filenames, 'alignment_tolerance': alignment_tolerance,
                                   'num_lyot_stops': num_lyot_stops},
                     'focal_plane_mask': {'radius': FPM1, 'num_pix': num_pix, 'grayscale': grayscale},
                     'image': {'contrast': contrast, 'iwa': iwa, 'owa': owa, 'bandwidth': bandwidth,
                               'num_wavelengths': num_wavelengths, 'resolution': resolution},
                     'method': {'starting_scale': starting_scale, 'ending_scale': ending_scale}}

# RUN DESIGN SURVEY
hicat = DesignParameterSurvey(APLC, survey_parameters,
                              'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, N, machine),
                              'masks/')
hicat.describe()

hicat.write_drivers(True)
hicat.run_optimizations(True)
hicat.run_analyses(True)

'''
Example for multiple design parameters NOT as a grid
########################################################

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

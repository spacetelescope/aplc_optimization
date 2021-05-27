'''
DO NOT RUN THIS SCRIPT - It is the launcher template for GPI design surveys.

Workflow:
- Make a copy of this script.
- rename the copy under the following naming schema: 'do_luvoir_<survey_name>_<machine>.py', designating the survey you
  are running and the name of the machine you will be running it on.
- run the launcher on the designated machine.
'''

import os

os.chdir('..')
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.GPI_Inputs_Generation import GPI_inputs_gen

# Survey information
N = 128  # number of pixels in input and final arrays
instrument = "gpi"  # instrument name
survey_name = "template_survey"  # survey name
machine = "telserv3"  # machine the survey is run on.

'''
Input (primary and Lyot stop) Mask Parameters
-----------------------------------------------
N: int
    The number of pixels in input (primary and lyot stop) arrays.
ap_spid: bool
    Whether to include the secondary supports primary (which are masked out in the Lyot planes regardless).
ap_sym: bool
    Whether to model the primary as precisely symmetric by treating all 4 struts as the same thickness.
satspots: bool
    Whether to include satellite spots grid.
ls_tabs: bool
    Whether to block out the bad actuators with tabs in the Lyot masks.
lyot_mask: string
    The name of the Lyot mask. Available lyot masks are: '080m12_02', '080m12_03', '080m12_04', '080m12_04_c', 
    '080m12_06', '080m12_06_03', '080mgit 12_07', '080m12_10', 'Open' and 'Blank'.
ls_spid: string
    Whether to include the secondary supports in the Lyot stop mask.
'''

# Aperture parameters
ap_spid = True
ap_sym = True
satspots = False

# Lyot stop parameters
ls_tabs = False
lyot_mask = '080m12_04'
ls_spid = True

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {'directory': 'GPI/', 'N': N,
                    'aperture': {'ap_spid': True, 'ap_sym': True},
                    'lyot_stop': {'lyot_mask': lyot_mask, 'ls_tabs': ls_tabs, 'ls_spid': ls_spid}}

# INPUT FILE GENERATION
pup_filename, ls_filename = GPI_inputs_gen(
    input_files_dict)  # generate input mask files according to the parameters defined above

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
ls_num: int
    The number of translated Lyot stops to optimize for simultaneously. If 1, design will be non-robust.
radius: float
    The radius of the focal plane mask in lambda_0/D.
num_pix: float
    The number of pixels in the focal plane mask array.
FPM_name: str
    Name of the focal plane mask. Available FPMs are: 'Y', 'J', 'H' or 'K1'
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
resolution: float
    The image resolution.
starting_scale: int
    The number of pixels per unit cell for the initial solution. This is used for the adaptive algorithm. 
    It must be a power of 2 times the `ending_scale`.
ending_scale: int
    The number of pixels per unit cell for the final solution. If this is the same as `starting_scale`,
    the adaptive algorithm is essentially turned off.
'''

# Focal plane mask parameters
nFPM = 80
FPM_name = 'K1'
radius = 3.476449131

# Optimization constraints
contrast = 8
IWA = 3
OWA = 22
bandwidth = 0.2
num_wavelengths = 1
resolution = 2

# Optimization method
starting_scale = 4
ending_scale = 1

# Robustness parameters
alignment_tolerance = 1
num_lyot_stops = 1

survey_parameters = {'instrument': {'inst_name': instrument.upper()},
                     'pupil': {'N': N, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filename, 'ls_tabs': True, 'alignment_tolerance': alignment_tolerance,
                                   'num_lyot_stops': num_lyot_stops},
                     'focal_plane_mask': {'FPM_name': FPM_name, 'radius': radius, 'num_pix': N},
                     'image': {'contrast': contrast, 'iwa': IWA, 'owa': OWA, 'bandwidth': bandwidth,
                               'num_wavelengths': num_wavelengths, 'resolution': resolution},
                     'method': {'starting_scale': starting_scale, 'ending_scale': ending_scale}}

# RUN DESIGN SURVEY
gpi = DesignParameterSurvey(APLC, survey_parameters,
                            'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, N, machine),
                            'masks/')
gpi.describe()

gpi.write_drivers(True)
gpi.run_optimizations(True)
gpi.run_analyses(True)

'''
Example for multiple design parameters NOT as a grid
########################################################

survey_parameters_2 = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':6.82, 'num_pix': 250, 'grayscale': True,},
                     'image': {'contrast':10,'iwa':6.72,'owa':23.72,'bandwidth':0.10,'num_wavelengths':5}, \
                     'method':{'starting_scale': 4}}

luvoir = DesignParameterSurvey(PorAPLC, survey_parameters_2, 'surveys/luvoir_{}_small_N{:04d}_{}/'.format(survey_name,n,machine), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)
'''

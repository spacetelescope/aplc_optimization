import os

os.chdir('..')
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.GPI_Inputs_Generation import GPI_inputs_gen

# Survey Information
N = 1168  # with 10µm pixels (so we have 11.68mm instead of 11.67mm diameter)
instrument = "gpi"  # instrument name
survey_name = "pathfinder_robust"  # survey name
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
ap_spid = True  # the secondary supports are masked out in the Lyot planes regardless
ap_sym = True  # symmetric supports
satspots = False  # whether to include satellite spots grid

# Lyot stop parameters
ls_tabs = False  # whether to block out the bad actuators with tabs in the Lyot masks
lyot_mask = '080m12_03'  # name of Lyot mask. - selected '080m12_03' because it has the thinnest spiders
# available lyot masks are: '080m12_02', '080m12_03', '080m12_04', '080m12_04_c',
# '080m12_06', '080m12_06_03', '080m12_07', '080m12_10', 'Open' and 'Blank'.
ls_spid = True

input_files_dict = {'directory': 'GPI/', 'N': N,
                    'aperture': {'ap_spid': True, 'ap_sym': True},
                    'lyot_stop': {'lyot_mask': lyot_mask, 'ls_tabs': ls_tabs, 'ls_spid': ls_spid}}

# GENERATE INPUT FILES
pup_filename, ls_filename = GPI_inputs_gen(input_files_dict)

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
    The inner edge of dark zone region in coronagraphic image in lambda_0/D (rho_i).
owa: float
    The outer edge of dark zone region in coronagraphic image in lambda0/D (rho_o).
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

# Focal plane mask parameters
num_pix = 80  # number of pixels across the focal plane mask
FPM_name = 'K1'  # name of FPM. Available FPMs are: 'Y', 'J', 'H' or 'K1'
radius = 3.476449131  # lambda_0/D: focal plane mask radius

# Optimization parameters
contrast = 8
IWA = 3
OWA = 22
bandwidth = 0.2
num_wavelengths = 3
resolution = 2

# Optimization method parameters
starting_scale = 4
ending_scale = 1

# Robustness parameters
alignment_tolerance = 1
num_lyot_stops = 9

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
                            'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, N, machine), 'masks/')
gpi.describe()

gpi.write_drivers(True)
gpi.run_optimizations(True)
gpi.run_analyses(True)

'''
DO NOT RUN THIS SCRIPT - this is the launcher template for LUVOIR design surveys.

Workflow:
- Make a copy of this script.
- rename the copy under the following naming schema: 'do_luvoir_<survey_name>_<machine>.py', designating the survey you
  are running and the name of the machine you will be running it on.
- Define the Survey Information, Input File and Design Survey Parameters, as desired.
- Run the script on the designated machine.
'''

import os

os.chdir('../..')
instrument = 'scda'  # do not edit
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.SCDA_Inputs_Generation import SCDA_inputs_gen

'''
Survey Information
------------------
survey_name: str
    The designated name of the design survey.
machine: str
    The name of the machine on which to run the design survey.
N: int
    The number of pixels in the input (TelAP, LS) and final (apodizer) arrays.
'''
survey_name = 'scda_hex_template'
machine = 'local'
N = 1000

'''
Input File Parameters
----------------------
Gap padding (seg_gap_pad) and grey levels (oversamp) are set according to number of input array pixels (nArray),  
configured in order to keep gap size as close to actual physical size of LUVOIR A, as possible. 

 - N = 1000, oversamp = 4, seg_gap_pad = 1
 - N = 500,  oversamp = 3, seg_gap_pad = 2
 - N = 300,  oversamp = 4, seg_gap_pad = 4
 - N = 200,  oversamp = 4, seg_gap_pad = 4
 - N = 100,  oversamp = 4, seg_gap_pad = 4

N: int
    The number of pixels in input (TelAP, LS) and final (apodizer) arrays
oversamp: int
    The oversampling factor (number of grey levels) - if set to 1 will return a bw pupil, for grey set to > 1.
gap_padding: int
    An arbitary padding of gap size to represent gaps on smaller arrays. This effectively makes the larger gaps larger 
    and the segments smaller to preserve the same segment pitch.
num_rings : int
    The number of rings of hexagons to include, not counting the central segment.
clipped : bool
    Remove corner segments to maximize the diameter of the inscribed circle.
ap_spid : bool
    Include the secondary mirror support structure in the aperture.
ap_obs : bool
    Include the secondary obstruction in the aperture.
lyot_ref_diam: float 
    The diameter used to reference the LS inner and outer diameter against.
ls_spid: bool 
    Include secondary support mirror structure in the aperture.
ls_spid_ov: int
    The factor by which to oversize the spiders compared to the LUVOIR-A aperture spiders.
LS_ID: float
    The Lyot stop inner diameter(s) relative to the `lyot_ref_diameter` (inscribed circle). This is re-normalized 
    against the circumscribed pupil in the `LUVOIR_inputs_gen` function.
LS_OD: float
    The Lyot stop outer diameter(s) relative to the `lyot_ref_diameter` (inscribed circle). This is re-normalized 
    against the circumscribed pupil in the `LUVOIR_inputs_gen` function.
'''
# Aperture parameters
pupil_diameter = 15.0  # m: actual LUVOIR A circumscribed pupil diameter
pupil_inscribed = 13.5  # m: actual LUVOIR A inscribed pupil diameter
oversamp = 4
seg_gap_pad = 1
num_rings = 10
clipped = True
ap_spid = False
ap_obs = False

# Lyot stop parameters
lyot_ref_diam = pupil_inscribed
ls_spid = False
ls_spid_ov = 2
LS_ID = 0.12
LS_OD = 0.982

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {'directory': 'SCDA/', 'N': N, 'oversamp': oversamp,
                    'aperture': {'num_rings': num_rings, 'clipped': clipped, 'ap_spid': ap_spid, 'ap_obs': ap_obs
                                 'seg_gap_pad': seg_gap_pad},
                    'lyot_stop': {'lyot_ref_diam': lyot_ref_diam, 'ls_spid': ls_spid, 'ls_spid_ov': ls_spid_ov,
                                  'LS_ID': [LS_ID], 'LS_OD': [LS_OD]}}

# INPUT FILE GENERATION
pup_filename, ls_filenames = SCDA_inputs_gen(input_files_dict)

'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple `survey_parameters` dictionaries
  (as shown in the commented block, at bottom of this script).

Parameters
----------  
radius: float  
    The radius of the FPM in lamda_0/D.
num_pix: float
    The number of pixels in the focal plane mask.
greyscale: bool
    Whether to model a grayscale focal plane mask, else black and white.
iwa: float
    The effective inner working angle (outer perimeter of the annular dark zone in coronagraphic image) in lam/D.
owa: float
    The effective outer working angle (inner perimeter of the annular dark zone in the coronagraphic image) in lam/D.
bandwidth: float
    The spectral bandwidth over which to optimize (fractional).
num_wavelengths: int
    The number of wavelengths spanning the bandpass.
contrast: int
    The contrast goal in the dark zone (exponent of 10). 
'''
# FPM Parameters
radius = 3.5
num_pix = 150
grayscale = True

# Optimization parameters (dark zone constraints)
iwa = 3.4
owa = 12.0
bandwidth = 0.1
num_wavelengths = 1
contrast = 10

# SURVEY PARAMS DICTIONARY
survey_parameters = {'pupil': {'N': N, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filenames},
                     'focal_plane_mask': {'radius': radius, 'num_pix': num_pix, 'grayscale': grayscale},
                     'image': {'contrast': contrast, 'iwa': iwa, 'owa': owa, 'bandwidth': bandwidth,
                               'num_wavelengths': num_wavelengths}}

# RUN DESIGN SURVEY
SCDA = DesignParameterSurvey(APLC, survey_parameters,
                               'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, N, machine), 'masks/')
SCDA.describe()

SCDA.write_drivers(False)
SCDA.run_optimizations(False)
SCDA.run_analyses(False)

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

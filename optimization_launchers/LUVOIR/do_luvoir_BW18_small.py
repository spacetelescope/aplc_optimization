import os

os.chdir('../..')
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

# Survey information
instrument = "luvoir"  # instrument name
survey_name = "BW18_small"  # survey name
machine = "local"  # machine the survey is run on.
N = 1000  # number of pixels in input (TelAP, LS) and final (apodizer) arrays

'''
Input (aperture and Lyot stop) Array Parameters
-----------------------------------------------
N: int
    The number of pixels in input (TelAP, LS) and final (apodizer) arrays
oversamp: int
    The oversampling factor (number of grey levels) - if set to 1 will return a bw pupil, for grey set to > 1.
gap_padding: int
    An arbitary padding of gap size to represent gaps on smaller arrays. This effectively makes the larger gaps larger 
    and the segments smaller to preserve the same segment pitch.
lyot_ref_diam: float 
    The diameter used to reference the LS inner and outer diameter against.
ls_spid: bool 
    Whether to include secondary support mirror structure in the aperture.
ls_spid_ov: int
    The factor by which to oversize the spiders compared to the LUVOIR-A aperture spiders.
LS_ID: float
    The Lyot stop inner diameter(s) as a fraction of `lyot_ref_diameter`
LS_OD: float
    The Lyot stop outer diameter as a fraction of `lyot_ref_diameter`.
'''

input_files_dict = {'directory': 'LUVOIR/', 'N': N, 'oversamp': 4,
                    'aperture': {'seg_gap_pad': 1},
                    'lyot_stop': {'lyot_ref_diam': 13.5, 'ls_spid': False, 'ls_spid_ov': 2, 'LS_ID': [0.195],
                                  'LS_OD': [0.965]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)

'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple survey_parameters dictionaries.

Parameters
----------
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

survey_parameters = {'pupil': {'N': N, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filenames},
                     'focal_plane_mask': {'radius': 3.82, 'num_pix': 150, 'grayscale': True, },
                     'image': {'contrast': 10, 'iwa': 3.00, 'owa': 12.00, 'bandwidth': 0.18, 'num_wavelengths': 8},
                     'method': {'starting_scale': 4}}

luvoir = DesignParameterSurvey(APLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, N, machine), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)

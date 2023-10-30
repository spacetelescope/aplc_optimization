'''This is the launcher template for the SCDA segment size survey.'''

import os

os.chdir('../..')
instrument = 'usort'  # do not edit
from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC

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
survey_name = 'onaxis_survey'
machine = 'telserv3'
N = 512

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
    An arbitrary padding of gap size to represent gaps on smaller arrays. This effectively makes the larger gaps larger 
    and the segments smaller to preserve the same segment pitch.
gap_size : float
    The gap size between the segments in meters.
num_rings : int
    The number of rings of hexagons to include, not counting the central segment.
segment_point_to_point : float, optional
    The point-to-point size of each of the hexagonal segments in meters.
obscuration_diameter : float, optional
    The diameter of the central obscuration in meters.
bottom_struts_width : float, optional
    The width of the bottom struts in meters. 
bottom_struts_connection_width : float, optional
    The width of the connection of the bottom struts to the central obscuration
    in meters. 
bottom_struts_connection_length : float, optional
    The length of the connection of the bottom struts to the origin, in meters.
bottom_struts_angle : float, optional
    The angle between the bottom struts in degrees. 
top_strut_width : float, optional
    The width of the bottom struts in meters. 
top_strut_connection_width : float, optional
    The width of the connection of the bottom struts to the central obscuration
    in meters. 
top_strut_connection_length : float, optional
    The length of the connection of the bottom struts to the origin, in meters.
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
pupil_diameter = 7.87
pupil_inscribed = 6.513
oversamp = 16
num_rings = 2
gap_size = 8e-3
segment_point_to_point = 1800e-3
obscuration_diameter = 787.1e-3
bottom_struts_width = 100.7e-3
bottom_struts_connection_width = 151.9e-3
bottom_struts_connection_length = 589.7e-3
bottom_struts_angle = 60
top_strut_width = 101.6e-3
top_strut_connection_width = 194.8e-3
top_strut_connection_length = 505.3e-3
offaxis = False
seg_gap_pad = 1

# Lyot stop parameters
circular = True
LS_ID = 0
LS_OD = 0.99

# INPUT FILE GENERATION
# gather the pre-generated input mask files (generated via
# notebooks/USORT/USORTApertureGenerationPython.ipynb)
pup_filename = 'USORT/TelAp_USORT_onaxis_ovsamp16_N0512.fits'
circLS_filename = f'USORT/LS_USORT_circ_ID0000_OD0990_ovsamp16_N0512.fits'
hexLS_filename = f'USORT/LS_USORT_hex_ID0000_OD0990_ovsamp16_N0512.fits'
ls_filename = circLS_filename, hexLS_filename

'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple `survey_parameters` dictionaries
  (as shown in the commented block, at bottom of this script).

Parameters
----------  
alignment_tolerance: int
    The Lyot stop alignment tolerance in pixels. 
num_lyot_stops: int
    The number of translated Lyot stops to optimize for simultaneously. If 1, design will be non-robust.
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
alignment_tolerance = 1 #0.2%
num_lyot_stops = 9

# FPM Parameters
radius = 3.4
num_pix = 150
grayscale = True
field_stop_radius = None, 14.5

# Optimization parameters (dark zone constraints)
iwa = radius - 0.1
owa = 14.0
bandwidth = 0.1, 0.15, 0.2
num_wavelengths = 5
contrast = 10

# Loop over FPM radius/IWA combinations, i.e. NOT as a grid
for i in range(5):
    radius += 0.1
    iwa = radius - 0.1

    # SURVEY PARAMS DICTIONARY
    survey_parameters = {'instrument': {'inst_name': instrument.upper()},
                        'pupil': {'N': N, 'filename': pup_filename},
                        'lyot_stop': {'filename': ls_filename, 'alignment_tolerance': alignment_tolerance,
                                      'num_lyot_stops': num_lyot_stops},
                        'focal_plane_mask': {'radius': radius, 'num_pix': num_pix, 'grayscale': grayscale, 'field_stop_radius': field_stop_radius},
                        'image': {'contrast': contrast, 'iwa': iwa, 'owa': owa, 'bandwidth': bandwidth,
                                   'num_wavelengths': num_wavelengths}}

    # RUN DESIGN SURVEY
    SCDA = DesignParameterSurvey(APLC, survey_parameters,
                                   'surveys/{}_{}_N{:04d}_{}/FPM{}_IWA{}'.format(instrument, survey_name, N, machine, int(10*radius), int(10*iwa)), 'masks/')
    SCDA.describe()
    SCDA.write_drivers(True)
    SCDA.run_optimizations(True)
    SCDA.run_analyses(True)


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

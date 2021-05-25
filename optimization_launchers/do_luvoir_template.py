from aplc_optimization.survey import DesignParameterSurvey
from aplc_optimization.aplc import APLC
from aplc_optimization.Inputs_Generation.LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

'''
DO NOT RUN THIS SCRIPT - this is the launcher template for LUVOIR. 
Default parameters for 1000 pix BW10 small design

Workflow: 
- Make a copy
- rename the copy under the following naming schema: 'do_luvoir_<survey_name>_<machine>.py', designating the survey you 
    are running and the name of the machine you will be running it on. 
    (e.g. 'do_luvoir_BW10_small_telserv3.py' designates a BW small design run on telserv3.)
- Edit input file parameters and survey parameters, as desired. 
- Once the survey is complete, move the launcher script to the 'survey/' directory..
'''

# Survey information
instrument = 'luvoir'                                   # instrument name;
survey_name = 'template'             # survey name;
machine = 'local'                                     # name of machine the survey is run on.

# Physical
pupil_diameter = 15.0 #m: actual LUVOIR A circumscribed pupil diameter
pupil_inscribed = 13.5 #m: actual LUVOIR A inscribed pupil diameter

'''
Input (aperture and Lyot stop) Array Parameters
-----------------------------------------------
Gap padding (seg_gap_pad) and grey levels (oversamp) are set according to number of input array pixels (nArray),  
configured in order to keep gap size as close to actual physical size of LUVOIR A, as possible. 

 - nArray = 1000, oversamp = 4, seg_gap_pad = 1
 - nArray = 500,  oversamp = 3, seg_gap_pad = 2
 - nArray = 300,  oversamp = 4, seg_gap_pad = 4
 - nArray = 200,  oversamp = 4, seg_gap_pad = 4
 - nArray = 100,  oversamp = 4, seg_gap_pad = 4
'''

# Input parameters
filepath = instrument.upper()+'/'   # directory in 'masks/' where the input files are stored
nArray = 100    # number of pixels in input (TelAP, LS) and final (apodizer) arrays
oversamp = 4    # oversampling factor (number of grey levels) - if set to 1 will return a bw pupil, for grey set to > 1

# Aperture parameters
gap_padding = 4 # arbitary padding of gap size to represent gaps on smaller arrays
                # effectively makes the larger gaps larger and the segments smaller to preserve the same segment pitch

# Lyot stop parameters
lyot_ref_diameter = pupil_inscribed # diameter used to reference the LS inner and outer diameter against
ls_spid = False# whether to include secondary support mirror structure in the aperture
spid_ov = 2     # factor by which to oversize the spiders compared to the LUVOIR-A aperture spiders
ls_id = 0.12    # LS inner diameter(s) as a fraction of `lyot_ref_diameter`
ls_od = 0.937   # LS outer diameter as a fraction of `lyot_ref_diameter`
        # Note: both ls_id and ls_od are re-normalized against circumscribed pupil diameter during LS generation

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {'directory': 'LUVOIR/', 'N': nArray, 'oversamp': oversamp,
                    'aperture': {'seg_gap_pad': gap_padding},
                    'lyot_stop':{'ls_spid':ls_spid, 'ls_spid_ov': spid_ov, 'lyot_ref_diam': lyot_ref_diameter, 'LS_ID':[ls_id], 'LS_OD':[ls_od]}}

# INPUT FILE GENERATION
pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)


'''
Survey Design Parameters
------------------------
- for multiple design parameters as a grid, input as list
- for multiple design parameters NOT as a grid, create multiple entries of below 
  (as shown in the commented block, at bottom of this script)

'''
# FPM Parameters
FPM_radius = 3.5 # lamda_0/D radius
nFPM = 150 # number of pixels in the focal plane mask
greyscale = True

# Optimization parameters (dark zone constraints)
IWA = 3.4   #lam/D: effective inner working angle (outer perimeter of the annular dark zone in coronagraphic image)
OWA = 12.0  #lam/D: effective outer working angle (inner perimeter of the annular dark zone in the coronagraphic image)
bandpass = 0.1  # dark zone bandpass (fractional)
nLams = 1   # number of wavelengths spanning bandpass
contrast = 10   # contrast goal in the dark zone

# SURVEY PARAMS DICTIONARY
survey_parameters = {'pupil': {'N': nArray, 'filename': pup_filename},
                     'lyot_stop': {'filename': ls_filenames},
                     'focal_plane_mask': {'radius':FPM_radius, 'num_pix': nFPM, 'grayscale': greyscale},
                     'image': {'contrast': contrast,'iwa': IWA,'owa': OWA,'bandwidth': bandpass,'num_wavelengths':nLams}}

# RUN DESIGN SURVEY
luvoir = DesignParameterSurvey(APLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument, survey_name, nArray, machine),
                               'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)


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
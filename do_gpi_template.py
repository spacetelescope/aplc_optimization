from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
#from GPI_Inputs_Generation import GPI_inputs_gen

'''
DO NOT RUN THIS SCRIPT - It is the launcher template. 
Default parameters for 1000 pix BW10 small design

Workflow: 
- Make a copy
- rename the copy to designate survey you are running and the machine name you will be 
  running it on (do_gpi_BW10_small_telserv3.py , for example, designates the BW10 small 
  design run on telserv3)

- set survey_name the same as survey name you use in the launcher file name
- Set the machine to be the name of the machine you're running on (need to write this to 
  grab machine name automatically)

- run the launcher
- after survey is done, move the launcher script into the survey directory
'''

nArray = 512  # number of pixels in input and final arrays
instrument = "gpi"  # instrument name
survey_name = "test"  # survey name
machine = "telserv3"  # machine the survey is run on.

#oversamp = 2 # oversampling factor for intermediate plane calculations.

ap_spid = False # the secondary supports are masked out in the Lyot planes regardless
satspots = False # whether to include satellite spots grid
ls_tabs = False # whether to include the bad-actuator masking tabs
lyot_mask = '080m12_03' # name of Lyot mask. Available lyot masks are:
                        #  '080m12_02', '080m12_03', '080m12_04', '080m12_04_c', '080m12_06',
                        #  '080m12_06_03', '080m12_07', '080m12_10', 'Open' and 'Blank'.

# Lyot mask geometry as projected (rescaled) onto the primary mirror.
ls_id = 7.57	# m: Lyot mask inner diameter
ls_od = 2.1954  # m: Lyot mask outer diameter
ls_spid_width = 0.24042 # m: support width


FPM = 156.2e-3 #arcsec: FPM diameter. {Y: 156.2, J: 184.7, H: 246.7, K1: 306.3, SCIENCE: 0.0}

input_files_dict = {'directory': 'GPI/', 'N': nArray,
                    'aperture': {'ap_spid': False},
                    'lyot_stop': {'lyot_mask':lyot_mask, 'ls_tabs': ls_tabs}}

#pup_filename, ls_filename = GPI_inputs_gen(input_files_dict)
pup_filename = 'GPI/GPI_primary_N0512.fits'
ls_filename = 'GPI/GPI_LS_N0512_080m12_03_notabs.fits'


ls_tabs = True # Whether to block out the bad actuators with tabs in the Lyot masks
satspots = True # Whether to include the satellite spots in the apodizer model

# Focal plane mask parameters
nFPM = nArray # number of pixels across the focal plane mask
FPM_name = 'Y' # Name of FPM. Available FPMs are: 'Y', 'J', 'H' or 'K1'
FPM_radius = 5.6/2  #lambda_0/D: focal plane mask radius

# Optimization parameters
contrast = 8 # 10-<value>: contrast goal in the dark zone of the coronagraphic image
IWA = 3.75 #lamda_0/D: inner edge of dark zone region in coronagraphic image
OWA = 15 #lambda_0/D: outer edge of dark zone region in coronagraphic image
spectral_bandwidth = 0.0 # fractional, dark zone bandpass
nLams = 1 # number of wavelengths spanning the design bandpass

# method parameters
starting_scale = 1 # number of pixels per unit cell for the initial solution.

survey_parameters = {'pupil': {'N': nArray, 'filename': pup_filename, 'satspots': True}, \
                     'lyot_stop': {'filename': ls_filename, 'ls_tabs': True}, \
                     'focal_plane_mask': {'FPM_name': FPM_name, 'radius': FPM_radius, 'num_pix': nFPM},
                     'image': {'contrast': contrast,'iwa': IWA,'owa': OWA,'bandwidth': spectral_bandwidth,'num_wavelengths':nLams}}




# RUN DESIGN SURVEY
gpi = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument,survey_name,nArray,machine), 'masks/')
gpi.describe()

gpi.write_drivers(True)
gpi.run_optimizations(True)
gpi.run_analyses(True)
from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

'''
DO NOT RUN THIS SCRIPT - It is the launcher template. 
Default parameters for 1000 pix BW10 small design

Workflow: 
- Make a copy
- rename the copy to designate survey you are running and the machine name you will be 
  running it on (do_luvoir_BW10_small_telserv3.py , for example, designates the BW10 small 
  design run on telserv3)
- run the launcher
- after survey is done, move the launcher script into the survey directory
'''
n = 1000 #number of pixels in input and final arrays

filename_info = os.path.basename(__file__)[:-3].split("_", 4) #extract info from the launcher file name
instrument = filename_info[1] #instrument name
survey_name = filename_info[2]+'_'+filename_info[3] #survey name (automatically set)
machine = filename_info[4] #name of machine the survey is run on (automatically set)

'''
input array (telap and ls) parameters
oversamp (grey levels) and seg_gap_pad set according to n, minimizing size of gaps to be 
as close to actual size i.e. minimizing the amount of padding needed

n = 1000, oversamp = 4, seg_gap_pad = 1
n = 500,  oversamp = 3, seg_gap_pad = 2
n = 300,  oversamp = 4, seg_gap_pad = 4
n = 200,  oversamp = 4, seg_gap_pad = 4
n = 100,  oversamp = 4, seg_gap_pad = 4
'''

input_files_dic = {'directory':'LUVOIR/', 'N':n, 'oversamp':4, \
                    'aperture': {'seg_gap_pad':1}, \
                    'lyot_stop':{'lyot_ref_diam':13.5,'ls_spid':False,'ls_spid_ov':2,'LS_ID':[0.12], 'LS_OD':[0.982]}}

pup_filename, ls_filenames = LUVOIR_inputs_gen(input_files_dict)


#design parameters
#for multiple design parameters as a grid, input as list
#for multiple design parameters NOT as a grid, make multiple entries of below (as shown in commented out block)
#resulting apodizer are packaged with their inputs in a fits cube

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.5, 'num_pix': 150, 'grayscale': True,},
                     'image': {'contrast':10,'iwa':3.4,'owa':12.0,'bandwidth':0.10,'num_wavelengths':5}, \
                     'method':{'starting_scale': 1}}


luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument,survey_name,n,machine), 'masks/')
luvoir.describe()

luvoir.write_drivers(True)
luvoir.run_optimizations(True)
luvoir.run_analyses(True)


'''
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

					
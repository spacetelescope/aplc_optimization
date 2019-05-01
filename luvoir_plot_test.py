from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

n = 1000
pup_filename = 'LUVOIR/TelAp_LUVOIR_gap_pad01_bw_ovsamp03_N1000.fits'
ls_filenames = 'LUVOIR/LS_LUVOIR_ID0190_OD0937_struts_pad02_gy_ovsamp3_N1000.fits'

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':3.82},
                     'image': {'contrast':10,'iwa':3,'owa':12,'bandwidth':0.18,'num_wavelengths':8}}

luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'survey/luvoir/', 'masks/')
luvoir.describe()


luvoir.run_analyses(overwrite=True)
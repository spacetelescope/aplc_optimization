from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from LUVOIR_Inputs_Generation import LUVOIR_inputs_gen

n = 200
pup_filename = 'LUVOIR/TelAp_LUVOIR_gap_pad01_bw_ovsamp03_N1000.fits'
ls_filenames = 'LUVOIR/LS_LUVOIR_ID0190_OD0937_no_struts_gy_ovsamp4_N0200.fits'

survey_parameters = {'pupil': {'N': n,'filename': pup_filename}, \
                     'lyot_stop': {'filename': ls_filenames}, \
                     'focal_plane_mask': {'radius':14.17},
                     'image': {'contrast':10,'iwa':15.05,'owa':60.20,'bandwidth':0.18,'num_wavelengths':8}}

luvoir = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/luvoir_N{:04d}_large_telserv3/'.format(n), 'masks/')
luvoir.describe()


luvoir.run_analyses(overwrite=True)
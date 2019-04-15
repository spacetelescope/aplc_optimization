from survey import DesignParameterSurvey
from por_aplc import PorAPLC

varying_parameters = {'pupil': {'filename': 'ehpor_hicat_apodizer_mask_512_gy.fits'}, 'lyot_stop': {'filename': 'ehpor_hicat_lyot_mask_512_gy_2.fits'}, 'image': {'owa': 12}}

survey = DesignParameterSurvey(PorAPLC, varying_parameters, 'magnfication_survey1/', 'masks/')
survey.describe()

survey.write_drivers(False)
survey.run_optimizations(False)
survey.run_analyses(True)

varying_parameters = {'pupil': {'filename': 'ehpor_hicat_apodizer_mask_512_gy.fits'}, 'lyot_stop': {'filename': 'ehpor_hicat_lyot_mask_512_gy_{:d}.fits', 'num_stops': 5}, 'image': {'owa': 12}}

survey = DesignParameterSurvey(PorAPLC, varying_parameters, 'magnfication_survey5/', 'masks/')
survey.describe()

survey.write_drivers(False)
survey.run_optimizations(False)
survey.run_analyses(True)

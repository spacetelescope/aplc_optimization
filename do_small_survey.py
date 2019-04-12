from survey import DesignParameterSurvey
from por_aplc import PorAPLC


varying_parameters = {'pupil': {'filename': 'ehpor_apodizer_mask_256_bw.fits'}, 'lyot_stop': {'filename': 'ehpor_lyot_mask_256_gy.fits'}, 'image': {'owa': 12}}



survey = DesignParameterSurvey(PorAPLC, varying_parameters, 'survey/', 'masks/')
survey.describe()

survey.write_drivers(False)
survey.run_optimizations(False)
survey.run_analyses(True)

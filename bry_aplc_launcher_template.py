from survey import DesignParameterSurvey
from por_aplc import PorAPLC
from astropy.io import fits
from HiCAT_Inputs_generation import HiCAT_inputs_gen

'''
    DO NOT EDIT: this is a launcher template to be used with the `creating_a_design_survey.ipynb` notebook.
'''

# SURVEY INFORMATION
instrument = instrument_name
survey_name = survey_name
machine = machine
N = N

# INPUT FILES PARAMETER DICTIONARY
input_files_dict = {}

# INPUT FILE GENERATION
pup_filename, ls_filenames =\
    HiCAT_inputs_gen(input_files_dict)

# SURVEY PARAMETER DICTIONARY
survey_parameters = {}

# RUN DESIGN PARAMETER SURVEY
Survey = DesignParameterSurvey(PorAPLC, survey_parameters, 'surveys/{}_{}_N{:04d}_{}/'.format(instrument,survey_name,N,machine), 'masks/')
Survey.describe()

Survey.write_drivers(True)
Survey.run_optimizations(True)
Survey.run_analyses(True)
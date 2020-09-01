from datetime import datetime

from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

def launch_script_gen(instrument, survey_name, machine, input_file_dictionary, survey_parameters):
    """Generates a launcher script for the parameters provided in 'Generate_a_launcher_script.ipynb'.

        Parameters
        ----------
        instrument: str
            name of the instrument
        survey: str
            designated name of the design survey
        machine: str
            name of the machine that the survey will be run on.
        input_files_dict: dict
            A dictionary of input parameters (defined in the launcher script).
        survey_parameters
        """

    launcher_fname = 'do_'+instrument+'_'+survey_name+'_'+machine+'.py'

    config = Path(launcher_fname)
    if config.is_file():
        print('A launcher script with the name "{0:s}" already exists.\n'.format(launcher_fname))
        while True:
            overwrite = input('Would you like to "{:s}"? Y = yes, N = No \n'.format(launcher_fname))
            if overwrite == 'Y':
                break
            elif overwrite == 'N':
                print('No new launcher scripts have been generated.')
                sys.exit()
            else:
                continue

    with open('bry_aplc_launcher_template.py') as template_file:
        template_content = template_file.readlines()
        template_file.close()

    new_content = template_content.copy()

    new_content[3] = 'from {0:s}_Inputs_generation import {0:s}_inputs_gen\n'.format(instrument)

    timestamp = datetime.today().strftime('%Y-%m-%d')
    new_content[6] = "    Launcher script created {:s} using `Creating_a_design_survey.ipynb` notebook.\n".format(
        timestamp)

    new_content[10] = "instrument = '" + instrument + "'\n"
    new_content[11] = "survey_name = '" + survey_name + "' \n"
    new_content[12] = "machine = '" + machine + "' \n"

    nArray = input_file_dictionary['N']
    new_content[13] = "N = '" + str(nArray) + "' \n"

    new_content[16] = "input_files_dict = " + str(input_file_dictionary) + "\n"
    new_content[20] = "    " + instrument + "_inputs_gen(input_files_dict)\n"

    new_content[23] = "survey_parameters = " + str(survey_parameters) + "\n"

    launcher_file = open(launcher_fname, 'w')
    launcher_file_content = "".join(new_content)
    launcher_file.write(launcher_file_content)
    launcher_file.close()

    print('\n New launcher script saved as "{:s}".\n'.format(launcher_fname))
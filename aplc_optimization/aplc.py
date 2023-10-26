import os
import shutil
import sys

from aplc_optimization import analysis
from aplc_optimization.survey import Coronagraph


class APLC(Coronagraph):
    """The class for the Apodized Pupil Lyot Coronagraph design.

    Attributes
    ----------
    identifier: str
        The APLC design identifier.
    parameters: dict
        The design parameters of the APLC design.
    file_organization: dict
        The file organization structure of the APLC design.
    """
    _default_parameters = {
        'instrument': {
            'inst_name': 'LUVOIR'
        },
        'pupil': {
            'filename': 'NoFilename',
            'N': None
        },
        'focal_plane_mask': {
            'radius': 4.0,
            'num_pix': 50,
            'grayscale': True,
            'field_stop_radius': None
        },
        'lyot_stop': {
            'filename': 'NoFilename',
            'alignment_tolerance': 0,
            'num_lyot_stops': 1
        },
        'image': {
            'contrast': 8.0,
            'iwa': 3.75,
            'owa': 15.0,
            'num_wavelengths': 3,
            'bandwidth': 0.1,
            'resolution': 2
        },
        'method': {
            'force_no_x_mirror_symmetry': False,
            'force_no_y_mirror_symmetry': False,
            'force_no_hermitian_symmetry': False,
            'starting_scale': 1,
            'ending_scale': 1,
            'edge_width_for_prior': 2,
            'num_throughput_iterations': 2,
            'initial_throughput_estimate': 1,
            'maximize_planet_throughput': True
        },
        'solver': {
            'num_threads': 0,
            'crossover': 0,
            'method': 2
        }
    }

    def __init__(self, identifier, parameters, file_organization):
        super(APLC, self).__init__(identifier, parameters, file_organization, analysis)

    def write_driver(self, overwrite=False):
        """Write a driver script for the design using the template ('driver_template.py') script and write to the
        survey's driver directory.

        Driver script is automatically written with "parameters" and "file_organization".

        Parameters
        ----------
        overwrite: bool
            Whether to overwrite the driver file if it already exists.
        """
        if self.check_driver() and not overwrite:
            print('Driver already exists and will not be overwritten.')
            return

        driver = 'parameters = {:s}\nfile_organization = {:s}\nsolution_fname = "{:s}"\n\n'.format(str(self.parameters),
                                                                                                   str(
                                                                                                       self.file_organization),
                                                                                                   self.solution_filename)

        with open('aplc_optimization/driver_template.py') as template_file:
            driver_template = template_file.read()
        driver += driver_template

        fname = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')
        with open(fname, 'w') as output_file:
            output_file.write(driver)

        fname_optimizer = os.path.join(self.file_organization['drivers_dir'], 'optimizer.py')
        if os.path.exists(fname_optimizer) and not overwrite:
            print('Optimizer already exists and will not be overwritten.')
        else:
            shutil.copy('aplc_optimization/optimizer.py', fname_optimizer)

    def check_input_files(self):
        """Check whether all of the necessary input files for the design exist."""
        pup_fname = self.parameters['pupil']['filename']
        ls_fname = self.parameters['lyot_stop']['filename']

        if not os.path.isabs(pup_fname):
            pup_fname = os.path.join(self.file_organization['input_files_dir'], pup_fname)
        if not os.path.isabs(ls_fname):
            ls_fname = os.path.join(self.file_organization['input_files_dir'], ls_fname)

        return os.path.exists(pup_fname) and os.path.exists(ls_fname)

    def check_driver(self):
        """Check whether the driver script for the design already exists."""
        fname_driver = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')
        fname_optimizer = os.path.join(self.file_organization['drivers_dir'], 'optimizer.py')

        return os.path.exists(fname_driver) and os.path.exists(fname_optimizer)

    def get_driver_command(self):
        """Get the bash command to run the driver script."""
        fname_driver = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')

        return '{:s} {:s} &> {:s}'.format(sys.executable, fname_driver, self.log_filename)

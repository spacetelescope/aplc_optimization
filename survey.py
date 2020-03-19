import collections
import copy
import csv
import datetime
import getpass
import inspect
import itertools
import os
import pprint
import socket
import warnings

import asdf
import matplotlib as mpl
import numpy as np
import six

mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages


# Check if the argument is iterable but not string_like
def is_iterable(arg):
    return isinstance(arg, collections.Iterable) and not isinstance(arg, six.string_types)


def mark_slow(func):
    if hasattr(func, 'marks'):
        func.marks.append('slow')
    else:
        func.marks = ['slow']
    return func


def is_marked(function, mark):
    if hasattr(function, 'marks'):
        return mark in function.marks
    else:
        return False


'''
File organization structure:
- survey_dir: The base directory for which all files related to this survey are written.
- solution_dir: The subdirectory where all solutions will be written.
- analysis_dir: The subdirectory where all analysis files will be written.
- drivers_dir: The subdirectory where all Python scripts will be written.
- log_dir: The subdirectory where all log files will be written.
- input_files_dir: The directory where all input files for the optimization are located.
'''


class DesignParameterSurvey(object):
    def __init__(self, coronagraph_class, parameter_sets, survey_dir, input_files_dir):
        self.coronagraph_class = coronagraph_class

        # Set up directories
        survey_dir = os.path.abspath(survey_dir)
        self.file_organization = {'survey_dir': survey_dir}
        self.file_organization['solution_dir']    = os.path.join(survey_dir, 'solutions')
        self.file_organization['analysis_dir']    = os.path.join(survey_dir, 'analysis')
        self.file_organization['drivers_dir']  	  = os.path.join(survey_dir, 'drivers')
        self.file_organization['log_dir']         = os.path.join(survey_dir, 'logs')
        self.file_organization['input_files_dir'] = os.path.abspath(input_files_dir)

        # Make sure all directories exist
        for key in self.file_organization:
            if not os.path.exists(self.file_organization[key]):
                os.makedirs(self.file_organization[key])

        self.varied_parameters = []
        self.varied_parameters_indices = []

        self.default_parameters = coronagraph_class._default_parameters.copy()

        # Find all fixed and varied parameters
        for category in parameter_sets:
            if category not in coronagraph_class._default_parameters:
                warnings.warn(
                    'Unrecognized parameter category "{0}". All parameters in this category will be ignored.'.format(
                        category))
                continue

            for key in parameter_sets[category]:
                if key not in coronagraph_class._default_parameters[category]:
                    warnings.warn(
                        'Unrecognized parameter name "{0}" in category "{1}". This parameter will be ignored.'.format(
                            key, category))
                    continue

                if is_iterable(parameter_sets[category][key]):
                    if len(parameter_sets[category][key]) == 1:
                        # It is a fixed parameter
                        self.default_parameters[category][key] = parameter_sets[category][key][0]
                    else:
                        # It is a varied parameter
                        self.varied_parameters.append(parameter_sets[category][key])
                        self.varied_parameters_indices.append((category, key))
                else:
                    # It is a fixed parameter
                    self.default_parameters[category][key] = parameter_sets[category][key]

        # Make list of fixed parameters
        self.fixed_parameters_indices = {}
        for category in self.default_parameters:
            for key in self.default_parameters[category]:
                if (category, key) not in self.varied_parameters_indices:
                    if category in self.fixed_parameters_indices:
                        self.fixed_parameters_indices[category].append(key)
                    else:
                        self.fixed_parameters_indices[category] = [key]

        # Create coronagraph objects with correct parameters
        self.coronagraphs = []
        self.parameter_sets = []
        num_parameter_sets = 1
        for p in self.varied_parameters:
            num_parameter_sets *= len(p)

        print(self.varied_parameters)

        # format_string = '{:0' + str(len(str(num_parameter_sets))) +  'd}'

        format_string = '{:0' + str(len(str(
            num_parameter_sets))) + 'd}' + '_LUVOIR_N{:}_FPM{:3d}M0{:d}_IWA{:04d}_OWA0{:04d}_C{:d}_BW{:d}_Nlam{:d}_LS_ID{:s}_OD{:s}_{:s}'

        params = list(itertools.product(*self.varied_parameters))
        if len(self.varied_parameters) == 0:
            params = [{}]

        for i, combo in enumerate(params):
            # Create parameter set for this coronagraph
            new_parameter_set = copy.deepcopy(self.default_parameters)
            for value, (category, key) in zip(combo, self.varied_parameters_indices):
                print(category, key, value)
                new_parameter_set[category][key] = value

            # LUVOIR/LS_LUVOIR_ID0190_OD0937_no_struts_gy_ovsamp2_N0050.fits

            N = new_parameter_set['pupil']['N']
            fpm = int(100 * new_parameter_set['focal_plane_mask']['radius'])
            m = new_parameter_set['focal_plane_mask']['num_pix']
            iwa = int(100 * new_parameter_set['image']['iwa'])
            owa = int(100 * new_parameter_set['image']['owa'])
            c = int(new_parameter_set['image']['contrast'])
            bw = int(100 * new_parameter_set['image']['bandwidth'])
            nlam = new_parameter_set['image']['num_wavelengths']
            ls_id = new_parameter_set['lyot_stop']['filename'][18:23]
            ls_od = new_parameter_set['lyot_stop']['filename'][26:30]

            ls_strut_check = new_parameter_set['lyot_stop']['filename'][31]

            if ls_strut_check == 'n':
                ls_strut_key = 'no_ls_struts'
            else:
                ls_strut_key = 'ls_' + new_parameter_set['lyot_stop']['filename'][31:43]

            # new_parameter_set['']['']

            # Create unique id
            # identifier = format_string.format(i)
            identifier = format_string.format(i, N, fpm, m, iwa, owa, c, bw, nlam, ls_id, ls_od, ls_strut_key)

            # Create coronagraph
            self.coronagraphs.append(coronagraph_class(identifier, new_parameter_set, self.file_organization))
            self.parameter_sets.append(new_parameter_set)

    def union(self, b):
        pass

    def describe(self):
        print('This survey has {:d} design parameter combinations.'.format(len(self.parameter_sets)))
        print('{:d} parameter are varied:'.format(len(self.varied_parameters)))
        for values, (category, key) in zip(self.varied_parameters, self.varied_parameters_indices):
            print('   {:s} - {:s}: {:s}'.format(category, key, str(values)))
        print('')
        print('File organization:')
        pprint.pprint(self.file_organization)
        print('')
        print('All input files exist? {}'.format(self.check_input_files()))
        print('All drivers exist? {}'.format(self.check_drivers()))
        print('All solutions exist? {}'.format(self.check_solutions()))

    def check_input_files(self):
        num_incomplete = 0

        for cor in self.coronagraphs:
            if not cor.check_input_files():
                num_incomplete += 1

        return num_incomplete == 0

    def check_drivers(self):
        num_incomplete = 0

        for cor in self.coronagraphs:
            if not cor.check_driver():
                num_incomplete += 1

        return num_incomplete == 0

    def check_solutions(self):
        num_incomplete = 0

        for cor in self.coronagraphs:
            if not cor.check_solution():
                num_incomplete += 1

        return num_incomplete == 0

    def write_drivers(self, overwrite=False):
        for cor in self.coronagraphs:
            cor.write_driver(overwrite)

    def write_serial_bash_script(self, overwrite=False):
        fname = os.path.join(self.file_organization['drivers_dir'], 'run.sh')

        if os.path.exists(fname) and not overwrite:
            print('Serial bash script already exists and is not overwritten.')
            return

        with open(fname, 'w') as f:
            for coronagraph in self.coronagraphs:
                f.write(coronagraph.get_driver_command() + '\n')

    def run_optimizations(self, force_rerun=False):
        for cor in self.coronagraphs:
            cor.run_optimization(force_rerun)

        return self.check_solutions()

    def write_spreadsheet(self, overwrite=False):
        fname_tail = "{1:s}_{2:s}.csv".format(getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
        fname = os.path.join(self.file_organization['survey_dir'], fname_tail)

        with open(fname, 'wb') as survey_spreadsheet:
            survey = csv.writer(survey_spreadsheet)

            # Write header
            survey.writerow(["Created by {:s} on {:s} at {:s}".format(getpass.getuser(), socket.gethostname(),
                                                                      datetime.datetime.now().strftime(
                                                                          "%Y-%m-%d %H:%M"))])
            survey.writerow([])
            survey.writerow(['File organization:'])
            survey.writerow(['Survey Directory', self.file_organization['survey_dir']])
            survey.writerow(['Solution Directory', self.file_organization['solution_dir']])
            survey.writerow(['Analysis Directory', self.file_organization['analysis_dir']])
            survey.writerow(['Drivers Directory', self.file_organization['drivers_dir']])
            survey.writerow(['Log Directory', self.file_organization['log_dir']])
            survey.writerow(['Input Files Directory', self.file_organization['iput_files_dir']])
            survey.writerow([])

            # Write fixed parameters
            survey.writerow(['Fixed design parameters:'])
            for category in self.fixed_parameters_indices:
                survey.writerow([category.upper()])
                for key in self.fixed_parameters_indices[category]:
                    survey.writerow([key.lower(), str(self.default_parameters[category][key])])
            survey.writerow([])

            # Write varied parameters
            survey.writerow(
                ['Varied design parameters ({:d} total combinations):'.format(len(self.varied_parameters_indices))])
            for category, key in self.varied_parameters_indices:
                survey.writerow([category.upper(), key.lower(), str(self.varied_parameters[category][key])])
            survey.writerow([])

            # Get keys for metrics to write out
            keys = set()
            for coronagraph in self.coronagraphs:
                keys.update(coronagraph.metrics.keys())

            # The keys to be written to the spreadsheet cannot be arrays.
            for coronagraph in self.coronagraphs:
                for key in keys:
                    if key in coronagraph.metrics:
                        if not np.isscalar(coronagraph.metrics[key]):
                            keys.discard(key)
            keys = sorted(list(keys))

            # Write individual optimizations
            survey.write('Individual optimizations:')
            survey.writerow([category.upper() for category, key in self.varied_parameters])
            header = [key.lower() for category, key in self.varied_parameters]
            header += ['input files?', 'Driver?', 'Solution?']
            header += keys
            survey.writerow(header)

            for coronagraph in self.coronagraphs:
                coro_row = [self.parameter_sets[category][key] for (category, key) in self.varied_parameters]
                coro_row.append('Y' if coronagraph.check_input_files() else 'N')
                coro_row.append('Y' if coronagraph.check_driver() else 'N')
                coro_row.append('Y' if coronagraph.check_solution() else 'N')
                for key in keys:
                    if key in coronagraph.metrics:
                        coro_row.append(str(coronagraph.metrics[key]))
                    else:
                        coro_row.append('')
                survey.writerow(coro_row)

        os.chmod(fname, 644)

    def run_analyses(self, overwrite=False, run_slow=True):
        for coronagraph in self.coronagraphs:
            coronagraph.run_analysis(overwrite, run_slow)


class Coronagraph(object):
    def __init__(self, identifier, parameters, file_organization, analysis_module=None):
        self._identifier = identifier
        self.parameters = parameters
        self.file_organization = file_organization
        self.analysis_module = analysis_module
        self.metrics = {}

    @property
    def identifier(self):
        return self._identifier

    def check_input_files(self):
        raise NotImplementedError()

    def check_driver(self):
        raise NotImplementedError()

    def check_solution(self):
        return os.path.exists(self.solution_filename)

    @property
    def solution_filename(self):
        return os.path.join(self.file_organization['solution_dir'], self.identifier + '.fits')

    @property
    def log_filename(self):
        return os.path.join(self.file_organization['log_dir'], self.identifier + '.log')

    def get_driver_command(self):
        raise NotImplementedError()

    def run_optimization(self, force_rerun=False):
        if (not self.check_driver()) or (not self.check_input_files()):
            print('Not all driver or input files are written.')
            return

        if self.check_solution() and not force_rerun:
            print('Solution already exists.')
            return

        os.system(self.get_driver_command())

    def run_analysis(self, overwrite=False, run_slow=True):
        if self.analysis_module is None:
            print('No analysis module was provided. No analysis will be performed.')
            return

        if not os.path.exists(self.solution_filename):
            print('The solution is not optimized yet. No analysis will be performed.')
            return

        print(self.solution_filename)

        # Read in functions from analysis_module and sort by name
        analysis = inspect.getmembers(self.analysis_module, inspect.isfunction)
        analysis = [x for _, x in sorted(zip([name for name, function in analysis], analysis))]

        analysis_summary_filename = os.path.join(self.file_organization['analysis_dir'], self.identifier + '.pdf')
        analysis_metrics_filename = os.path.join(self.file_organization['analysis_dir'], self.identifier + '.asdf')

        if os.path.exists(analysis_summary_filename) and os.path.exists(analysis_metrics_filename) and not overwrite:
            print('Analysis already exists and is not overwritten.')
            return

        self.metrics = {}

        with PdfPages(analysis_summary_filename) as pdf:
            # Run all analysis functions
            for name, function in analysis:
                if not name.startswith('analyze_'):
                    continue

                if is_marked(function, 'slow') and not run_slow:
                    continue

                res = function(self.solution_filename, pdf)
                try:
                    self.metrics.update(res)
                except:
                    pass

        from asdf import AsdfFile
        # Write out metrics to a file
        f = AsdfFile(self.metrics)
        f.write_to(analysis_metrics_filename)

import collections
import six
import logging
import copy
import getpass
import asdf
import warnings
import itertools
import pprint
import os
import csv
import datetime
import socket
import asdf
import inspect

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

import por_aplc_analysis

# Check if the argument is iterable but not string_like
def is_iterable(arg):
	return isinstance(arg, collections.Iterable) and not isinstance(arg, six.string_types)

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
		self.file_organization['solution_dir'] = os.path.join(survey_dir, 'solutions')
		self.file_organization['analysis_dir'] = os.path.join(survey_dir, 'analysis')
		self.file_organization['drivers_dir'] = os.path.join(survey_dir, 'drivers')
		self.file_organization['log_dir'] = os.path.join(survey_dir, 'logs')
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
				warnings.warn('Unrecognized parameter category "{0}". All parameters in this category will be ignored.'.format(category))
				continue
			
			for key in parameter_sets[category]:
				if key not in coronagraph_class._default_parameters[category]:
					warnings.warn('Unrecognized parameter name "{0}" in category "{1}". This parameter will be ignored.'.format(key, category))
					continue
				
				if is_iterable(parameter_sets[category][key]):
					# It is a varied parameter
					self.varied_parameters.append(parameter_sets[category][key])
					self.varied_parameters_indices.append((category, key))
				else:
					# It is a fixed parameter
					self.default_parameters[category][key] = parameter_sets[category][key]
		
		# Make list of fixed parameters
		self.fixed_parameters_indices = {}
		for category in self.default_parameters:
			for key in self.default_parameters[key]:
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
		format_string = '{:.' + str(len(str(num_parameter_sets))) + 'd}'

		params = list(itertools.product(self.varied_parameters))
		if num_parameter_sets == 1:
			params = [{}]
		
		for i, combo in enumerate(params):
			# Create parameter set for this coronagraph
			new_parameter_set = self.default_parameters.copy()
			for value, (category, key) in zip(combo, self.varied_parameters_indices):
				new_parameter_set[category][key] = value
			
			# Create unique id
			identifier = format_string.format(i)
		
			# Create coronagraph
			self.coronagraphs.append(coronagraph_class(identifier, new_parameter_set, self.file_organization))
			self.parameter_sets.append(new_parameter_set)
	
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
		print('All drivers exist? {}')
		print('All output files exist? {}'.format(self.check_output_files()))
	
	def check_input_files(self):
		res = True

		for cor in self.coronagraphs:
			if not cor.check_input_files():
				res = False
		
		return res
	
	def check_output_files(self):
		res = True

		for cor in self.coronagraphs:
			if not cor.check_output_files():
				res = False
			
		return res
	
	def write_drivers(self, overwrite=False):
		for cor in self.coronagraphs:
			cor.write_driver(overwrite)
		self.coronagraph_class.write_static_files()
	
	def write_serial_bash_script(self, overwrite=False):
		fname = os.path.join(self.file_organization['drivers_dir'], 'run.sh')

		if os.path.exists(fname) and not overwrite:
			print('Serial bash script already exists and is not overwritten.')
			return

		with open(fname, 'w') as f:
			for coronagraph in self.coronagraphs:
				f.write(coronagraph.get_driver_command() + '\n')
	
	def run_serial(self, force_rerun=False):
		for cor in self.coronagraphs:
			if not cor.check_output_files():
				cor.run_optimization()
		
		return self.check_output_files()
	
	def write_spreadsheet(self, overwrite=False):
		fname_tail = "{1:s}_{2:s}.csv".format(getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
		fname = os.path.join(self.file_organization['survey_dir'], fname_tail)

		with open(fname, 'wb') as survey_spreadsheet:
			survey = csv.writer(survey_spreadsheet)

			# Write header
			survey.writerow(["Created by {:s} on {:s} at {:s}".format(getpass.getuser(), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))])
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
			survey.writerow(['Varied design parameters ({:d} total combinations):'.format(len(self.varied_parameters_indices))])
			for category, key in self.varied_parameters_indices:
				survey.writerow([category.upper(), key.lower(), str(self.varied_parameters[category][key])])
			survey.writerow([])

			# Write individual optimizations
			survey.write('Individual optimizations:')
			survey.writerow([category.upper() for category, key in self.varied_parameters])
			header = [key.lower() for category, key in self.varied_parameters] + ['input files?', 'Driver?', 'Solution?']
			survey.writerow(header)
			for coronagraph in self.coronagraphs:
				coro_row = [self.parameter_sets[category][key] for (category, key) in self.varied_parameters]
				coro_row.append('Y' if coronagraph.check_input_files() else 'N')
				coro_row.append('Y' if coronagraph.check_driver() else 'N')
				coro_row.append('Y' if coronagraph.check_solution() else 'N')
				survey.writerow(coro_row)

		os.chmod(fname, 644)

class PorAPLC(object):
	_default_parameters = {
		'pupil': {
			'filename': 'NoFilename'
			},
		'focal_plane_mask': {
			'radius': 4.0,
			'num_pix': 50,
			'grayscale': True,
			'field_stop_radius': -1.0
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
			'initial_throughput_estimate': 1
			},
		'solver': {
			'num_threads': 0,
			'crossover': 0,
			'method': 2
			}
		}
	
	def __init__(self, identifier, parameters, file_organization):
		self.identifier = identifier
		self.parameters = parameters
		self.file_organization = file_organization
	
	def write_driver(self, overwrite=False):
		if self.check_driver() and not overwrite:
			print('Driver already exists and will not be overwritten.')
			return

		with open('por_aplc_driver_template.py') as template_file:
			driver_template = template_file.read()

		fname = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')
		with open(fname, 'w') as output_file:
			output_file.write(driver_template.format())

	@classmethod
	def write_static_files(self):
		pass
	
	def check_input_files(self):
		return False
	
	def check_driver(self):
		fname_driver = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')
		return os.path.exists(fname_driver)
	
	def check_solution(self):
		fname_sol = os.path.join(self.file_organization['solution_dir'], self.get_identifier() + '.fits')
		return os.path.exists(fname_sol)
	
	def get_driver_command(self):
		fname_driver = os.path.join(self.file_organization['drivers_dir'], self.identifier + '.py')
		fname_log = os.path.join(self.file_organization['log_dir'], self.identifier + '.log')

		return '{:s} {:s} &> {:s}'.format(os.__file__, fname_driver, fname_log)

	def get_identifier(self):
		return self.identifier
	
	def run_optimization(self):
		os.system(self.get_driver_command())
	
	def run_analysis(self, overwrite=False):
		analysis = inspect.getmembers(por_aplc_analysis, inspect.isfunction)
		analysis = [x for _, x in sorted(zip([name for name, function in analysis], analysis)) if _.startswith('analyze')]
		
		fname = os.path.join(self.file_organization['analysis_dir'], self.get_identifier() + '.pdf')
		if os.path.exists(fname):
			print('Analysis already exists and is not overwritten.')
			return
		
		pdf = PdfPages(fname)
		self.metrics = {}

		for name, function in analysis:
			if not name.startswith('analyze_'):
				continue

			fname_sol = os.path.join(self.file_organization['solution_dir'], self.get_identifier() + '.fits')
			
			res = function(fname_sol, pdf)
			self.metrics.update(res)
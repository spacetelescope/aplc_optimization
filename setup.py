from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(name='aplc_optimization',
      version='1.0',
      description='An Apodized Pupil Lyot Coronagraph design survey toolkit.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      license="BSD 3-Clause License",
      url='https://github.com/spacetelescope/aplc_optimization',
      author='Space Telescope Science Institute, Segmented Coronagraph Design and Analysis team',
      author_email='help@stsci.edu',
      project_urls={
            'Documentation': 'aplc_optimization.readthedocs.io',
            'Funding': 'https://exoplanets.nasa.gov/exep/technology/SCDA/',
            'Source': 'https://github.com/spacetelescope/aplc_optimization',
            'Tracker': 'https://github.com/spacetelescope/aplc_optimization/issues'
      },
      packages=find_packages())
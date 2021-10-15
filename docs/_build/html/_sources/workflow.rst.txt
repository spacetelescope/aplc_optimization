.. _workflow:

Design Survey Workflow
========================

Using the ``aplc_optimization`` toolkit, the workflow for producing a design survey proceeds as follows:

1. Create a new launcher script
-------------------------
Within ``aplc_optimization``, *Launcher scripts* are the format by which coronagraph design surveys are defined and executed. Template launcher scripts are provided
for HiCAT-, LUVOIR- and GPI-like instruments (named :mod:`do_hicat_template.py`, :mod:`do_luvoir_template.py` and :mod:`do_gpi_template.py``, respectively).

In order to initiate a new survey, make a copy of the appropriate launcher template and re-name the file with the following naming convention:

    - do_<*instrument*>_<*survey*>_<*machine*>.py

where *<instrument>* is the name of the APLC instrument for which you are are performing the survey; '*<survey>*' is a name
descriptive of the survey that you intend to run; and '*<machine>*' is the name of the machine on which the survey will be run on. For example,
a launcher file named ':mod:`do_luvoir_BW10_small_local.py`' designates the BW10 small design run on your local machine for a LUVOIR-like instrument.

2. Define a set of design parameters to survey
-------------------------------------------
Inside the launcher script define a set of design parameters to survey. Any unspecified parameters are set to reasonable default values.


3. Run the launcher script
------------------------

Run the launcher script using the following command in terminal::

    (aplc_optimization) $ python do_<instrument>_<survey>_<machine>.py

.. _output:

If all goes well, this should provide you with the following output::

    This survey has # design parameter combinations.
    # parameter are varied:

    File organization:
    {'analysis_dir': '../aplc_optimization/surveys/<instrument>_<survey>_<machine>/analysis',
     'drivers_dir': '../aplc_optimization/surveys/<instrument>_<survey>_<machine>/drivers',
     'input_files_dir': '../aplc_optimization/masks',
     'log_dir': '../aplc_optimization/surveys/<instrument>_<survey>_<machine>/logs',
     'solution_dir': '../aplc_optimization/surveys/<instrument>_<survey>_<machine>',
     'survey_dir': '../aplc_optimization/surveys/<instrument>_<survey>_<machine>'}


    All input files exist? False
    All drivers exist? False
    All solutions exist? False


.. _file-struct:

4. Inspect the products
-----------------------

- **Input masks:** Once launched, the toolkit will first run a "check" routine to verify whether the necessary static input files defining the telescope aperture and lyot stop masks are in place. If they are not, the program will call the corresponding input file generation script and write them to disk in the :mod:`.../aplc_optimization/masks/` directory.
- **Drivers:** After checking that the necessary input files are in place, the toolkit will write a batch of driver scripts (one for each design parameter combination) and stores them in the :mod:`drivers/` sub-directory, inside the survey's base directory (:mod:`.../aplc_optimization/surveys/<instrument>_<survey>_<machine>/`). Each driver script is a re-usable python script storing the parameters and file organization defined in the launcher, which in turn calls the optimizer.
- **Logs:** While the optimizer runs, an automatically produced :mod:`.log` file containing a record of events is written to file and stored in the :mod:`logs/` sub-directory (:mod:`.../aplc_optimization/surveys/<instrument>_<survey>_<machine>/logs/`).
- **Solutions:** Once each linear optimization program completes, each apodizer solution is written to file in the :mod:`solutions/` sub-directory (:mod:`.../aplc_optimization/surveys/<instrument>_<survey>_<machine>/solutions`).
- **Analysis files:** For each subsequent apodizer solution the toolkit calls an analysis script that produces a number of analysis products and writes them to a :mod:`.pdf` file in the :mod:`analysis` sub-directory (:mod:`.../aplc_optimization/surveys/<instrument>_<survey>_<machine>/analysis`).


5. Run the Analysis notebook
-----------------------------
In addition to the analysis :mod:`.pdf` file created automatically created for each design survey, the toolkit also
provides a Python notebook interface with which more in-depth analyses can be performed for each apodizer solution.
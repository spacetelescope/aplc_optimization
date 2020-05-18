.. _workflow:

Workflow
============================

.. warning::

   Pardon our dust, this space is under construction.

Create a launcher file
----------------------------

The launcher file is a script that creates
set the survey parameters,
launcher templates are provided for LUVOIR and HICAT, these scripts are not designed to be run,m but are provided as templates only.

~~~ About launcher file. ~~~

Make a copy of the appropriate launcher template---depending on the desired instrument/ telescope (e.g. if you are optimizing an apodizer for
a LUVOIR-type telescope, make a copy of the '``do_luvoir_template.py``' file---and rename the copy with the following naming convention::

    do_<instrument_name>_<survey_name>_<machine_name>.py

where <instrument_name> is the name of the instrument/ telescope for which you are apodizing; <survey_name> is the name
of the survey that you are running; and <machine_name> is the name of the machine that the survey will be run on. For example,
a launcher file named '``do_luvoir_BW1o_small_telserv3.py``' designates the BW10 small design run on telserv3 for a LUVOIR-type telescope.

define the input file generation parameters
--------------------------------------------
The next step is to define the parameters for the generation of the input files.

.. py:function:: input_files_dict[]

   :py:obj:`input_files_dict` is a python dictionary containing the input parameters for the telescope aperture and lyot stop file generation---these parameters are to be defined as follows:

   :param str directory: (*value pre-set and not to be altered*) the name of the folder located in *'aplc-optmization/MASKS/'* within which the input files are stored.
   :param int N: set to the number of pixels in the input and output arrays
   :param int oversamp: set the oversampling level (i.e. the number of grey levels).
   .. py:function:: aperture[]

        :py:obj:`aperture` is a nested dictionary containing the input parameters for the generation of the telescope aperture. These parameters should be defined as follows:

        :param int seg_gap_pad: set the arbitrary padding of gap size to represent gaps on smaller arrays - this effectively makes the gaps larger and the segments smaller to preserve the same segment pitch.
   .. py:function:: lyot_stop[]

    :py:obj:`lyot_stop` is a nested dictionary containing the input parameters for the Lyot stop generation.  These parameters should be defined as follows:

    :param float lyot_ref_diam: define the Lyot reference diameter.
    :param bool ls_spid: opt whether to include secondary mirror support structure in the aperture.
    :param int ls_spid_ov: set to the factor by which to oversize the spiders.
    :param float LS_ID: set to the fractional size of the Lyot stop's circular central obstruction as fraction of the reference diameter (**LS_ref_diam**)
    :param float LS_OD: set to the fractional size of the circular outer edge of the Lyot stop as fraction of the reference diameter (**LS_ref_diam**).

The values set in :py:obj:'input_files_dict' dictionary are passed to the :py:func:`LUVOUR_input_gen()` function to create the input telescope pupil and lyot stop, which are saved as fits files with the
filenames ``pup_filename`` and ``ls_filenames``, located in the ``'MASKS/<directory>'`` folder.

Define the design survey parameters
-------------------------------------

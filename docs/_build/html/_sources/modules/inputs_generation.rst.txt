Inputs Generation
#############################################

Pupil and Lyot stop generation scripts, for use as inputs to ``aplc``.

HiCAT Inputs Generation
========================




LUVOIR Inputs Generation
========================

Functions
LUVOIR_inputs_gen


.. py:function:: LUVOIR_inputs_gen(input_files_dict)

    Called by a *launcher template*, the function creates a LUVOIR telescope aperture and lyot stop using the
    input parameters provided in the ``input_files_dict`` dictionary, and saves them to file.

Using the ``aperture`` parameters from the input_files_dict prodided

LUVOIR aperture and Lyot stop generation

    :param input_files_dict: input file dictionary
    :type input_files_dict: dict
    :return: ``pup_fname`` and ``ls_fname``, the names of the aperture and Lyot stop fits file, respetively.
    :rtype: list
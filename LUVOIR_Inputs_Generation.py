from pathlib import Path

import numpy as np
from astropy.io import fits
from hcipy import *


def LUVOIR_inputs_gen(input_files_dict):
    """ Generate LUVOIR-A aperture and Lyot stop(s) for the APLC coronagraph.

    Parameters
    ----------
    input_files_dict: dict
        A dictionary of input parameters provided in `do_luvoir` launcher script.

    Returns
    -------
    pup_filename: str
        The name of the generated LUVOIR A aperure FITS file.
    ls_filename: list
        The name(s) of the generated Lyot stop FITS file(s).
    """


    filepath =      input_files_dict['directory']   # directory in which the constructed FITS files are stored
    N =             input_files_dict['N']   # number of pixels in the input (aperture and lyot stop) arrays

    # Aperture parameters
    oversamp = 		input_files_dict['oversamp']    # oversampling factor in `evaluate_supersampled()` [hcipy/field/util.py].
    gap_padding = 	input_files_dict['aperture']['seg_gap_pad'] # arbitratry padding of gap size to represent gaps on small arrays

    # Lyot stop parameters
    lyot_ref_diam = input_files_dict['lyot_stop']['lyot_ref_diam']  # diameter used to reference LS_ID and LS_OD against
    LS_SPID = 		input_files_dict['lyot_stop']['ls_spid']    # flag for inclusion or exclusion of struts in lyot stop
    ls_spid_ov = 	input_files_dict['lyot_stop']['ls_spid_ov'] # spider oversize scale
    LS_OD = 		input_files_dict['lyot_stop']['LS_OD']  # Lyot stop inner diameter(s), relative to inscribed circle
    LS_ID = 		input_files_dict['lyot_stop']['LS_ID']  # Lyot stop outer diameter(s), relative to inscribed circle

    pup_filename = filepath + 'TelAp_LUVOIR_gap_pad{0:02d}_bw_ovsamp{1:02d}_N{2:04d}.fits'.format(gap_padding, oversamp,
                                                                                                  N)
    pup_filename_indexed = filepath + 'TelAp_LUVOIR_gap_pad{0:02d}_bw_ovsamp{1:02d}_N{2:04d}_indexed.fits'.format(
        gap_padding, oversamp, N)

    grid = make_pupil_grid(N)

    config = Path('masks/' + pup_filename)

    '''
    Checks if aperture file already exists, otherwise creates LUVOIR aperture by calling 
    make_a_luvoir_aperture() [hcipy/aperture/realistic.py]
    '''
    if config.is_file():
        print('{0:s} exists'.format('masks/' + pup_filename))
    else:
        LUVOIR_ap, header = make_luvoir_a_aperture(gap_padding, return_header=True)  # header = True
        LUVOIR_ap_indexed, _, segment_positions = make_luvoir_a_aperture(gap_padding, return_header=True,
                                                                         segment_transmissions=np.arange(1, 121),
                                                                         return_segments=True)
                                                                        # return_segment_positions=True, header=True
        pupil = evaluate_supersampled(LUVOIR_ap, grid, oversamp)
        pupil_indexed = evaluate_supersampled(LUVOIR_ap_indexed, grid, 1)

        header['EFF_GAP'] = header['GAP_PAD'] * header['SEG_GAP']

        hdr = fits.Header()
        hdr.set('TELESCOP', header['TELESCOP'])
        hdr.set('D_CIRC', header['D_CIRC'], 'm: circumscribed diameter')
        hdr.set('D_INSC', header['D_INSC'], 'm: inscribed diameter')
        hdr.set('SEG_F2F', header['SEG_F2F_D'], 'm: actual segment flat-to-flat diameter')
        hdr.set('SEG_GAP', header['SEG_GAP'], 'm: actual gap size between segments')
        hdr.set('GAP_PAD', header['GAP_PAD'], 'arbitrary multiplicative padding of gaps')
        hdr.set('EFF_GAP', header['EFF_GAP'], 'm: effective gap size after padding')
        hdr.set('STRUT_W', header['STRUT_W'], 'm: actual support strut width')
        hdr.set('STRUT_ST', header['STRUT_ST'], 'm: lower spider starting point d from center')
        hdr.set('STRUT_AN', header['STRUT_AN'], 'deg: angle lower spiders offset from vertical')
        hdr.set('NORM', header['NORM'], 'normalization keyword, OD scaled to 1 by Dcirc')
        hdr.set('SEG_TRAN', header['SEG_TRAN'], 'The transmission for each of the segments')
        hdr.set('EDGE', 'bw', 'black and white, or grey pixels')
        hdr.set('PROV', header['PROV'])

        hdr_indexed = hdr.copy()
        seg_xys = segment_positions.points
        for segment in range(seg_xys.shape[0]):
            hdr_indexed.set('seg' + str(segment + 1) + '_x', segment_positions.x[segment], 'x-position segment ' + str(segment + 1))
            hdr_indexed.set('seg' + str(segment + 1) + '_y', segment_positions.y[segment], 'y-position segment ' + str(segment + 1))

        fits.writeto('masks/' + pup_filename, pupil.shaped, hdr, overwrite=True)
        fits.writeto('masks/' + pup_filename_indexed, pupil_indexed.shaped, hdr_indexed, overwrite=True)
        print('{0:s} has been written to file'.format('masks/' + pup_filename))

    ls_filenames = []

    for ls_od in LS_OD:
        for ls_id in LS_ID:
            # ls_id = (0.19/0.937)*ls_od

            # filename key for struts or no struts
            # black and white or grey pixels
            if ls_id == 0:
                LS_SPID = False

            if oversamp == 1:
                edge = 'bw'
            elif oversamp > 1:
                edge = 'gy'

            if LS_SPID:
                strut_key = 'struts'
            else:
                strut_key = 'no_struts'

            if LS_SPID:
                ls_filename = filepath + 'LS_LUVOIR_ID{0:04d}_OD{1:04d}_{2:s}_pad{3:02d}_{4:s}_ovsamp{5:d}_N{6:04d}.fits'.format(
                    int(ls_id * 1000), \
                    int(ls_od * 1000), \
                    strut_key, ls_spid_ov, \
                    edge, oversamp, N)
            else:
                ls_filename = filepath + 'LS_LUVOIR_ID{0:04d}_OD{1:04d}_{2:s}_{3:s}_ovsamp{4:d}_N{5:04d}.fits'.format(
                    int(ls_id * 1000), \
                    int(ls_od * 1000), \
                    strut_key, edge, \
                    oversamp, N)

            config = Path('masks/' + ls_filename)
            if config.is_file():
                print('{0:s} exists'.format('masks/' + ls_filename))

            else:

                LUVOIR_ls, ls_header = make_luvoir_lyot_stop(inner_diameter_fraction=ls_id,
                                                             outer_diameter_fraction=ls_od, lyot_reference_diameter=lyot_ref_diam,
                                                             spider_oversize=ls_spid_ov, with_spiders=LS_SPID,
                                                             return_header=True)
                ### previously make_luvour_a_lyot_stop(ls_id, ls_od, lyot_ref_diam, spid_oversize=ls_spid_ov, spiders=LS_SPID, header = True)
                lyot_stop = evaluate_supersampled(LUVOIR_ls, grid, oversamp)

                # header.update(ls_header)
                hdr = fits.Header()
                ls_header['OVERSAMP'] = oversamp
                ls_header['EDGE'] = edge

                hdr.set('TELESCOP', ls_header['TELESCOP'])
                hdr.set('D_CIRC', ls_header['D_CIRC'], 'm: circumscribed diameter')
                hdr.set('D_INSC', ls_header['D_INSC'], 'm: inscribed diameter')
                hdr.set('LS_REF_D', ls_header['LS_REF_D'], 'm: used to reference given LS id and od')
                hdr.set('LS_ID', ls_header['LS_ID'], 'LS inner d, fraction of LS_REF_D')
                hdr.set('LS_OD', ls_header['LS_OD'], 'LS outer d, fraction of LS_REF_D')

                if LS_SPID:
                    hdr.set('STRUT_W', ls_header['STRUT_W'], 'm: actual support strut width')
                    hdr.set('STRUT_ST', ls_header['STRUT_ST'], 'm: lower spider starting point d from center')
                    hdr.set('STRUT_AN', ls_header['STRUT_AN'], 'deg: angle lower spiders offset from vertical')
                    hdr.set('STRUT_P', ls_header['STRUT_P'], 'spider padding factor')

                hdr.set('NORM', ls_header['NORM'], 'normalization keyword, OD scaled to 1 by Dcirc')
                hdr.set('EDGE', ls_header['EDGE'], 'black and white, or grey pixels')
                hdr.set('OVERSAMP', ls_header['OVERSAMP'], 'oversampling factor, # grey levels')

                fits.writeto('masks/' + ls_filename, lyot_stop.shaped, hdr, overwrite=True)
                print('{0:s} has been written to file'.format('masks/' + ls_filename))

            ls_filenames.append(ls_filename)
            # remove duplicates
            s = set(ls_filenames)
            ls_filenames = list(s)

        return pup_filename, ls_filenames


def make_luvoir_lyot_stop(normalized=False, with_spiders=False, spider_oversize=1, lyot_reference_diameter=13.5, inner_diameter_fraction=0.2,
                            outer_diameter_fraction=0.9, return_header=False):
    '''Make a LUVOIR-A Lyot stop for the APLC coronagraph.

    Parameters
    ----------
    normalized : boolean
        If this is True, the pupil diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 15.0 meters.
    with_spiders : boolean
        Include the secondary mirror support structure in the aperture.
    lyot_reference_diameter : scaler
        The diameter used to reference LS id and od against, by default equal to the pupil inscribed diameter.
    inner_diameter_fraction : scalar
        The fractional size of the lyot stop inner diameter(s) as a fraction of the inscribed circle diameter.
    outer_diameter_fraction : scalar
        The fractional size of the lyot stop outer diameter(s) as a fraction of the inscribed circle diameter.
    spider_oversize : scalar
        The factor by which to oversize the spiders compared to the LUVOIR-A aperture spiders.
    return_header : boolean
        If this is True, a header will be returned giving all important values for the
        created aperture for reference.

    Returns
    -------
    lyot_stop : Field generator
        A field generator for the Lyot stop.
    header : dict
        A dictionary containing all important values for the created aperture. Only returned
        if `return_header` is True.
    '''
    pupil_diameter = 15.0  # m actual circumscribed diameter, used for lam/D calculations, other measurements normalized by this diameter
    pupil_diameter_inscribed = 13.5
    spider_width = 0.150  # m actual strut size
    lower_spider_angle = 12.7  # deg angle at which lower spiders are offset from vertical
    spid_start = 0.30657  # m spider starting point offset from center of aperture

    outer_D = pupil_diameter * outer_diameter_fraction #re-normalize the LS OD against circumscribed pupil diameter
    inner_D = pupil_diameter * inner_diameter_fraction #re-normalize the LS ID against circumscribed pupil diameter
    pad_spid_width = spider_width * spider_oversize

    #lyot_reference_diameter = pupil_diameter

    ls_header = {'TELESCOP': 'LUVOIR A', 'D_CIRC': pupil_diameter, 'D_INSC': pupil_diameter_inscribed,
                 'LS_ID': inner_diameter_fraction, 'LS_OD': outer_diameter_fraction,
                 'LS_REF_D': lyot_reference_diameter, 'NORM': normalized, 'STRUT_ST': spid_start}

    if with_spiders:
        ls_header['STRUT_W'] = spider_width
        ls_header['STRUT_AN'] = lower_spider_angle
        ls_header['STRUT_P'] = spider_oversize

    if normalized:
        outer_D /= pupil_diameter
        inner_D /= pupil_diameter
        pad_spid_width /= pupil_diameter
        spid_start /= pupil_diameter

    outer_diameter = circular_aperture(outer_D)
    central_obscuration = circular_aperture(inner_D)

    if with_spiders:
        spider1 = make_spider_infinite([0, 0], 90, pad_spid_width)
        spider2 = make_spider_infinite([spid_start, 0], 270 - lower_spider_angle, pad_spid_width)
        spider3 = make_spider_infinite([-spid_start, 0], 270 + lower_spider_angle, pad_spid_width)

    def aper(grid):
        result = outer_diameter(grid) - central_obscuration(grid)

        if with_spiders:
            result *= spider1(grid) * spider2(grid) * spider3(grid)

        return result

    if return_header:
        return aper, ls_header

    return aper
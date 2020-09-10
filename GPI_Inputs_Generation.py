from pathlib import Path
from astropy.io import fits
import gpipsfs


def GPI_inputs_gen(input_files_dict):
    """ Generate GPI Primary and Lyot stop(s) for the APLC coronagraph.

    Parameters
    ----------
    input_files_dict: dict
        A dictionary of input parameters provided in `do_gpi...` launcher script.

    Returns
    -------
    pup_filename: str
        The name of the generated GPI aperure FITS file.
    ls_filename: list
        The name(s) of the generated Lyot stop FITS file(s).
    """

    #   filepath =      input_files_dict['directory']   # directory in which the constructed FITS files are stored
    N = input_files_dict['N']  # number of pixels in the input (aperture and lyot stop) arrays
    ls_name = input_files_dict['lyot_stop']['ls_name']
    ls_tabs = input_files_dict['lyot_stop']['ls_tabs']

    pup_filename = 'GPI_primary_N{0:04d}.fits'.format(N)
    ls_filename = 'LS_GPI_{0:s}_N{1:04d}.fits'.format(ls_name, N)

    config = Path(pup_filename)
    '''
    Checks if aperture file already exists, otherwise creates GPI aperture by calling 
    GeminiPrimary() [gpipsfs/gpipsfs/main.py]
    '''

    if config.is_file():
        print('{0:s} exists'.format(pup_filename))
    else:
        GPI_primary = gpipsfs.GeminiPrimary().sample(npix=N)

        hdr = fits.Header()
        hdr.set('TELESCOP', 'GPI')
        hdr.set('NAME', 'GEMINI SOUTH PRIMARY')
        hdr.set('D_CIRC', 7.701, 'm: projected diameter of baffle on M2')
        hdr.set('D_INSC', 1.2968, 'm: projected diameter of M2 inner oversized hole')
        hdr.set('STRUT_W', '0.01, 0.014', 'm: actual support strut widths (laser vane is slightly thicker)')
        hdr.set('STRUT_ST', -0.2179, 'm: angle lower spiders offset from y')
        hdr.set('STRUT_AN', 226.9, 'deg: angle lower spiders offset from vertical')

        fits.writeto(pup_filename, GPI_primary, header=hdr)
        print('{0:s} has been written to file'.format(pup_filename))

    config = Path(ls_filename)
    '''
    Checks if lyot stop file already exists, otherwise creates GPI aperture by calling 
    GPI_LyotMask() [gpipsfs/gpipsfs/main.py]
    '''

    if config.is_file():
        print('{0:s} exists'.format(ls_filename))
    else:
        Lyot_mask = gpipsfs.GPI_LyotMask(ls_name).sample(npix=N)

        ls_atr = gpipsfs.GPI_LyotMask(ls_name)
        hdr = fits.Header()
        hdr.set('TELESCOP', 'GPI')
        hdr.set('LS_MASK', ls_atr.name, "name of Lyot mask")
        hdr.set('LS_ID', ls_atr.inner_radius, 'm: Lyot mask inner radius')
        hdr.set('LS_OD', ls_atr.outer_radius, 'm: Lyot mask outer radius')
        if ls_tabs:
            hdr.set('NTABS', ls_atr.ntabs, "Number of bad actuator tabs")
            hdr.set('TAB_DIAM', ls_atr.tabradius * 2, "Tab diameter")
        hdr.set('PUP_DIAM', 9.825e-3, 'm: Gemini pupil size in Lyot plane (without undersizing)')
        hdr.set('MAG', 7.7701 / .009825, 'm: magnification between primary & lyot')
        hdr.set('ST_ANG', str(ls_atr.support_angles), 'deg: spider angles')

        fits.writeto(ls_filename, Lyot_mask, header=hdr, overwrite=True)
        print('{0:s} has been written to file'.format(ls_filename))

    return pup_filename, ls_filename

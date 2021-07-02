import os
import re
import time

import asdf
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from hcipy import *


def create_coronagraph(solution_filename):
    """Create APLC object from solution file.

    Parameters
    ----------
    solution_filename: string
        The file path to the solution file.

    Returns
    -------
    pupil : Field
        The telescope pupil.
    focal_plane_mask : Field
        The focal plane mask of the APLC. The grid is assumed to be in lambda_0/D.
    lyot_stops : list of Fields
        A list of the Lyot stops used in the optimization problem.
    apodizer: Field
        The optimized apodizer for the APLC.
    parameters: dict
        A dictionary of the survey parameters.
    file_organization: dict
        A dictionary of the file organization structure of the survey.
    """

    # Open the solution as ASDF file and read the tree of data
    solution = asdf.open(solution_filename)
    parameters = solution.tree['parameters']
    file_organization = solution.tree['file_organization']
    apodizer = solution.tree['apodizer']
    pup_fname = parameters['pupil']['filename']
    fpm_radius = parameters['focal_plane_mask']['radius']
    fpm_num_pix = parameters['focal_plane_mask']['num_pix']
    fpm_grayscale = parameters['focal_plane_mask']['grayscale']
    ls_fname = parameters['lyot_stop']['filename']
    ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
    ls_num_stops = parameters['lyot_stop']['num_lyot_stops']

    if not os.path.isabs(pup_fname):
        pup_fname = os.path.join(file_organization['input_files_dir'], pup_fname)
    if not os.path.isabs(ls_fname):
        ls_fname = os.path.join(file_organization['input_files_dir'], ls_fname)

    pupil = read_fits(pup_fname)
    num_pix = pupil.shape[0]
    pupil_grid = make_uniform_grid(num_pix, [1, 1])
    pupil = Field(pupil.ravel(), pupil_grid)
    apodizer = Field(apodizer.ravel(), pupil_grid)

    try:
        lyot_stops = [Field(read_fits(ls_fname.format(i)).ravel(), pupil_grid) for i in range(ls_num_stops)]
    except:
        lyot_stop = Field(read_fits(ls_fname).ravel(), pupil_grid)

        # Building Lyot stops according to alignment tolerance
        if ls_num_stops in [1, 5, 9]:
            lyot_stops = [lyot_stop]
        else:
            lyot_stops = []

        if ls_num_stops in [4, 5, 8, 9]:
            lyot_stop_pos_x = np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1).ravel()
            lyot_stop_neg_x = np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1).ravel()
            lyot_stop_pos_y = np.roll(lyot_stop.shaped, ls_alignment_tolerance, 0).ravel()
            lyot_stop_neg_y = np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 0).ravel()

            lyot_stops.extend([lyot_stop_pos_x, lyot_stop_neg_x, lyot_stop_pos_y, lyot_stop_neg_y])

        if ls_num_stops in [8, 9]:
            lyot_stop_pos_x_pos_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1),
                                            ls_alignment_tolerance, 0).ravel()
            lyot_stop_pos_x_neg_y = np.roll(np.roll(lyot_stop.shaped, ls_alignment_tolerance, 1),
                                            -ls_alignment_tolerance, 0).ravel()
            lyot_stop_neg_x_pos_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1),
                                            ls_alignment_tolerance, 0).ravel()
            lyot_stop_neg_x_neg_y = np.roll(np.roll(lyot_stop.shaped, -ls_alignment_tolerance, 1),
                                            -ls_alignment_tolerance, 0).ravel()

            lyot_stops.extend(
                [lyot_stop_pos_x_pos_y, lyot_stop_pos_x_neg_y, lyot_stop_neg_x_pos_y, lyot_stop_neg_x_neg_y])

        # Build focal plane mask
    q_foc = fpm_num_pix / (fpm_radius * 2)
    x_foc = (np.arange(fpm_num_pix) + 0.5 - fpm_num_pix / 2) / q_foc
    focal_mask_grid = CartesianGrid(RegularCoords(1.0 / q_foc, [fpm_num_pix, fpm_num_pix], x_foc.min()))

    if fpm_grayscale:
        focal_plane_mask = 1 - evaluate_supersampled(circular_aperture(2 * fpm_radius), focal_mask_grid, 8)
    else:
        focal_plane_mask = 1 - circular_aperture(2 * fpm_radius)(focal_mask_grid)

    return pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization


def analyze_aplc_design_summary(solution_filename, pdf=None):
    """Generate a summary table for the APLC design.

    Parameters
    ----------
    solution_filename: string
        The path to the apodizer solution file.
    pdf: bool
        Whether to save the table to file or display on screen.

    Returns
    -------
    An empty dictionary.
    """

    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    core_throughput = analyze_throughput(solution_filename, pdf=None)
    # FPM
    fpm_radius = parameters['focal_plane_mask']['radius']
    fpm_num_pix = parameters['focal_plane_mask']['num_pix']
    fpm_grayscale = parameters['focal_plane_mask']['grayscale']
    # Image
    img_contrast = parameters['image']['contrast']
    img_iwa = parameters['image']['iwa']
    img_owa = parameters['image']['owa']
    img_num_wavelengths = parameters['image']['num_wavelengths']
    img_bandwidth = parameters['image']['bandwidth'] * 100
    # Lyot stop
    ls_fname = parameters['lyot_stop']['filename']
    ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
    # pupil
    pup_fname = parameters['pupil']['filename']
    N = parameters['pupil']['N']

    # Extract solution filename without pathing
    sol_fname = solution_filename[solution_filename.find('solutions/') + 10:]

    # Extract instrument name, gap padding and oversampling info from the pupil filename
    if 'LUVOIR' in pup_fname:
        instrument = 'LUVOIR'
    elif 'HiCAT' in pup_fname:
        instrument = 'HiCAT'
    elif 'GPI' in pup_fname:
        instrument = 'GPI'

    regex = re.compile(r'\d+')
    pup_vals = regex.findall(pup_fname)
    if instrument == 'LUVOIR':
        seg_gap_pad, oversamp = int(pup_vals[0]), int(pup_vals[1])
    elif instrument == 'HiCAT':
        seg_gap_pad = 'n/a'
        oversamp = int(pup_vals[0])
    elif instrument == 'GPI':
        seg_gap_pad = 'n/a'
        oversamp = 'n/a'

    # Extract lyot stop inner and outer diameter from lyot stop filename
    ls_vals = regex.findall(ls_fname)
    ls_id, ls_od = int(ls_vals[0]) / 1000, int(ls_vals[1]) / 1000

    if fpm_grayscale is True:
        gs = ' (grayscale)'
    else:
        gs = ''

    fname_pup = parameters['pupil']['filename']
    fname_ls = parameters['lyot_stop']['filename']

    if not os.path.isabs(fname_pup):
        fname_pup = os.path.join(file_organization['input_files_dir'], fname_pup)
    if not os.path.isabs(fname_ls):
        fname_ls = os.path.join(file_organization['input_files_dir'], fname_ls)
    fname_apod = solution_filename

    TelAp = fits.getdata(fname_pup)
    LS = fits.getdata(fname_ls)
    A = fits.getdata(fname_apod)

    maximum_integrated_throughput = ((TelAp * LS * A).sum() / (TelAp * LS).sum()) ** 2

    fig = plt.figure(dpi=100)

    col_labels = ['$\mathbf{APLC \ Design \ Summary}$ \n\n', '']

    ax = fig.add_subplot(1, 1, 1)

    table_data = [
        ['Instrument', '$\mathtt{' + str(instrument) + '}$'],
        ['nPup', '$\mathtt{' + str(N) + ' \ x \ ' + str(N) + ' \ pixels}$'],
        ['Coronagraphic throughput (transmitted energy)',
         '$\mathtt{' + str(round(float(maximum_integrated_throughput), 4)) + '}$'],
        ['Core throughput (encircled energy)', '$\mathtt{' + str(round(float(core_throughput), 4)) + '}$'],
        ['Lyot stop inner diamater (% of inscribed circle)', '$\mathtt{' + str(ls_id) + '}$'],
        ['Lyot stop outer diameter (% of inscribed circle)', '$\mathtt{' + str(ls_od) + '}$'],
        ['Bandpass', '$\mathtt{' + str(img_bandwidth) + '\%}$'],
        ['# wavelengths', '$\mathtt{' + str(img_num_wavelengths) + '}$'],
        ['FPM radius' + str(gs), '$\mathtt{' + str(round(fpm_radius, 4)) + ' \ \lambda/D}$'],
        ['nFPM', '$\mathtt{' + str(fpm_num_pix) + ' \ pixels}$'],
        ['IWA ' + u'\u2014' + ' OWA',
         '$\mathtt{' + str(float(img_iwa)) + '\u2014' + str(float(img_owa)) + ' \ \lambda/D}$'],
        ['Contrast constraint', '$\mathtt{10^{-' + str(img_contrast) + '}}$'],
        ['Lyot Stop alignment tolerance', '$\mathtt{' + str(ls_alignment_tolerance) + ' \ pixels}$'],
        ['$\mathit{Input \ Files:}$', ''],
        [r'       $\triangleright \ \ Pupil \ file:$         ' + str(pup_fname), ''],
        [r'       $\triangleright \ \ Lyot \ stop \ file:$ ' + str(ls_fname), ''],
        ['$\mathit{Solution \ File:}$', ''],
        [r'$\triangleright$ ' + sol_fname,
         time.ctime(os.path.getmtime(solution_filename))]]  # Last modified, date & time

    table = ax.table(cellText=table_data, colLabels=col_labels, colWidths=[0.8, 0.3], cellLoc='left', edges='open',
                     loc='center')
    table.set_fontsize(20)
    table.scale(1.2, 1.5)
    plt.axis('off')

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {}


def analyze_contrast(solution_filename, pdf=None):
    """Monochromatic PSF evaluation for the APLC design.

    Parameters
    ----------
    solution_filename:
        The location of the apodizer solution file.
    pdf: bool
        Whether to save the table to file or display on screen.

    Returns
    -------
    The monochromatic normalized irradiance image and radial profile.
    """
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    iwa = parameters['image']['iwa']
    owa = parameters['image']['owa']
    contrast = parameters['image']['contrast']
    radius_fpm = parameters['focal_plane_mask']['radius']
    num_wavelengths = parameters['image']['num_wavelengths']
    bandwidth = parameters['image']['bandwidth']
    lyot_stop = lyot_stops[0]

    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
    focal_grid = make_focal_grid(8, owa * 1.2)  # make_focal_grid(q, fov) - grid for a focal plane
    prop = FraunhoferPropagator(pupil.grid, focal_grid)

    if num_wavelengths == 1:
        monochromatic = True
    else:
        monochromatic = False

    fig, ax = plt.subplots()
    if monochromatic:
        wf = Wavefront(pupil * apodizer)
        img = prop(coro(wf)).intensity
        img_ref = prop(Wavefront(apodizer * pupil * lyot_stop)).intensity

        # For the first subplot
        plt.title('Monochromatic Normalized Irradiance')
        caption = '\n \n $\mathit{ Monochromatic \ on-axis \ PSF \ in \ log \ irradiance,}$ \n ' \
                  '$\mathit{normalized \ to \ the \ peak \ irradiance \ value.}$'
    else:
        img = 0
        img_ref = 0
        imgs = []
        wavelengths = np.linspace(1 - bandwidth / 2, 1 + bandwidth / 2, 11)

        for wl in wavelengths:
            wf = Wavefront(apodizer * pupil, wl)
            psf = prop(coro(wf)).power
            psf_ref = prop(Wavefront(apodizer * pupil * lyot_stop, wl)).power

            img += psf
            img_ref += psf_ref

            imgs.append(psf / psf_ref.max())

        img /= 11
        img_ref /= 11

        # For the first subplot
        plt.title('Broadband Normalized Irradiance')
        caption = '\n \n $\mathit{ On-axis \ PSF \ in \ log \ irradiance,}$ \n ' \
                  '$\mathit{normalized \ to \ the \ peak \ irradiance \ value.}$'

    plt.figure(dpi=100)
    imshow_field(np.log10(img / img_ref.max()), vmin=-contrast - 0.5, vmax=-contrast + 3.5, cmap='inferno')
    cbar = plt.colorbar()
    cbar.set_label('log irradiance')
    plt.xlabel('Angular coordinate ($\lambda_0/D$)' + caption, fontsize=10)
    plt.ylabel('Angular coordinate ($\lambda_0/D$)', fontsize=10)
    plt.tick_params(axis='both', labelsize=8)
    plt.tight_layout()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    plt.figure(figsize=(10, 8), dpi=100)
    r, profile, std_profile, n_profile = radial_profile(img / img_ref.max(), 0.2)
    plt.plot(r, profile, c='k', lw=2)

    if monochromatic:
        plt.title('\n Monochromatic Normalized Irradiance (Radial Average)')
        caption_radial = '\n \n \n $\mathit{Monochromatic \ on-axis \ PSF \ azimuthally \ averaged \ over \ angular}$' \
                         '\n $\mathit{seperations \ }$' + str(round(min(r), 4)) + '-' + str(round(max(r), 4)) + \
                         '$\mathit{ \ λ/D, \ normalized \ to \ the \ peak \ irradiance. \ The \ vertical,}$ \n' \
                         ' $\mathit{solid \ black \ line \ at \ separation \ ' + str(radius_fpm) + \
                         '\ λ/D \ marks \ the \ radius \ of \ the \ FPM \ occulting}$ \n $\mathit{spot. \ The \ ' \
                         'vertical, \ red \ lines \ at \ ' + str(float(iwa)) + ' \ and \ ' + str(float(owa)) + \
                         ' \ λ/D \ respectively \ indicate \ the}$ \n $\mathit{radii \ of \ the \ inner \ and \ ' \
                         'outermost \ constraints \ applied \ during \ the \ apodizer}$\n $\mathit{optimization.}$'
    else:
        cmap = mpl.cm.get_cmap('Spectral_r')
        cols = cmap(np.linspace(0, 1, 11))
        norm = mpl.colors.Normalize(vmin=1 - bandwidth / 2, vmax=1 + bandwidth / 2)

        for contrast_at_wl, col in zip(imgs, cols):
            r, y, _, _ = radial_profile(contrast_at_wl, 0.1)
            plt.plot(r, y, c=col)

        plt.colorbar(mpl.cm.ScalarMappable(norm, cmap)).set_label('Relative wavelength')

        plt.title('\n Normalized Irradiance (Radial Average)')
        caption_radial = '\n \n Radial intensity profile for the broadband APLC design at 11 simulated wavelengths' \
                         'centered \n around $\lambda_0/D$ and equally spatially sampled over the ' \
                         + str(
            bandwidth * 100) + '% bandpass. The black curve shows the average \n intensity across the 11 wavelength samples. The dashed red vertical lines delimit' \
                               'the high-contrast dark\n zone (between ' + str(iwa) + ' and ' + str(
            owa) + ' $\lambda_0/D$).' \
                   ' The blue dotted line delimits the FPM radius, set to ' + str(round(radius_fpm, 2)) + \
                         ' $\lambda_0/D$. '

    plt.ylabel('Normalized intensity in log scale')
    plt.xlabel(r'Angular separation ($\lambda_0/D$)' + caption)
    plt.title('')
    plt.tight_layout()

    plt.axvline(iwa, color=colors.red, linestyle='-.')
    plt.axvline(owa, color=colors.red, linestyle='--')
    plt.axvline(radius_fpm, color='k')
    plt.axhline(10 ** (-contrast), xmin=0, xmax=owa * 1.2, linewidth=1, color='k', linestyle='--')
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=3, fancybox=True, shadow=True)
    plt.yscale('log')
    plt.xlim(0, owa * 1.2)
    plt.ylim(5e-12, 2e-5)
    plt.ylabel('Normalized irradiance')
    plt.xlabel(r'Angular separation ($\lambda_0/D$)' + caption_radial)
    plt.tight_layout()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    if not monochromatic:
        wavelengths = np.linspace(0.5, 1.5, 151)

        contrasts = []
        for wl in wavelengths:
            psf = prop(coro(Wavefront(apodizer * pupil, wl))).power
            psf_ref = prop(Wavefront(apodizer * pupil * lyot_stop, wl)).power

            r, y, _, _ = radial_profile(psf / psf_ref.max(), 0.1)
            contrasts.append(y)

        grid = CartesianGrid(SeparatedCoords([r, wavelengths]))
        contrasts = Field(np.array(contrasts).ravel(), grid)

        fig, (ax, cax) = plt.subplots(ncols=2, figsize=(10, 8), dpi=100,
                                      gridspec_kw={"width_ratios": [1, 0.05]})
        imshow_field(np.log10(contrasts), aspect='auto', vmin=-contrast - 0.5, vmax=-contrast + 3.5, cmap='inferno',
                     interpolation='bilinear',
                     ax=ax)
        cs = contour_field(np.log10(contrasts), levels=[-7, -6, -5, -4], colors='w', linestyles='-', ax=ax)
        ax.axhline(1 - bandwidth / 2, c='w', ls=':')
        ax.axhline(1 + bandwidth / 2, c='w', ls=':')
        ax.axvline(radius_fpm, c='w', ls=':')
        ax.clabel(cs, fmt='1e%d')
        plt.colorbar(cax=cax).set_label('log10(normalized irradiance)')
        ax.set_xlabel('Angular separation [$\lambda_0/D$]')
        ax.set_ylabel('Wavelength [$\lambda/\lambda_0$]')
        ax.set_xlim(0, 20)

        def ang_sep_to_arcsec(x):
            return x * (1.593e-6 / 8) * 206265 * 1000

        def arcsec_to_ang_sep(x):
            return x / 206265 / (1.593e-6 / 8) / 1000

        ax2 = ax.secondary_xaxis('top', functions=(ang_sep_to_arcsec, arcsec_to_ang_sep))
        ax2.set_xlabel('Angular separation [mas]')

        def frac_wl_to_wl(x):
            return x * 1.593

        def wl_to_frac_wl(x):
            return x / 1.593

        ax3 = ax.secondary_yaxis('right', functions=(frac_wl_to_wl, wl_to_frac_wl))
        ax3.set_ylabel('Wavelength [um]')

        if pdf is not None:
            pdf.savefig()
            plt.close()
        else:
            plt.show()

    return {'normalized_irradiance_image': img / img_ref.max(),
            'normalized_irradiance_radial': (r, profile, std_profile, n_profile)}


def analyze_max_throughput(solution_filename, pdf=None):
    """Calculate the the sum of the PSF in a fixed circular photometric aperture of radius 0.7 λ0/D.

    Parameters
    ----------
    solution_filename:
        The location of the apodizer solution file.
    pdf: bool
        Whether to save to file or display on screen.

    Returns
    -------
    The Band-averaged r=.7 lam/D throughput.
    """
    Pupil, Apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)

    fname_pup = parameters['pupil']['filename']
    fname_ls = parameters['lyot_stop']['filename']

    if not os.path.isabs(fname_pup):
        fname_pup = os.path.join(file_organization['input_files_dir'], fname_pup)
    if not os.path.isabs(fname_ls):
        fname_ls = os.path.join(file_organization['input_files_dir'], fname_ls)

    fname_apod = solution_filename

    TelAp = fits.getdata(fname_pup)
    LS = fits.getdata(fname_ls)
    A = fits.getdata(fname_apod)

    maximum_integrated_throughput = ((TelAp * LS * A).sum() / (TelAp * LS).sum()) ** 2

    # Account for the ratio of the diameter of the square enclosing the aperture to the circumscribed circle
    D = 0.982  # for gpi
    N = parameters['pupil']['N']
    rho_out = parameters['image']['owa']
    fp2res = parameters['focal_plane_mask']['num_pix']
    bw = parameters['image']['bandwidth']
    Nlam = parameters['image']['num_wavelengths']

    dx = (D / 2) / N
    dy = dx
    xs = np.matrix(np.linspace(-N + 0.5, N - 0.5, N) * dx)
    ys = xs.copy()
    M_fp2 = int(np.ceil(rho_out * fp2res))  # rounded to next intiger (owa * fpm_num_pix)
    dxi = 1. / fp2res
    xis = np.matrix(np.linspace(-M_fp2 + 0.5, M_fp2 - 0.5, 2 * M_fp2) * dxi)
    etas = xis.copy()
    wrs = np.linspace(1.00 - bw / 2, 1.00 + bw / 2, Nlam)
    XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
    YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
    RRs = np.sqrt(XXs ** 2 + YYs ** 2)
    p7ap_ind = np.less_equal(RRs, 0.7)

    intens_D_0_polychrom = np.zeros((Nlam, 2 * M_fp2, 2 * M_fp2))
    intens_D_0_peak_polychrom = np.zeros((Nlam, 1))
    intens_TelAp_polychrom = np.zeros((Nlam, 2 * M_fp2, 2 * M_fp2))
    intens_TelAp_peak_polychrom = np.zeros((Nlam, 1))

    for wi, wr in enumerate(wrs):
        Psi_D_0 = dx * dy / wr * np.dot(
            np.dot(np.exp(-1j * 2 * np.pi / wr * np.dot(xis.T, xs)), TelAp * A * LS[::-1, ::-1]),
            np.exp(-1j * 2 * np.pi / wr * np.dot(xs.T, xis)))
        intens_D_0_polychrom[wi] = np.power(np.absolute(Psi_D_0), 2)

    intens_D_0 = np.mean(intens_D_0_polychrom, axis=0)

    p7ap_sum_APLC = np.sum(intens_D_0[p7ap_ind]) * dxi * dxi
    p7ap_circ_thrupt = p7ap_sum_APLC / (np.pi / 4)

    print("\n max integrated throughput: {}".format(maximum_integrated_throughput))

    fits.setval(solution_filename, 'MAX_THR', value=maximum_integrated_throughput,
                comment='Maximum integrated throughput')

    return {'P7APTH_throughput': p7ap_circ_thrupt}


def analyze_offaxis_throughput(solution_filename, pdf=None):
    pass


def analyze_summary(solution_filename, pdf=None):
    """Display an example of light propagation through the APLC design, showing the light in each critical plane.

    Parameters
    ----------
    solution_filename:
        The location of the apodizer solution file.
    pdf: bool
        Whether to save the figure to file or display on screen.

    Returns
    -------
    An empty dictionary.
    """
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    lyot_stop = lyot_stops[0]
    owa = parameters['image']['owa']
    bandwidth = parameters['image']['bandwidth']
    radius_fpm = parameters['focal_plane_mask']['radius']
    contrast = parameters['image']['contrast']
    num_wavelengths = parameters['image']['num_wavelengths']

    if num_wavelengths == 1:
        monochromatic = True
    else:
        monochromatic = False
    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
    coro_without_lyot = LyotCoronagraph(pupil.grid, focal_plane_mask)
    focal_grid = make_focal_grid(8, owa * 1.2)
    prop = FraunhoferPropagator(pupil.grid, focal_grid)
    wavelengths = np.linspace(-bandwidth / 2, bandwidth / 2, 11) + 1

    focal_plane_mask_large = 1 - circular_aperture(2 * radius_fpm)(focal_grid)

    if monochromatic:
        wf = Wavefront(pupil * apodizer)
        img = prop(coro(wf)).intensity
        img_foc = prop(wf).intensity
        img_ref = prop(Wavefront(apodizer * pupil * lyot_stop)).intensity
        lyot = coro_without_lyot(wf).intensity

    else:
        img = 0
        img_ref = 0
        img_foc = 0
        lyot = 0

        for wl in wavelengths:
            wf = Wavefront(pupil * apodizer, wl)
            img += prop(coro(wf)).intensity
            img_foc += prop(wf).intensity
            img_ref += prop(Wavefront(pupil * apodizer * lyot_stop)).intensity
            lyot += coro_without_lyot(wf).intensity

        img /= 11
        img_ref /= 11
        img_foc /= 11
        lyot /= 11

    font = {'family': 'DejaVu Sans', 'weight': 'medium', 'size': 9}
    matplotlib.rc('font', **font)
    plt.figure(dpi=100)
    # apodizer and telescope aperture
    plt.suptitle('Analysis Summary \n \n ', fontsize=12, fontweight='bold', style='italic')
    plt.subplot(2, 3, 1)
    imshow_field(pupil * apodizer, cmap='Greys_r')
    plt.title('\n Apodizer & \n Telescope Aperture')
    plt.axis('off')

    # image plane
    plt.subplot(2, 3, 2)
    im = imshow_field(np.log10(img_foc / img_foc.max()), vmin=-5, vmax=0, cmap='inferno')
    plt.title('\n Image plane')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # image plane masked by focal plane mask
    plt.subplot(2, 3, 3)
    im = imshow_field(np.log10(img_foc / img_foc.max() * (1e-20 + focal_plane_mask_large)), vmin=-5, vmax=0,
                      cmap='inferno')
    plt.title('\n Image plane \n w/FPM')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # lyot plane
    plt.subplot(2, 3, 4)
    im = imshow_field(np.log10(lyot / lyot.max()), vmin=-3, vmax=0, cmap='inferno')
    plt.title('Lyot plane')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # lyot plane masked by lyot stop
    plt.subplot(2, 3, 5)
    im = imshow_field(np.log10(lyot / lyot.max() * (1e-20 + lyot_stop)), vmin=-3, vmax=0, cmap='inferno')
    plt.title('Lyot plane \n w/lyot stop')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # final image plane
    plt.subplot(2, 3, 6)
    im = imshow_field(np.log10(img / img_ref.max()), vmin=-contrast - 0.5, vmax=-contrast + 3.5, cmap='inferno')
    plt.title('Final image plane')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')
    plt.tight_layout()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {}


def analyze_lyot_robustness(solution_filename, pdf=None):
    """Plot a grid of the coronagraphic image for different translations of the Lyot stop mask.

    Parameters
    ----------
    solution_filename:
        The location of the apodizer solution file.
    pdf: bool
        Whether to save the table to file or display on screen.

    Returns
    -------
    An empty dictionary.
    """
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    owa = parameters['image']['owa']
    bandwidth = parameters['image']['bandwidth']
    contrast = parameters['image']['contrast']
    num_wavelengths = parameters['image']['num_wavelengths']
    ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']

    focal_grid = make_focal_grid(8, owa * 1.2)
    prop = FraunhoferPropagator(pupil.grid, focal_grid)
    coro_no_lyot = LyotCoronagraph(pupil.grid, focal_plane_mask)
    lyot_stop = lyot_stops[0]

    def get_image(coronagraphic=True, lyot_dx=0, lyot_dy=0):
        lyot_stop_shifted = np.roll(lyot_stop.shaped, (lyot_dy, lyot_dx), (1, 0)).ravel()

        img = 0
        for wl in np.linspace(1 - bandwidth / 2, 1 + bandwidth / 2, 11):
            wf = Wavefront(apodizer * pupil, wl)
            if coronagraphic:
                wf = coro_no_lyot(wf)
            wf.electric_field *= lyot_stop_shifted
            img += prop(wf).power
        return img

    if num_wavelengths == 1:
        monochromatic = True
        label = 'monochromatic'
    else:
        monochromatic = False
        label = 'broadband'

    # dxs = np.array([-8, -6, -4, -2, 0, 2, 4, 6, 8])
    dxs = np.array(range(-ls_alignment_tolerance - 1, +ls_alignment_tolerance + 1 + 1, 1))

    dither_grid = CartesianGrid(SeparatedCoords((dxs, dxs)))

    plt.figure(figsize=(8, 8), dpi=100)

    for i, (dx, dy) in enumerate(dither_grid.points):

        if monochromatic:
            lyot_stop_shifted = np.roll(lyot_stop.shaped, (dy, dx), (1, 0)).ravel()

            wf = Wavefront(apodizer * pupil)
            img_wf = coro_no_lyot(wf)
            img_wf.electric_field *= lyot_stop_shifted
            wf.electric_field *= lyot_stop_shifted

            img = prop(img_wf).power
            img_ref = prop(wf).power

        else:
            img = get_image(lyot_dx=dx, lyot_dy=dy)
            img_ref = get_image(coronagraphic=False, lyot_dx=dx, lyot_dy=dy)

        x = i % len(dxs)
        y = i // len(dxs)

        plt.subplot(len(dxs), len(dxs), x + (len(dxs) - y - 1) * len(dxs) + 1)

        imshow_field(np.log10(img / img_ref.max()), vmin=-contrast - 0.5, vmax=-contrast + 3.5, cmap='inferno')

        frame1 = plt.gca()
        frame1.axes.xaxis.set_ticklabels([])
        frame1.axes.yaxis.set_ticklabels([])

        # frame1.spines['bottom'].set_color('w')
        # frame1.spines['top'].set_color('w')
        # frame1.spines['right'].set_color('w')
        # frame1.spines['left'].set_color('w')
        # frame1.tick_params(axis='x', colors='w')
        # frame1.tick_params(axis='y', colors='w')

        if dx == dxs.min():
            plt.ylabel('%.1f%%' % (dy / pupil.grid.shape[1] * 100), fontsize=8)

        if dy == dxs.min():
            plt.xlabel('%.1f%%' % (dx / pupil.grid.shape[0] * 100), fontsize=8)

    plt.subplots_adjust(top=0.98, bottom=0.02, left=0.02, right=0.98, wspace=0, hspace=0)

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    dxs_tol = np.array(range(-ls_alignment_tolerance, +ls_alignment_tolerance + 1, 1))
    dither_grid_tol = CartesianGrid(SeparatedCoords((dxs_tol, dxs_tol)))

    imgs = 0
    img_refs = 0

    for i, (dx, dy) in enumerate(dither_grid_tol.points):

        if monochromatic:
            lyot_stop_shifted = np.roll(lyot_stop.shaped, (dy, dx), (1, 0)).ravel()

            wf = Wavefront(apodizer * pupil)
            img_wf = coro_no_lyot(wf)
            img_wf.electric_field *= lyot_stop_shifted
            img = prop(img_wf).power

            wf.electric_field *= lyot_stop_shifted
            img_ref = prop(wf).power

        else:
            img = get_image(lyot_dx=dx, lyot_dy=dy)
            img_ref = get_image(coronagraphic=False, lyot_dx=dx, lyot_dy=dy)

        imgs += img
        img_refs += img_ref

    imgs /= len(dither_grid_tol)
    img_refs /= len(dither_grid_tol)

    plt.figure(figsize=(6, 4), dpi=100)
    imshow_field(np.log10(imgs / img_refs.max()), vmin=-contrast - 0.5, vmax=-contrast + 3.5, cmap='inferno')
    plt.colorbar()
    plt.title('Average ' + label + ' normalized irradiance \n within ' + str(
        ls_alignment_tolerance) + ' pixels of Lyot stop misalignment.')

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {}


def analyze_throughput(solution_filename, pdf=None):
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    bandwidth = parameters['image']['bandwidth']
    fpm_radius = parameters['focal_plane_mask']['radius']
    num_wavelengths = parameters['image']['bandwidth']

    if num_wavelengths == 1:
        monochromatic = True
    else:
        monochromatic = False

    lyot_stop = lyot_stops[0]
    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)

    focal_grid_hires = make_focal_grid(64, 2)
    prop_hires = FraunhoferPropagator(pupil.grid, focal_grid_hires)

    wf_pup = Wavefront(pupil)

    mask = circular_aperture(1.4)(focal_grid_hires)

    tips = np.linspace(0, 10, 51)
    core_throughputs = []

    def get_image_offaxis(tip):
        img = 0

        for wl in np.linspace(1 - bandwidth / 2, 1 + bandwidth / 2, 11):
            wf = Wavefront(apodizer * pupil * np.exp(2j * np.pi * tip / wl * pupil.grid.x), wl)
            wf.total_power /= wf_pup.total_power

            wf_post_lyot = coro(wf)
            wf_post_lyot.electric_field *= np.exp(-2j * np.pi * tip / wl * pupil.grid.x)
            img += prop_hires(wf_post_lyot).power

        return img / 11

    imgs = []
    for tt in tips:

        if monochromatic:
            wf = Wavefront(apodizer * pupil * np.exp(2j * np.pi * tip / wl * pupil.grid.x))
            wf.total_power /= wf_pup.total_power

            wf_post_lyot = coro(wf)
            wf_post_lyot.electric_field *= np.exp(-2j * np.pi * tip / wl * pupil.grid.x)
            img = prop_hires(wf_post_lyot).power

        else:
            img = get_image_offaxis(tt)

        imgs.append(img)

        core_throughput = (img * mask).sum()
        core_throughputs.append(core_throughput)

    # Calculate limits monochromatically
    wf = Wavefront(pupil)
    wf.total_power = 1
    img = prop_hires(wf).power
    core_throughput_pupil = (img * mask).sum()

    wf.electric_field *= lyot_stop
    img = prop_hires(wf).power
    core_throughput_lyot = (img * mask).sum()

    wf = Wavefront(apodizer * pupil * lyot_stop)
    wf.total_power /= wf_pup.total_power
    img = prop_hires(wf).power
    core_throughput_max = (img * mask).sum()

    core_throughputs = np.array(core_throughputs)

    half_throughput = 0.5 * core_throughput_max
    iwa_throughput = np.interp(half_throughput, core_throughputs, tips)

    plt.figure(figsize=(8, 6), dpi=100)
    plt.axhline(core_throughput_pupil, label='Limit given by the telescope pupil', c='k', ls='--')
    plt.axhline(core_throughput_lyot, label='Limit given by the Lyot stop', c=colors.orange, ls='--')
    plt.plot(tips, core_throughputs, label='Core throughput of the coronagraph', c=colors.blue, lw=2)
    plt.axhline(core_throughput_max, c=colors.blue, ls='--')
    plt.plot([0, iwa_throughput, iwa_throughput], [half_throughput, half_throughput, 0], ls=':', c=colors.blue)
    plt.plot([iwa_throughput], [half_throughput], 'o', c=colors.blue)
    plt.axvline(fpm_radius, ls='--', c='k')
    plt.xlabel('Angular separation [$\lambda/D$]')
    plt.ylabel('Core throughput')
    plt.xlim(0, 10)
    plt.ylim(0, core_throughput_pupil * 1.05)
    plt.legend(loc='lower right')
    plt.grid(ls=':', c='0.7')

    table_data = [['Pupil core throughput: ', str(core_throughput_pupil)],
                  ['Lyot stop core throughput: ', str(core_throughput_lyot)],
                  ['Maximum core throughput: ', str(core_throughput_max)],
                  ['Maximum core throughput w.r.t. pupil core throughput: ',
                   str(core_throughput_max / core_throughput_pupil)],
                  ['Maximum core throughput w.r.t. Lyot stop core throughput: ',
                   str(core_throughput_max / core_throughput_lyot)],
                  ['Inner working angle: ', str(iwa_throughput) + r' $\lambda_0/D$']]

    plt.table(cellText=table_data, colWidths=[0.7, 0.3], cellLoc='right', edges='open',
              bbox=[0.0, -0.45, 1, .3])
    plt.subplots_adjust(bottom=0.3)

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    seps = np.arange(200) * 0.25
    index = np.argmin(np.abs(seps - iwa_throughput))
    seps = seps[index - 4:index + 4]

    img_max = get_image_offaxis(10).max()

    plt.figure(dpi=100)
    for i, sep in enumerate(seps):
        img = get_image_offaxis(sep)
        img.grid = img.grid.shifted([sep, 0])

        plt.subplot(2, 4, i + 1)
        plt.title('%0.2f $\lambda_0/D$' % sep)
        imshow_field(img / img_max, cmap='inferno', vmax=1)

        circle = plt.Circle((0, 0), fpm_radius, edgecolor='w', linestyle='--', fill=False)
        plt.gca().add_patch(circle)

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return core_throughput_max / core_throughput_pupil


def analyze_tt_jitter(solution_filename, pdf=None):
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    iwa = parameters['image']['iwa']
    owa = parameters['image']['owa']
    bandwidth = parameters['image']['bandwidth']
    fpm_radius = parameters['focal_plane_mask']['radius']
    img_contrast = parameters['image']['contrast']
    pupil_grid = pupil.grid
    lyot_stop = lyot_stops[0]

    focal_grid = make_focal_grid(8, owa * 1.2)

    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
    prop = FraunhoferPropagator(pupil.grid, focal_grid)

    imgs_tt = []

    tt_rmss = [0.01, 0.03, 0.1, 0.3]
    N = 100

    plt.figure(figsize=(10, 8), dpi=100)
    for i, tt_rms in enumerate(tt_rmss):
        tips = np.random.randn(N) * tt_rms / np.sqrt(2)
        tilts = np.random.randn(N) * tt_rms / np.sqrt(2)
        wls = np.random.uniform(1 - bandwidth / 2, 1 + bandwidth / 2, N)

        img = 0
        img_ref = 0

        for j in range(N):
            wf = Wavefront(apodizer * pupil, wls[j])
            wf.electric_field *= np.exp(2j * np.pi * (tips[j] * pupil_grid.x + tilts[j] * pupil_grid.y))

            wf_post_lyot = coro(wf)
            img += prop(wf_post_lyot).power

            wf_nocoro = wf.copy()
            wf_nocoro.electric_field *= lyot_stop
            img_ref += prop(wf_nocoro).power

        img /= N
        img_ref /= N

        contrast = img / img_ref.max()
        imgs_tt.append(contrast)

        plt.subplot(2, 2, i + 1)
        plt.title('$\sigma$_rms = %.2f $\lambda/D$' % tt_rms)
        imshow_field(np.log10(contrast), vmin=-img_contrast - 0.5, vmax=-img_contrast + 3.5, cmap='inferno')
        plt.suptitle(
            '\n\n Broadband normalized irradiance for four representative levels of residual pointing jitter.',
            fontsize=10)
        plt.colorbar()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    plt.figure(figsize=(8, 6), dpi=100)
    for i, img in enumerate(imgs_tt):
        r, y, _, _ = radial_profile(img, 0.1)

        plt.plot(r, y, label='$\sigma$_rms = %.2f $\lambda/D$' % tt_rmss[i])

    plt.yscale('log')
    plt.xlim(0, owa * 1.2)
    plt.ylim(5e-12, 2e-3)
    plt.legend(loc='upper right')

    caption = '\n\n Azimuthally averaged raw contrast for four representative levels of rms residual pointing jitter.'

    plt.xlabel('Angular separation [$\lambda/D$]' + caption)
    plt.ylabel('Normalized irradiance')
    plt.grid(ls=':', c='0.7')
    plt.axvline(fpm_radius, ls='--', c='k')
    plt.axvline(iwa, color=colors.red, linestyle='-.')
    plt.axvline(owa, color=colors.red, linestyle='--')
    plt.axhline(10 ** (-7), xmin=0, xmax=owa * 1.2, linewidth=1, color='k', linestyle='--')
    plt.axvline()
    plt.tight_layout()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

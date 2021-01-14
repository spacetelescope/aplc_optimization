import matplotlib as mpl
import numpy as np
from astropy.io import fits
from hcipy import *

#mpl.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
import asdf
import os
import re


def create_coronagraph(solution_filename):
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


def analyze_pdf_summary(solution_filename, pdf=None):
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)

    #FPM
    fpm_radius = parameters['focal_plane_mask']['radius']
    fpm_num_pix = parameters['focal_plane_mask']['num_pix']
    fpm_grayscale = parameters['focal_plane_mask']['grayscale']
    fpm_field_radius = parameters['focal_plane_mask']['field_stop_radius']
    #Image
    img_contrast = 10 ** (-parameters['image']['contrast'])
    img_iwa = parameters['image']['iwa']
    img_owa = parameters['image']['owa']
    img_num_wavelengths = parameters['image']['num_wavelengths']
    img_bandwidth = 100 ** parameters['image']['bandwidth']
    img_resolution = parameters['image']['resolution']
    #Lyot stop
    ls_fname = parameters['lyot_stop']['filename']
    ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']
    ls_num_stops = parameters['lyot_stop']['num_lyot_stops']
    #pupil
    pup_fname = parameters['pupil']['filename']
    N = parameters['pupil']['N']


    # Extract telescope name, gap padding and oversampling info from the pupil filename
    if "LUVOIR" in pup_fname:
        telescope = "LUVOIR"
    elif "HiCAT" in pup_fname:
        telescope = "HiCAT"
    elif "GPI" in pup_fname:
        telescope = "GPI"


    regex = re.compile(r'\d+')
    pup_vals = regex.findall(pup_fname)
    if telescope is "LUVOIR":
        seg_gap_pad, oversamp = int(pup_vals[0]), int(pup_vals[1])
    elif telescope is "HiCAT":
        seg_gap_pad = ""
        oversamp = int(pup_vals[0])
    elif telescope is "GPI":
        seg_gap_pad = ""
        oversamp = ""

    # Extract lyot stop inner and outer diameter from lyot stop filename
    ls_vals = regex.findall(ls_fname)
    ls_id, ls_od = int(ls_vals[0]) / 1000, int(ls_vals[1]) / 1000


    if fpm_grayscale is True:
        gs = ' (grayscale)'
    else:
        gs = ''

    fig = plt.figure(dpi=160)

    col_labels = ['$\mathbf{APLC \ Analysis \ Summary}$' + ' \n \n', ' ']

    ax = fig.add_subplot(1, 1, 1)

    table_data = [
        ["APLC design", "$\mathtt{" + str(img_bandwidth) + "\%}$"],
        ["nPup", "$\mathtt{" + str(N) + " \ x \ " + str(N) + " \ pixels}$"],
        ["Gap padding (multiplicative)", "$\mathtt{" + str(seg_gap_pad) + "}$"],
        ["Oversampling (grey levels)", "$\mathtt{" + str(oversamp) + "}$"],
        ["Telescope", "$\mathtt{" + str(telescope) + "}$"],
        ["Lyot stop inner diamater (% of inscribed circle)", "$\mathtt{" + str(ls_id) + "}$"],
        ["Lyot stop outer diameter (% of inscribed circle)", "$\mathtt{" + str(ls_od) + "}$"],
        ["Bandpass", "$\mathtt{" + str(img_bandwidth) + "\%}$"],
        ["# wavelengths", "$\mathtt{" + str(img_num_wavelengths) + "}$"],
        ["FPM radius" + str(gs), "$\mathtt{" + str(fpm_radius) + ' \ \lambda/D}$'],
        ["nFPM", "$\mathtt{" + str(fpm_num_pix) + ' \ pixels}$'],
        ["IWA " + u"\u2014" + " OWA",
         "$\mathtt{" + str(float(img_iwa)) + "\u2014" + str(float(img_owa)) + ' \ \lambda/D}$'],
        ["", ""],
        ["$\mathit{Optimizer \ called \ with \ the \ following \ parameters:}$", ""],
        [r"       $\triangleright \ \ Pupil \ file:$         " + str(pup_fname), ""],
        [r"       $\triangleright \ \ Lyot \ stop \ file:$ " + str(ls_fname), ""]

    ]
    table = ax.table(cellText=table_data, colLabels=col_labels, colWidths=[0.8, 0.4], cellLoc='left', edges='open',
                     loc='center')
    table.set_fontsize(14)
    table.scale(1,1)
    plt.axis('off')



    #    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.set_title('Optimizer called with the following parameters:', fontsize=10)
#    txt = ax.text(0.01, 1.0,
#                  '\n Focal Plane Mask: \n'
#                    '     - Field stop radius: '+ str(fpm_field_radius) +'\n'
#                    '     - Grayscale: '+str(fpm_grayscale)+'\n'
#                    '     - Number of pixels: ' + str(fpm_num_pix) + '\n'
#                    '     - Radius'+gs+': '+ str(fpm_radius) +'\n'
#                  'Image: \n'
#                    '     - Bandpass: '+ str(img_bandwidth) +'% \n'
#                    '     - Contrast: '+ str(img_contrast) +'\n'
#                    '     - IWA - OWA: '+ str(img_iwa)+' - '+str(img_owa)+' '+r'$\lambda$'+'/D \n'
#                    '     - # wavelengths: '+ str(img_num_wavelengths) +'\n'
#                    '     - Resolution: ' + str(img_resolution) +'\n'
#                  'Lyot stop:\n'
#                    '     - Alignment tolerance: '+ str(ls_alignment_tolerance) +'\n'
#                    '     - Filename: '+ str(ls_fname) +'\n'
#                    '     - Number of Lyot stops: '+ str(ls_num_stops) +'\n'
#                  'Pupil: \n'
#                    '     - nPup: '+str(N)+' x '+ str(N)+'\n'
#                    '     - Filename: '+ str(pup_fname), verticalalignment='top', horizontalalignment='left', fontsize=9, clip_on=True)
#    ax.axis('off')

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {}


def analyze_contrast_monochromatic(solution_filename, pdf=None):
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    iwa = parameters['image']['iwa']
    owa = parameters['image']['owa']
    contrast = parameters['image']['contrast']
    radius_fpm = parameters['focal_plane_mask']['radius']

    lyot_stop = lyot_stops[0]

    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
    focal_grid = make_focal_grid(8, owa * 1.2)  # make_focal_grid(q, fov) - grid for a focal plane
    prop = FraunhoferPropagator(pupil.grid, focal_grid)

    wf = Wavefront(pupil * apodizer)
    img = prop(coro(wf)).intensity
    img_ref = prop(Wavefront(apodizer * lyot_stop)).intensity

    fig, ax = plt.subplots()

    plt.title('Monochromatic Normalized Irradiance')
    imshow_field(np.log10(img / max(img_ref)), vmin=-contrast - 1, vmax=-contrast + 4, cmap='inferno')
    cbar = plt.colorbar()
    cbar.set_label('log irradiance')

    # For the figure caption
    caption = '\n \n $\mathit{\mathbf{Figure \ 1:} \ Monochromatic \ on-axis \ PSF \ in \ log \ irradiance,}$ \n ' \
              '$\mathit{normalized \ to \ the \ peak \ irradiance \ value.}$'

    plt.xlabel('Angular coordinate ($\lambda_0/D$)'+caption, fontsize=10)
    plt.ylabel('Angular coordinate ($\lambda_0/D$)', fontsize=10)
    plt.tick_params(axis='both', labelsize=8)
    plt.tight_layout()

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    r, profile, std_profile, n_profile = radial_profile(img / max(img_ref), 0.2)

    plt.figure(figsize=(6, 6))

    plt.title('\n Monochromatic Normalized Irradiance (Radial Average)')
    plt.plot(r, profile)
    iwa_line = plt.axvline(iwa, color=colors.red, linestyle = '-.')
    owa_line = plt.axvline(owa, color=colors.red, linestyle = '--')
    radius_line = plt.axvline(radius_fpm, color='k')
    plt.axhline(10 ** (-contrast), xmin=0, xmax=owa * 1.2, linewidth=1, color='k', linestyle='--')
    plt.legend([iwa_line, owa_line, radius_line],
               [r'$\rho_o$ = '+str(float(iwa))+r' $\lambda_0/D$', r'$\rho_i$ = '+str(float(owa))+r' $\lambda_0/D$',
                'FPM radius = '+str(float(radius_fpm))+r' $\lambda_0/D$'])


    caption_radial ='\n \n \n Figure 2: ' \
                    '$\mathit{Monochromatic \ on-axis \ PSF \ azimuthally \ averaged \ over \ angular}$ \n $\mathit{seperations \ }$' \
                    +str(round(min(r),4))+'-'+str(round(max(r),4))+'$\mathit{ \ λ/D, \ normalized \ to \ the \ peak \ irradiance. \ ' \
                    'The \ vertical,}$ \n $\mathit{solid \ black \ line \ at \ separation \ '+str(radius_fpm)+ \
                    '\ λ/D \ marks \ the \ radius \ of \ the \ FPM \ occulting}$ \n $\mathit{spot. \ The \ vertical, \ red ' \
                    ' lines \ at \ '+str(float(iwa))+' \ and \ '+str(float(owa))+ \
                    ' \ λ/D \ respectively \ indicate \ the}$ \n $\mathit{radii \ of \ the \ inner \ and \ outermost \ ' \
                    'constraints \ applied \ during \ the \ apodizer}$\n $\mathit{optimization.}$'

    plt.yscale('log')
    plt.xlim(0, owa * 1.2)
    plt.ylim(5e-12, 2e-5)
    plt.ylabel('Normalized irradiance')
    plt.xlabel(r'Angular separation ($\lambda_0/D$)'+caption_radial)
    plt.tight_layout()


    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {'normalized_irradiance_image': img / max(img_ref),
            'normalized_irradiance_radial': (r, profile, std_profile, n_profile)}


def analyze_max_throughput(solution_filename, pdf=None):
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
    D = 0.982
    N = parameters['pupil']['N']
    rho_out = parameters['image']['owa']
    fp2res = parameters['focal_plane_mask']['num_pix']
    bw = parameters['image']['bandwidth']
    Nlam = parameters['image']['num_wavelengths']

    dx = (D / 2) / N
    dy = dx
    xs = np.matrix(np.linspace(-N + 0.5, N - 0.5, N) * dx)
    ys = xs.copy()
    M_fp2 = int(np.ceil(rho_out * fp2res))
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

    print(p7ap_circ_thrupt)

    fits.setval(solution_filename, 'P7APTH', value=p7ap_circ_thrupt, comment='Band-averaged r=.7 lam/D throughput')

    return {'P7APTH_throughput': p7ap_circ_thrupt}


def analyze_offaxis_throughput(solution_filename, pdf=None):
    pass


def analyze_summary(solution_filename, pdf=None):
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
    lyot_stop = lyot_stops[0]
    owa = parameters['image']['owa']
    bandwidth = parameters['image']['bandwidth']
    radius_fpm = parameters['focal_plane_mask']['radius']
    contrast = parameters['image']['contrast']

    coro = LyotCoronagraph(pupil.grid, focal_plane_mask, lyot_stop)
    coro_without_lyot = LyotCoronagraph(pupil.grid, focal_plane_mask)
    focal_grid = make_focal_grid(8, owa * 1.2)
    prop = FraunhoferPropagator(pupil.grid, focal_grid)
    wavelengths = np.linspace(-bandwidth / 2, bandwidth / 2, 11) + 1

    focal_plane_mask_large = 1 - circular_aperture(2 * radius_fpm)(focal_grid)

    img = 0
    img_ref = 0
    img_foc = 0
    lyot = 0

    for wl in wavelengths:
        wf = Wavefront(pupil * apodizer, wl)
        img += prop(coro(wf)).intensity
        img_foc += prop(wf).intensity
        img_ref += prop(Wavefront(pupil * lyot_stop)).intensity
        lyot += coro_without_lyot(wf).intensity

    font = {'family': 'normal', 'weight': 'medium', 'size': 9}
    matplotlib.rc('font', **font)

    # apodizer and telescope aperture
    plt.suptitle('Analysis Summary \n \n ', fontsize=12, fontweight='bold', style='italic')
    plt.subplot(2, 3, 1)
    imshow_field(pupil * apodizer, cmap='Greys_r')
    plt.title('\n Apodizer & \n Telescope Aperture')
    plt.axis('off')

    # image plane
    plt.subplot(2, 3, 2)
    im = imshow_field(np.log10(img_foc / max(img_foc)), vmin=-5, vmax=0, cmap='inferno')
    plt.title('\n Image plane')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # image plane masked by focal plane mask
    plt.subplot(2, 3, 3)
    im = imshow_field(np.log10(img_foc / max(img_foc) * (1e-20 + focal_plane_mask_large)), vmin=-5, vmax=0,
                      cmap='inferno')
    plt.title('\n Image plane \n w/FPM')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # lyot plane
    plt.subplot(2, 3, 4)
    im = imshow_field(np.log10(lyot / max(lyot)), vmin=-3, vmax=0, cmap='inferno')
    plt.title('Lyot plane')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # lyot plane masked by lyot stop
    plt.subplot(2, 3, 5)
    im = imshow_field(np.log10(lyot / max(lyot) * (1e-20 + lyot_stop)), vmin=-3, vmax=0, cmap='inferno')
    plt.title('Lyot plane \n w/lyot stop')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    plt.axis('off')

    # final image plane
    plt.subplot(2, 3, 6)
    im = imshow_field(np.log10(img / max(img_ref)), vmin=-contrast - 1, vmax=-contrast + 4, cmap='inferno')
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
    pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)

    num_pix = parameters['pupil']['N']  # px
    iwa = parameters['image']['iwa']
    owa = parameters['image']['owa']
    bandwidth = parameters['image']['bandwidth']
    radius_fpm = parameters['focal_plane_mask']['radius']
    contrast = parameters['image']['contrast']

    ls_alignment_tolerance = parameters['lyot_stop']['alignment_tolerance']

    wavelengths = np.linspace(-bandwidth / 2, bandwidth / 2, 11) + 1

    focal_grid = make_focal_grid(8, owa * 1.2)
    prop = FraunhoferPropagator(pupil.grid, focal_grid)

    lyot_stop = lyot_stops[0]

    #dxs = np.array([-8, -6, -4, -2, 0, 2, 4, 6, 8])
    dxs = np.array(range(-ls_alignment_tolerance-1,+ls_alignment_tolerance+1+1,2))



    dither_grid = CartesianGrid(SeparatedCoords((dxs, dxs)))

    E = []
    mean_intensity = []
    plt.figure(figsize=(8, 8))

    for i, (dx, dy) in enumerate(dither_grid.points):
        shifted_lyot_stop = np.roll(np.roll(lyot_stop.shaped, dx, 1), dy, 0).ravel()
        coro = LyotCoronagraph(pupil.grid, focal_plane_mask, shifted_lyot_stop)

        img = 0
        img_ref = 0

        for wl in wavelengths:
            wf = Wavefront(pupil * apodizer, wl)
            img += prop(coro(wf)).intensity
            img_ref += prop(Wavefront(pupil * lyot_stop)).intensity

        x = i % len(dxs)
        y = i // len(dxs)

        plt.subplot(len(dxs), len(dxs), x + (len(dxs) - y - 1) * len(dxs) + 1)

        imshow_field(np.log10(img / max(img_ref)), vmin=-contrast - 1, vmax=-contrast + 4, cmap='inferno')

        frame1 = plt.gca()
        frame1.axes.xaxis.set_ticklabels([])
        frame1.axes.yaxis.set_ticklabels([])

    plt.subplots_adjust(top=0.98, bottom=0.02, left=0.02, right=0.98, wspace=0, hspace=0)

    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

    return {}


'''
TODO:#

def analyze_calculate_throughput(solution_filename):
	pupil, apodizer, focal_plane_mask, lyot_stops, parameters, file_organization = create_coronagraph(solution_filename)
	
'''

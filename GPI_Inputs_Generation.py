from pathlib import Path
from astropy.io import fits
#import gpipsfs
import poppy
import numpy as np
import astropy.units as u

#_grayscale_pixels=4
_default_pupil_array_size = 8*u.m


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
    
    filepath = input_files_dict['directory']
    N = input_files_dict['N']
    lyot_mask = input_files_dict['lyot_stop']['lyot_mask']
    ls_tabs = input_files_dict['lyot_stop']['ls_tabs']
    ls_spid = input_files_dict['lyot_stop']['ls_spid']
    ap_spid = input_files_dict['aperture']['ap_spid']
    ap_sym = input_files_dict['aperture']['ap_sym']


    if ap_sym:
        pup_filename = filepath + 'Primary_GPI_Symmetric_Spid{0:}_N{1:04d}.fits'.format(ap_spid, N)
    else:
        pup_filename = filepath + 'Primary_GPI_Spid{0:}_N{1:04d}.fits'.format(ap_spid, N)

    if ls_spid:
        strut_key = 'struts'
    else:
        strut_key = 'no_struts'

    if ls_tabs:
        ls_filename = filepath +'LS_GPI_{0:s}_ovsamp01_{1:s}_tabs_N{2:04d}.fits'.format(lyot_mask, strut_key, N)
    else:
        ls_filename = filepath + 'LS_GPI_{0:s}_ovsamp01_{1:s}_no_tabs_N{2:04d}.fits'.format(lyot_mask, strut_key, N)

    
    config = Path('masks/' + pup_filename)
    '''
        Checks if aperture file already exists, otherwise creates GPI aperture by calling
        GeminiPrimary() 
        '''

    if config.is_file():
        print('{0:s} exists'.format('masks/' + pup_filename))
    else:
        GPI_primary = GeminiPrimary(force_symmetric=ap_sym).sample(npix=N, grid_size=8.0)

        hdr = fits.Header()
        hdr.set('TELESCOP', 'GPI')
        hdr.set('NAME', 'GEMINI SOUTH PRIMARY')
        hdr.set('D_CIRC', 7.701, 'm: projected diameter of baffle on M2')
        hdr.set('D_INSC', 1.2968, 'm: projected diameter of M2 inner oversized hole')
        #hdr.set('STRUT_W', '0.01, 0.014', 'm: actual support strut widths (laser vane is slightly thicker)')
        #hdr.set('STRUT_ST', -0.2179, 'm: angle lower spiders offset from y')
        #hdr.set('STRUT_AN', 226.9, 'deg: angle lower spiders offset from vertical')

        fits.writeto('masks/'+pup_filename, GPI_primary, header=hdr)
        print('{0:s} has been written to file'.format('masks/' + pup_filename))

    config = Path('masks/' + ls_filename)
    '''
        Checks if lyot stop file already exists, otherwise creates GPI aperture by calling
        GPI_LyotMask() [gpipsfs/gpipsfs/main.py]
        '''

    if config.is_file():
        print('{0:s} exists'.format('masks/' + ls_filename))
    else:
        Lyot_mask = GPI_LyotMask(name=lyot_mask, tabs=ls_tabs).sample(npix=N, grid_size=8.0)
        
        ls_atr = GPI_LyotMask(lyot_mask)
        hdr = fits.Header()
        hdr.set('TELESCOP', 'GPI')
        hdr.set('LS_MASK', ls_atr.name, "name of Lyot mask")
        hdr.set('LS_ID', ls_atr.inner_radius, 'm: Lyot mask inner radius')
        hdr.set('LS_OD', ls_atr.outer_radius, 'm: Lyot mask outer radius')
        hdr.set('NTABS', ls_atr.ntabs, "Number of bad actuator tabs")
        if ls_tabs:
            hdr.set('TAB_DIAM', ls_atr.tabradius * 2, "Tab diameter")
        hdr.set('', ls_atr.support_width, 'm: spider size')
        hdr.set('PUP_DIAM', 9.825e-3, 'm: Gemini pupil size in Lyot plane (without undersizing)')
        hdr.set('MAG', 7.7701 / 0.009825, 'm: magnification between primary & lyot')
        
        fits.writeto('masks/'+ls_filename, Lyot_mask, header=hdr, overwrite=True)
        print('{0:s} has been written to file'.format('masks/' + ls_filename))

    return pup_filename, ls_filename


class GeminiPrimary(poppy.CompoundAnalyticOptic):
    ''' taken from [geminiplanetimager/gpipsfs] '''

    primary_diameter = 7.701  # projected diameter of baffle on M2
    obscuration_diameter = 1.2968  # projected diameter of M2 inner oversized hole
    support_angles = [90 - 43.10, 90 + 43.10, 270 - 43.10, 270 + 43.10]
    support_widths = [0.014, 0.01, 0.01, 0.01]  # laser vane is slightly thicker
    support_offset_y = [0.2179, -0.2179, -0.2179, 0.2179]
    _default_display_size = _default_pupil_array_size  # choose reasonable size for display and sampling.

    def __init__(self, name='Gemini South Primary', undersized=False, force_symmetric=True):
        """Gemini Primary

            :param symmetric: if True, treat all 4 struts as identical. If False, model the fact that
            one of the struts is slightly thicker due to the laser AO beam transfer optics.
            """

        self.force_symmetric = force_symmetric

        outer = poppy.CircularAperture(radius=self.primary_diameter / 2)
        outer._default_display_size = self._default_display_size
        outer.pupil_diam = _default_pupil_array_size  # slightly oversized array

        # Secondary obscuration from pupil diagram provided by Gemini

        sr = self.obscuration_diameter / 2
        if undersized:
            sr = 1.02375 / 2  # SM outer diameter (vs inner hole projected diameter)

        support_width = np.max(self.support_widths) if self.force_symmetric else self.support_widths

        # FIXME could antialias using same technique as used for apodizer grids
        obscuration = poppy.AsymmetricSecondaryObscuration(
            secondary_radius=sr,
            support_angle=self.support_angles,
            support_width=support_width,
            support_offset_y=self.support_offset_y)

        return super(GeminiPrimary, self).__init__(opticslist=[outer, obscuration], name=name)


class GPI_LyotMask(poppy.AnalyticOpticalElement):
    ''' taken from [geminiplanetimager/gpipsfs] '''

    magnification = 7.7701/.009825   # meters at primary/meters at Lyot

    # Locations of bad actuator tabs:
    #   tab locations in millimeters on the physical Lyot mask offset from center
    #   3rd coordinate defines type of linear extension 1:radial, 2:to spider, 0: none
    #  REFERENCE:  CAD drawings and 'mask design notes 2012.pptx' in [GPIMAIN/Instrument docs/coronagraph masks]
    bad_actuator_tab_locations = [ (3.977, 1.932,  1),
                                   (-2.614, 3.750, 1),
                                   (2.614, -3.977, 1),
                                   (-1.932,-1.250, 2),
                                   (1.023, -1.136, 0)]  # coupled to its neighbor, not fully broken
            #Name           R_out  R_in   spider_width  ntabs
    lyot_table = {
            'Blank':       (0,     0,     0,     0, 0 ),
            'Open':        (5.1,   0,     0,     0, 0 ),
            '080m12_02':   (4.786, 1.388, 0.204, 4, 0.328),     # **not installed in GPI**
            '080m12_03':   (4.786, 1.388, 0.304, 4, 0.378),     # used with Y
            '080m12_03_06':(4.786, 1.388, 0.304, 4, 0.528),     # **not installed in GPI**
            '080_04':      (4.786, 1.388, 0.404, 3, 0.428),     # Older generation, does not have 4th tab
            '080m12_04':   (4.786, 1.388, 0.404, 4, 0.428),     # used with J and H
            '080m12_04_c': (4.786, 1.388, 0.604, 5, 0.378),     # extra 5th tab
            '080m12_06':   (4.786, 1.388, 0.604, 4, 0.528),     #
            '080m12_06_03':(4.786, 1.388, 0.604, 4, 0.528),     # exact duplicate of above mask?
            '080m12_07':   (4.686, 1.488, 0.604, 4, 0.528),     #
            '080m12_10':   (4.377, 1.772, 1.004, 4, 0.728) }

    def __init__(self, name='080m12_04', tabs=True, ls_spid=True):

        super(GPI_LyotMask,self).__init__(planetype=poppy.poppy_core._PUPIL, name='GPI Lyot '+name)
        self.pupil_diam=8.0 # default size for display
        # Read in some geometry info from the primary mirror class
        self.support_angles=GeminiPrimary.support_angles
        self.support_offset_y = GeminiPrimary.support_offset_y

        #'name': (primary radius,  secondary radius,  spiderwidth, ntabs, tabwidth)
        # All sizes are in mm on the physical parts; needs to be converted by
        # magnification to get the size reimaged onto the primary
        for k in self.lyot_table.keys():
            if k.upper()==name.upper():
                cname=k
                break
        else:
            raise ValueError('Unknown Lyot stop name: '+name)
        params = self.lyot_table[cname]

        # convert params from mm to meters, then magnify onto the primary
        self.outer_radius = params[0]*1e-3*self.magnification
        self.inner_radius = params[1]*1e-3*self.magnification
        self.support_width = params[2]*1e-3*self.magnification
        self.ntabs= params[3]
        self.tabradius= params[4]*1e-3*self.magnification

        if ls_spid==False: self.support_width = 0
        if tabs==False: self.ntabs=0  # allow user to turn tabs off if desired

        self.wavefront_display_hint = 'intensity' # preferred display for wavefronts at this plane
        #self.wavefront_display_vmax_hint = 3e-9  # Hacky to hard code this, but this value works OK with H_coron sims.

    def get_transmission(self, wave):
        """ Compute the transmission inside/outside of the obscuration

        Based on poppy.AsymmetricSecondaryObscuration but with the bad actuator tabs added.
        """
        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._PUPIL)

        self.transmission = np.ones(wave.shape)

        y, x = self.get_coordinates(wave)

        y *= -1 # Flip Y coordinate convention to match
                # Lyot bad actuator tabs to AOWFS display bad actuators
        r = np.sqrt(x ** 2 + y ** 2)  #* wave.pixelscale

        self.transmission[r < self.inner_radius] = 0
        self.transmission[r > self.outer_radius] = 0

        for angle_deg, offset_y in zip(self.support_angles,
                self.support_offset_y):
            angle = np.deg2rad(angle_deg + 90)  # 90 deg offset is to start from the +Y direction
            width = self.support_width

            # calculate rotated x' and y' coordinates after rotation by that angle.
            # and application of offset
            xp =  np.cos(angle) * (x) + np.sin(angle) * (y-offset_y)
            yp = -np.sin(angle) * (x) + np.cos(angle) * (y-offset_y)

            self.transmission[(xp > 0) & (np.abs(yp) < width / 2)] = 0
            # TODO check here for if there are no pixels marked because the spider is too thin.
            # In that case use a grey scale approximation

        for itab in range(self.ntabs):

            offset_x = self.bad_actuator_tab_locations[itab][0]*1e-3*self.magnification
            offset_y = self.bad_actuator_tab_locations[itab][1]*1e-3*self.magnification
            #print (offset_x, offset_y)

            xo = x-offset_x
            yo = y-offset_y
            r = np.sqrt( xo**2+yo**2)
            self.transmission[r < self.tabradius] = 0

            if self.bad_actuator_tab_locations[itab][2] == 1:
                # Extend tab radially outwards
                angle = np.arctan2(offset_y, offset_x) #-(np.pi/2)
                xp =  np.cos(angle) * xo + np.sin(angle) * yo
                yp = -np.sin(angle) * xo + np.cos(angle) * yo
                self.transmission[(xp > 0) & (np.abs(yp) < self.tabradius )] = 0
            elif self.bad_actuator_tab_locations[itab][2] == 2:
                # special case the one that hangs off the secondary support

                angle=np.deg2rad(self.support_angles[1]+90)
                angle += np.pi/2
                xp =  np.cos(angle) * xo + np.sin(angle) * yo
                yp = -np.sin(angle) * xo + np.cos(angle) * yo
                self.transmission[(xp > 0) & (xp < 0.5) & (np.abs(yp) < self.tabradius )] = 0

        return self.transmission

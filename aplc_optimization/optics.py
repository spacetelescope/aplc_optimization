import hcipy
import numpy as np


class LyotCoronagraphWithFieldStop(hcipy.OpticalElement):
    '''A Lyot coronagraph with a small focal-plane mask and potentially a field stop.

    The area outside of this focal-plane mask is assumed to be fully transmisive. The
    method for propagation is based on [Soummer2007]_.

    .. [Soummer2007] Soummer et al. 2007, "Fast computation of Lyot-style
        coronagraph propagation".

    Parameters
    ----------
    input_grid : Grid
        The grid on which the incoming wavefront is defined.
    focal_plane_mask : Field or OpticalElement
        The (complex) transmission of the focal-plane mask. If this is an :class:`OpticalElement`,
        this will be used instead. This allows for more realistic implementations of focal-plane
        masks.
    lyot_stop : Field or OpticalElement or None
        The (complex) transmission of the Lyot stop. If this is an :class:`OpticalElement`,
        this will be used instead. This allows for more realistic implementations of Lyot stops.
    field_stop : Field or OpticalElement or None
        The (complex) transmission of the field stop. If this is an :class:`OpticalElement`,
        this will be used isntead. This allows for more realistic implementations of the field stop.
        Note: the effective focal plane mask of the coronagraph will be the compound propagation through
        both the FPM and field stop. This means that the field stop should _not_ contain the FPM dot itself.
    focal_length : scalar
        The internal focal length of the Lyot system.
    '''
    def __init__(self, input_grid, focal_plane_mask, lyot_stop=None, field_stop=None, focal_length=1):
        if hasattr(focal_plane_mask, 'input_grid'):
            # Focal plane mask is an optical element.
            fpm_grid = focal_plane_mask.input_grid
        else:
            # Focal plane mask is a field.
            fpm_grid = focal_plane_mask.grid
            focal_plane_mask = hcipy.Apodizer(focal_plane_mask)

        self.prop_fpm = hcipy.FraunhoferPropagator(input_grid, fpm_grid, focal_length=focal_length)

        if lyot_stop is not None and not hasattr(lyot_stop, 'input_grid'):
            lyot_stop = hcipy.Apodizer(lyot_stop)

        if field_stop is not None:
            if hasattr(field_stop, 'input_grid'):
                fs_grid = field_stop.input_grid
            else:
                fs_grid = field_stop.grid
                field_stop = hcipy.Apodizer(field_stop)

            self.prop_fs = hcipy.FraunhoferPropagator(input_grid, fs_grid, focal_length=focal_length)

        self.focal_plane_mask = focal_plane_mask
        self.lyot_stop = lyot_stop
        self.field_stop = field_stop

    def forward(self, wavefront):
        '''Propagate the wavefront through the Lyot coronagraph.

        Parameters
        ----------
        wavefront : Wavefront
            The wavefront to propagate. This wavefront is assumed to be in the pupil plane.

        Returns
        -------
        Wavefront
            The Lyot-plane wavefront.
        '''
        wf_foc = self.prop_fpm.forward(wavefront)
        wf_foc.electric_field -= self.focal_plane_mask.forward(wf_foc).electric_field

        lyot = self.prop_fpm.backward(wf_foc)

        if self.field_stop is None:
            pup = wavefront
        else:
            pup = self.prop_fs.backward(self.field_stop(self.prop_fs.forward(wavefront)))

        # The next line is equal to
        #     lyot.electric_field = pup.electric_field - lyot.electric_field
        # but doesn't copy a numpy array.
        np.subtract(pup.electric_field, lyot.electric_field, out=lyot.electric_field)

        if self.lyot_stop is not None:
            lyot = self.lyot_stop.forward(lyot)

        return lyot

    def backward(self, wavefront):
        '''Propagate the wavefront from the Lyot plane to the pupil plane.

        Parameters
        ----------
        wavefront : Wavefront
            The wavefront to propagate.

        Returns
        -------
        Wavefront
            The pupil-plane wavefront.
        '''
        if self.lyot_stop is not None:
            wf = self.lyot_stop.backward(wavefront)
        else:
            wf = wavefront

        wf_foc = self.prop_fpm.forward(wf)
        wf_foc.electric_field -= self.focal_plane_mask.backward(wf_foc).electric_field

        pup = self.prop_fpm.backward(wf_foc)

        if self.field_stop is not None:
            wf = self.prop_fs.backward(self.field_stop.backward(self.prop_fs.forward(wf)))

        # The next line is equal to
        #     pup.electric_field = wf.electric_field - pup.electric_field
        # but doesn't copy a numpy array.
        np.subtract(wf.electric_field, pup.electric_field, out=pup.electric_field)

        return pup

from astropy.io import fits
import numpy as np


fname_pup  = "masks/Neil_TelAp.fits"
fname_ls   = "masks/Neil_LS.fits"
fname_apod = "masks/Neil_Apod.fits"

TelAp = fits.getdata(fname_pup)
LS    = fits.getdata(fname_ls)
A     = fits.getdata(fname_apod)


# Account for the ratio of the diameter of the square enclosing the aperture to the circumscribed circle
square2circ_ratio = {'luvoir2017novAp05': 0.982, 'HiCAT':1.000}

D       = 0.982
N       = 256
rho_out = 12.0
fp2res  = 50
bw      = 0.10
Nlam    = 4



dx = (D/2)/N
dy = dx
xs = np.matrix(np.linspace(-N+0.5, N-0.5, N)*dx)
ys = xs.copy()
M_fp2 = int(np.ceil(rho_out*fp2res))
dxi = 1./fp2res
xis = np.matrix(np.linspace(-M_fp2+0.5, M_fp2-0.5, 2*M_fp2)*dxi)
etas = xis.copy()
wrs = np.linspace(1.00 - bw/2, 1.00 + bw/2, Nlam)
XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
RRs = np.sqrt(XXs**2 + YYs**2)
p7ap_ind = np.less_equal(RRs, 0.7)

intens_D_0_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
intens_D_0_peak_polychrom = np.zeros((Nlam, 1))
intens_TelAp_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
intens_TelAp_peak_polychrom = np.zeros((Nlam, 1))


for wi, wr in enumerate(wrs):
	
	Psi_D_0 = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp*A*LS[::-1,::-1]), np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
	
	
	
	intens_D_0_polychrom[wi] = np.power(np.absolute(Psi_D_0), 2)
	intens_D_0_peak_polychrom[wi] = (np.sum(TelAp*A*LS[::-1,::-1])*dx*dy/wr)**2
	Psi_TelAp = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp), np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
	intens_TelAp_polychrom[wi] = np.power(np.absolute(Psi_TelAp), 2)
	intens_TelAp_peak_polychrom[wi] = (np.sum(TelAp)*dx*dy/wr)**2

intens_D_0 = np.mean(intens_D_0_polychrom, axis=0)
intens_D_0_peak = np.mean(intens_D_0_peak_polychrom)
intens_TelAp = np.mean(intens_TelAp_polychrom, axis=0)
intens_TelAp_peak = np.mean(intens_TelAp_peak_polychrom)

fwhm_ind_APLC = np.greater_equal(intens_D_0, intens_D_0_peak/2)
fwhm_ind_TelAp = np.greater_equal(intens_TelAp, intens_TelAp_peak/2)

fwhm_sum_TelAp = np.sum(intens_TelAp[fwhm_ind_TelAp])*dxi*dxi
fwhm_sum_APLC = np.sum(intens_D_0[fwhm_ind_APLC])*dxi*dxi
p7ap_sum_TelAp = np.sum(intens_TelAp[p7ap_ind])*dxi*dxi
p7ap_sum_APLC = np.sum(intens_D_0[p7ap_ind])*dxi*dxi

inc_energy = np.sum(np.power(TelAp,2)*dx*dx)
tot_thrupt = np.sum(intens_D_0*dxi*dxi)/np.sum(np.power(TelAp,2)*dx*dx)
fwhm_thrupt = fwhm_sum_APLC/np.sum(np.power(TelAp,2)*dx*dx)
fwhm_circ_thrupt = fwhm_sum_APLC/(np.pi/4)
p7ap_thrupt = p7ap_sum_APLC/np.sum(np.power(TelAp,2)*dx*dx)
p7ap_circ_thrupt = p7ap_sum_APLC/(np.pi/4)
rel_fwhm_thrupt = fwhm_sum_APLC/fwhm_sum_TelAp

#think threl_p7ap thrupt is is the relevant thrpughput metric, double check SCDA survey spreadsheets and paper
rel_p7ap_thrupt = p7ap_sum_APLC/p7ap_sum_TelAp
fwhm_area = np.sum(fwhm_ind_APLC)*dxi*dxi
	
print("////////////////////////////////////////////////////////")
print("{:s}".format(fname_apod))
print("Incident energy on aperture (dimensionless): {:.3f}".format(inc_energy))
print("Non-binary residuals, as a percentage of clear telescope aperture area: {:.2f}%".format(100*tot_thrupt))
print("Band-averaged total throughput: {:.2f}%".format(100*tot_thrupt))
print("Band-averaged half-max throughput: {:.2f}%".format(100*fwhm_thrupt))
print("Band-averaged half-max throughput, circ. ref.: {:.2f}%".format(100*fwhm_circ_thrupt))
print("Band-averaged r=.7 lam/D throughput: {:.2f}%".format(100*p7ap_thrupt))
print("Band-averaged r=.7 lam/D throughput, circ. ref.: {:.2f}%".format(100*p7ap_circ_thrupt))
print("Band-averaged relative half-max throughput: {:.2f}%".format(100*rel_fwhm_thrupt))
print("Band-averaged relative r=0.7 lam/D throughput: {:.2f}%".format(100*rel_p7ap_thrupt))
print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(fwhm_area))	
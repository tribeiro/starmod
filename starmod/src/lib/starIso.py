'''
	Class for isolated stars. Spherical stars.
'''

from baseAstroObj import *
import healpy as H


################################################################################

class StarIsolated(BaseAstroObj):
    ############################################################################

    def __init__(self, i, m, r, prot, nside=16, **kwargs):
        BaseAstroObj.__init__(self, i=i,
                              m=m,
                              r=r,
                              prot=prot,
                              nside=nside,
                              **kwargs)

        # Make surface
        self.npix = H.nside2npix(self.nside)
        self.sresize(self.npix)

        for ipix in range(self.npix):
            self.phi[ipix], self.theta[ipix] = H.pix2ang(self.nside, ipix)
            self.phi[ipix] -= np.pi / 2.0

        self.nx = np.sin(self.phi)
        self.ny = np.cos(self.phi) * np.cos(self.theta)
        self.nz = np.cos(self.phi) * np.sin(self.theta)

        self.px = self.r * self.nx
        self.py = self.r * self.ny
        self.pz = self.r * self.nz

        self.vx = r / prot * np.cos(self.phi)
        self.vy = r / prot * -np.sin(self.phi) * np.cos(self.theta)
        self.vz = r / prot * -np.sin(self.phi) * np.sin(self.theta)

        self.I += 1

    ############################################################################

################################################################################

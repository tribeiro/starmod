'''
	Class for stars in binary. The surface of the star is defined by Roche equipontential surfaces.
'''

from baseAstroObj import *
import healpy as H
from scipy import constants


################################################################################

class StarBinary(BaseAstroObj):
    ############################################################################

    def __init__(self, i, q, ffac=1.0, m=1.0, porb=10., nside=16, gravDark=0.05,
                 **kwargs):
        '''
Initialize StarBinary. User must specify at least inclination and mass ratio
of the system. The definition of the parameters are:
i			:	inclination, in degree (90 -> 0).
q			:	mass ratio (M2/M1)
ffac		:	Roche lobe filling factor, multiply Roche potential by this 
                number. 1.0 for the Roche lobe, larger values for inside the RL.
                Must be >=1.0. default=1.0
m			:	Mass of the star, in Msun. default=1.0Msun
porb		:	Orbital period in hours. default=10.0hours
nside		:	Pixelization number. Must be power of 2, default=16
gravDark	:	Gravity darkening coefficients. default=0.05
[limb_law]	:	Integer selecting a limb-darkening law. [optional]
[limb_coeff]:	Limb-darkening coeficients. [optional]
        '''
        BaseAstroObj.__init__(self, i=i,
                              q=q,
                              ffac=ffac,
                              m=m,
                              porb=porb,
                              nside=nside,
                              gravDark=gravDark,
                              **kwargs)

        self.__MAX_RL_ITER__ = 100
        # Make surface
        self.npix = H.nside2npix(self.nside)
        self.sresize(self.npix)
        self.RL1 = self.rl1()
        self.r = self.RL1

        for ipix in range(self.npix):
            self.phi[ipix], self.theta[ipix] = H.pix2ang(self.nside, ipix)
            self.phi[ipix] -= np.pi / 2.0

        # self.nx = np.sin(self.phi)
        # self.ny = np.cos(self.phi)*np.cos(self.theta)
        # self.nz = np.cos(self.phi)*np.sin(self.theta)

        # rlob = np.zeros(self.npix)
        # niter = np.zeros(self.npix,dtype=int)

        # for ipix in range(self.npix):
        #	rlob[ipix],niter[ipix] = self.rocheLobe(ipix)

        rlob, niter = self.rocheLobeAll()

        # print '# - NITER = %i'%(niter)
        # print len(rlob)
        # print '# - Std. NITER = %i'%(np.std(niter))
        # print '# - Min. NITER = %i'%(np.min(niter))
        # print '# - Max. NITER = %i'%(np.max(niter))

        self.px = rlob * np.sin(self.phi)
        self.py = rlob * np.cos(self.phi) * np.cos(self.theta)
        self.pz = rlob * np.cos(self.phi) * np.sin(self.theta)

        nr = np.sqrt(self.px * self.px + self.py * self.py + self.pz * self.pz)

        n_r = np.sqrt((1. - self.px) * (1. - self.px) + self.py * self.py + self.pz * self.pz)
        nr = nr * nr * nr
        n_r = n_r * n_r * n_r

        self.nx = -self.px / nr + (self.q * (1. - self.px) / n_r) + (self.q + 1) * (self.px) - self.q;
        self.ny = -self.py * (1. / nr + self.q / n_r - (self.q + 1.));
        self.nz = -self.pz * (1. / nr + self.q / n_r);

        n_p = np.sqrt(self.nx * self.nx + self.ny * self.ny + self.nz * self.nz)

        ww = 1. / (self.porb * 3600.)  # frequencia angular
        # sep. orbital
        #
        G_CGS = constants.physical_constants['Newtonian constant of gravitation'][0] * 1e3
        MSOL_CGS = 1.989e33
        #
        a = (G_CGS * MSOL_CGS * self.m * (1 + self.q) / (4. * np.pi * np.pi * ww * ww)) ** (1. / 3.)
        r2 = self.q / (self.q + 1.)  # Dist. orbital da estrela em unidade de a

        self._vrad = 2. * ww * np.pi * a * r2
        self._setvrad = True

        self.vx = -2. * ww * self.py * np.pi * a
        self.vy = 2. * ww * np.pi * a * (self.px - r2)
        self.vz = 0.

        self.I += (n_p / n_p.max()) ** self.gravDark

        #
        # Setup limbdarkeing law
        #
        if hasattr(self, 'limb_law') and hasattr(self, 'limb_coeff'):
            self.limbdarkening = limbDarkeningLaws(self.limbdarkening, self.limb_law, self.limb_coeff)
        elif hasattr(self, 'limb_law') and not hasattr(self, 'limb_coeff'):
            print '[WARNING] - Limb darkenning needs law and coefficients...'
            print '[WARNING] - Law(%i) = %s, needs %i coefficients...' % (self.limb_law,
                                                                          limbDarkeningLawsNames[self.limb_law],
                                                                          limbDarkeningLawsNCoeff[self.limb_law])

    ############################################################################

    def rocheLobe(self, v0, u0):
        '''
        Find radius of the Roche lobe for specified angles.
        '''

        R1 = 0.
        R2 = self.RL1

        a = 2. / ((1. + self.q) * R2)
        b = 2. * self.q / ((1. + self.q) * np.sqrt((R2 * R2) - (2. * R2) + 1.0))
        c = (R2 - (self.q / (1.0 + self.q))) ** 2.0

        fimax = (a + b + c) * self.ffac

        for i in range(self.__MAX_RL_ITER__):

            R = (R2 + R1) / 2.
            r_2 = np.sqrt(R ** 2. - 2. * R * np.sin(u0) + 1.)
            x = R * np.sin(u0)
            y = R * np.cos(u0) * np.cos(v0)

            a = 2. / ((1. + self.q) * R)
            b = 2. * self.q / ((1. + self.q) * r_2)
            c = (x - self.q / (1. + self.q)) ** 2.
            d = y ** 2.

            fi = a + b + c + d;

            # print '%i %f %f %f %f %f %f'%(i,R,R1,R2,fi,fimax,fi-fimax)

            if (fi < fimax):
                R2 = R
            elif (fi > fimax):
                R1 = R

            if (np.abs(fi - fimax) < 1e-5):
                return R

        raise RuntimeError('''Maximum number of interations reached (NITER=%i).
Could not find Roche lobe radius (pix = %i). Pot = %+8.2e''' % (self.__MAX_RL_ITER__, ipix, fi))

    ############################################################################

    def rocheLobeAll(self):
        '''
        Find radius of the Roche lobe for specified angles.
        '''

        v0 = self.theta
        u0 = self.phi

        R1 = np.zeros(self.npix)
        R2 = self.RL1 + np.zeros(self.npix)
        rNotFound = np.zeros(self.npix) == 0  # Mask for pixels that were not found yet

        a = 2. / ((1. + self.q) * R2)
        b = 2. * self.q / ((1. + self.q) * np.sqrt((R2 * R2) - (2. * R2) + 1.0))
        c = (R2 - (self.q / (1.0 + self.q))) ** 2.0

        fimax = (a + b + c) * self.ffac

        for i in range(self.__MAX_RL_ITER__):

            R = (R2 + R1) / 2.
            r_2 = np.sqrt(R ** 2. - 2. * R * np.sin(u0) + 1.)
            x = R * np.sin(u0)
            y = R * np.cos(u0) * np.cos(v0)

            a = 2. / ((1. + self.q) * R)
            b = 2. * self.q / ((1. + self.q) * r_2)
            c = (x - self.q / (1. + self.q)) ** 2.
            d = y ** 2.

            fi = a + b + c + d;

            # print '%i %f %f %f %f %f %f'%(i,R,R1,R2,fi,fimax,fi-fimax)
            m2 = fi <= fimax
            m1 = fi > fimax
            R2[m2] = R[m2]
            R1[m1] = R[m1]

            rNotFound = np.abs(fi - fimax) > 1e-5

            if not rNotFound.any():
                return R, i

        raise RuntimeError('''Maximum number of interations reached (NITER=%i).
Could not find Roche lobe radius.''' % (self.__MAX_RL_ITER__))

    ############################################################################

    def rl1(self):
        '''
        Find distance to internal lagrange point (RL1).
        '''

        x1 = 0.
        x2 = 1.

        for i in range(self.__MAX_RL_ITER__):

            x = ((x2 - x1) / 2.0) + x1
            a = -2. * x / ((1. + self.q) * (x * x) ** (3. / 2.))
            b = -2. * self.q * (x - 1.) / ((1. + self.q) * ((x - 1) ** 2.) ** (3. / 2.))
            c = 2. * (x - (self.q / (1.0 + self.q)))
            dfi = a + b + c

            if dfi >= 0.:
                x2 = x
            elif dfi < 0.:
                x1 = x

            if (np.abs(dfi) < 1.e-5):
                return x

        raise RuntimeError('Could not find RL1...')

    ############################################################################

    def gpole(self):

        qq = self.q
        rlob = self.rocheLobe(np.pi / 2., 0.)
        x = rlob
        n_r = np.abs(1 - x)
        nr = x * x * x
        n_r = n_r * n_r * n_r
        nx = -x / nr + (qq * (1. - x) / n_r) + (qq + 1.) * (x) - qq
        n_p = np.abs(nx);
        return n_p;

    ############################################################################

    def set_radv(self, choice):
        '''
Set/unset radial velocity to the star motion. Usefull when you have 
radial-velocity corrected spectra.
        '''
        if choice:  # and not self._setvrad:
            # RV not set... Set it...
            self.vy -= self._vrad
            self._setvrad = True
        else:  # if not choice and self._setvrad:
            # RV set... unset...
            self.vy += self._vrad
            self._setvrad = False

        ############################################################################

################################################################################

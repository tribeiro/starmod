'''
	Defines basic class for modeling stellar objects.
'''

import numpy as np
import auxiliar as aux

limbDarkeningLawsNames = ['Square-root', 'Linear', 'Quadratic', 'Logarithmic', 'Claret2000']
limbDarkeningLawsNCoeff = [2, 1, 2, 2, 4]


################################################################################

class Surface():
    '''
    Defines the required quantites to describe an object surface.
    '''

    def __init__(self):
        self.theta = np.array([])
        self.phi = np.array([])
        self.px = np.array([])  # posicao x
        self.py = np.array([])  # posicao y
        self.pz = np.array([])  # posicao z
        self.nx = np.array([])  # versor x
        self.ny = np.array([])  # versor y
        self.nz = np.array([])  # versor z
        self.I = np.array([])  # intensidade
        self.vx = np.array([])  # velocide x
        self.vy = np.array([])  # velocide y
        self.vz = np.array([])  # velocide z

    def sresize(self, size, Itype='float'):
        self.theta = np.zeros(size)
        self.phi = np.zeros(size)
        self.px = np.zeros(size)  # posicao x
        self.py = np.zeros(size)  # posicao y
        self.pz = np.zeros(size)  # posicao z
        self.nx = np.zeros(size)  # versor x
        self.ny = np.zeros(size)  # versor y
        self.nz = np.zeros(size)  # versor z
        self.I = np.zeros(size, dtype=Itype)  # intensidade
        self.vx = np.zeros(size)  # velocide x
        self.vy = np.zeros(size)  # velocide y
        self.vz = np.zeros(size)  # velocide z


################################################################################

class BaseAstroObj(Surface):
    '''
    Defines basic astronomical stellar object.
    '''

    def __init__(self, **kwargs):

        Surface.__init__(self)
        for par in kwargs.keys():
            self.__dict__[par] = kwargs[par]

    ############################################################################

    def flux(self, phase):
        '''
        Return the total flux at phase zero.
        '''

        u1 = phase * 2. * np.pi;

        z1 = np.cos(self.i * aux.deg2rad)
        z2 = np.sin(self.i * aux.deg2rad)
        f1 = np.sin(u1)
        f2 = np.cos(u1)

        x1 = z2 * f2
        y1 = f1 * z2

        n_s = self.nx * x1 + self.ny * y1 + self.nz * z1
        mask = n_s > 0.
        return np.sum(self.I[mask] * self.limbdarkening(n_s[mask]))

    ############################################################################

    def vprof(self, phase, vaxis=[]):
        '''
Return the velocity profile for the specified phase and velocity axis. If no axis
is given will try to guess a good one from the information available.
        '''

        u1 = phase * 2. * np.pi;

        z1 = np.cos(self.i * aux.deg2rad)
        z2 = np.sin(self.i * aux.deg2rad)
        f1 = np.sin(u1)
        f2 = np.cos(u1)

        x1 = z2 * f2
        y1 = f1 * z2

        n_s = self.nx * x1 + self.ny * y1 + self.nz * z1
        mask = n_s > 0.

        if len(vaxis) < 3:
            vmax = np.max(np.sqrt(self.vy ** 2. + self.vx ** 2.)) * 1.15
            vaxis = np.arange(-vmax, vmax + vmax / 25., vmax / 25.)

        vflux = np.zeros_like(vaxis)
        vv = x1 * self.vx + y1 * self.vy

        for index in range(len(vaxis) - 1):
            maskv = np.bitwise_and(mask,
                                   np.bitwise_and(vv > vaxis[index],
                                                  vv < vaxis[index + 1]))
            vflux[index] = np.sum(self.I[maskv] * self.limbdarkening(n_s[maskv]))

        return vaxis, vflux

    ############################################################################

    def limbdarkening(self, nv):

        return nv

    ############################################################################

    def makeSpot(self, stheta, sphi, size, val):

        if not (-90. <= sphi <= 90.):
            print '[WARNING] Phi out of range [-90:90]...'
            return -1

        if not (0. <= stheta <= 360.):
            print '[WARNING] Theta out of range [0:360]...'
            return -1

        if size > 90.:
            print '[WARNING] Sport too big [<90.]...'
            return -1

        mdist = self.r * np.abs(np.sin(size * aux.deg2rad))

        px = self.r * np.sin(sphi * aux.deg2rad)
        py = self.r * np.cos(sphi * aux.deg2rad) * np.cos(stheta * aux.deg2rad)
        pz = self.r * np.cos(sphi * aux.deg2rad) * np.sin(stheta * aux.deg2rad)

        dist = np.sqrt((px - self.px) ** 2. + (py - self.py) ** 2. + (pz - self.pz) ** 2.)

        mask = dist < mdist

        self.I[mask] = val

        return 0

        ############################################################################


################################################################################

def limbDarkeningLaws(ldlfunc, law, ldcoeff):
    '''
    Definitions of limb-darkening laws. This is used as a decorator.
    '''
    if len(ldcoeff) != limbDarkeningLawsNCoeff[law] and law >= 0:
        raise IOError('%s law needs %i coeficients. %i Given.' % (limbDarkeningLawsNames[law],
                                                                  limbDarkeningLawsNCoeff[law],
                                                                  len(ldcoeff)))

    if law == 0:
        def ldl_sqrt(nv):
            return (1.0 - ldcoeff[0] * (1.0 - nv) - ldcoeff[1] * (1.0 - np.sqrt(nv)))

        return ldl_sqrt
    elif law == 1:
        def ldl_lin(nv):
            return 1.0 - ldcoeff[0] * (1.0 - nv)

        return ldl_lin
    elif law == 2:
        def ldl_quad(nv):
            return 1. - ldcoeff[0] * (1. - nv) - ldcoeff[1] * (1. - nv) * (1. - nv)

        return ldl_quad
    elif law == 3:
        def ldl_log(nv):
            return 1. - ldcoeff[0] * (1. - nv) - ldcoeff[1] * nv * np.log10(nv)

        return ldl_log
    elif law == 4:
        def ldl_claret2000(nv):
            return 1. - ldcoeff[0] * (1. - nv ** 0.5) - ldcoeff[1] * (1. - nv) - ldcoeff[2] * (1. - nv ** 1.5) - \
                   ldcoeff[3] * (1. - nv * nv)

        return ldl_claret2000
    else:
        return ldlfunc

################################################################################


################################################################################


################################################################################

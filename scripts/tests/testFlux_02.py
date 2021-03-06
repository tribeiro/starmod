#! /usr/bin/env python

import healpy as H
import numpy as np
import pylab as py
from starmod.src.lib.starBin import StarBinary


def main():
    star = StarBinary(90.0, 0.5, nside=64, limb_law=-1, limb_coeff=[0.8])
    Flux0 = star.flux(0.0)
    phases = np.linspace(-0.5, 0.5, 100)

    flux1 = np.zeros(len(phases))
    for i in range(len(phases)):
        flux1[i] = star.flux(phases[i]) / Flux0

    # for tt in np.arange(0,360,80):
    # star.makeSpot(0,-45,10.,0.8)
    # star.makeSpot(45,+00,10.,0.8)
    # star.makeSpot(00, -90, 20., 0.)
    # for theta in range(0,360,45):
    #	star.makeSpot(theta,65,10.,0.8)

    # star.makeSpot(00,-90,20.,0.)
    # star.makeSpot(0,+45,10.,0.8)
    # star.makeSpot(270,-10.,10.,0.8)
    # star.makeSpot(180,-45.,10.,0.8)
    #	star.makeSpot(tt,65.,10.,0.5)

    flux2 = np.zeros(len(phases))
    for i in range(len(phases)):
        flux2[i] = star.flux(phases[i]) / Flux0

    H.mollview(star.I, sub=211, rot=(-90, 90))

    #ff = np.loadtxt('/tmp/cl.dat', unpack=True)

    py.subplot(212)
    py.plot(phases, flux1, '-')
    # py.plot(phases,flux2,'-')
    #py.plot(ff[0], ff[1], '.')
    py.show()


if __name__ == "__main__":
    main()

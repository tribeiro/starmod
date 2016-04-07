#! /usr/bin/env python

import healpy as H
import numpy as np
import pylab as py


def main():
    star = StarBinary(i=80.0, q=0.5, m=1.0, nside=64)
    Flux0 = star.flux(0.0)
    phases = np.linspace(-0.5, 0.5, 50)
    star.set_radv(False)
    # flux1 = np.zeros(len(phases))
    # for i in range(len(phases)):
    #	flux1[i] = star.flux(phases[i])/Flux0


    vv01, pp01 = star.vprof(-0.5)
    py.plot(vv01, -pp01 / Flux0, ls='steps-mid')
    # py.plot(vv02,-pp02/Flux0,ls='steps-mid')
    # py.plot(vv03,-pp03/Flux0,ls='steps-mid')

    py.show()
    return 0
    vprof1 = np.zeros(len(vv01) * len(phases)).reshape(len(phases), len(vv01))
    vprof1[0] = pp01
    for i in range(1, len(phases)):
        print '->', i
        vv01, pp01 = star.vprof(phases[i])
        vprof1[i] = pp01

    # star.makeSpot(00,+90,20.,0.)
    # star.makeSpot(00,-90,20.,0.)
    star.makeSpot(90, 0, 20., 0.)

    vv01, pp01 = star.vprof(-0.5)
    vprof2 = np.zeros(len(vv01) * len(phases)).reshape(len(phases), len(vv01))
    vprof2[0] = pp01
    for i in range(1, len(phases)):
        print '->', i
        vv01, pp01 = star.vprof(phases[i])
        vprof2[i] = pp01

    # vv02,pp02 = star.vprof(0.25)
    # vv03,pp03 = star.vprof(0.75)
    # for tt in np.arange(0,360,80):
    # star.makeSpot(0,-45,10.,0.8)
    # star.makeSpot(45,+00,10.,0.8)

    # for theta in range(0,360,45):
    #	star.makeSpot(theta,65,10.,0.8)

    # star.makeSpot(00,-90,20.,0.)
    # star.makeSpot(0,+45,10.,0.8)
    # star.makeSpot(270,-10.,10.,0.8)
    # star.makeSpot(180,-45.,10.,0.8)
    #	star.makeSpot(tt,65.,10.,0.5)

    H.mollview(star.I, sub=211, rot=(-90, 90))

    # ff = np.loadtxt('/tmp/cl.dat',unpack=True)

    py.subplot(212)

    vprof = vprof1 - vprof2
    # py.imshow(vprof.T,aspect='auto',interpolation='nearest')

    py.plot(vv01, -pp01 / Flux0, ls='steps-mid')
    # py.plot(vv02,-pp02/Flux0,ls='steps-mid')
    # py.plot(vv03,-pp03/Flux0,ls='steps-mid')

    py.show()


if __name__ == "__main__":
    main()

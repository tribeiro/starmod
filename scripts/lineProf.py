#! /usr/bin/env python

'''
	lineProf.py - Generate rotationaly broadened line profiles for Binary star components.
'''

import ConfigParser
import sys
import traceback
from optparse import OptionParser

from astropy.io import fits as pyfits
import pylab as py


################################################################################

def main(argv):
    '''
    Get parameters from a file and generate simulated line profiles for one of
    the components of a binary star. The input file is a
    '''

    parser = OptionParser()
    parser.add_option('-f', '--filename',
                      help='A .cfg file containing the parameters to generate the line profiles.',
                      type='string')
    parser.add_option('-o', '--output',
                      help='The name of the output file. This can be either the root name for ascii files or a single fits files.',
                      type='string')
    parser.add_option('-s', '--show',
                      help='Show results.', action='store_true', default=False)

    opt, args = parser.parse_args(argv)

    # Print header
    print '''
	LINEPROF -	Generate rotationaly broadened line profiles for a binary
				star component.
			
	c - Tiago Ribeiro

	Reading in configuration from %s
''' % (opt.filename)
    #
    # Read input file
    #
    config = ConfigParser.RawConfigParser()

    config.read(opt.filename)

    # Parameters sections
    oPar = 'OrbitalParameters'
    aPar = 'AtmosphericParameters'
    pPar = 'Pixelization'
    iPar = 'InstrumentParameters'
    obsP = 'ObservationParameters'

    # Limb Darkening translator dictionary
    ldLaws = {'None': -1,
              'Square-root': 0,
              'Linear': 1,
              'Quadratic': 2,
              'Logarithmic': 3,
              'Claret2000': 4}
    # Reading in limbdarkening laws and coeficients
    limb_law = ldLaws[config.get(aPar, 'limb_law')]
    ldcoeff = []
    if limb_law > -1:
        ncoeff = limbDarkeningLawsNCoeff[limb_law]
        ldcoeff = np.zeros(ncoeff)
        for lc in range(ncoeff):
            ldcoeff[lc] = config.getfloat(aPar, 'gravDark')

    # Build star
    star = StarBinary(i=config.getfloat(oPar, 'i'),
                      q=config.getfloat(oPar, 'q'),
                      ffac=config.getfloat(oPar, 'ffac'),
                      m=config.getfloat(oPar, 'm'),
                      porb=config.getfloat(oPar, 'porb'),
                      gravDark=config.getfloat(aPar, 'gravDark'),
                      limb_law=limb_law,
                      limb_coeff=ldcoeff,
                      nside=config.getint(pPar, 'nside'))
    # Get flux @ phase zeros for reference
    Flux0 = star.flux(0.0)

    # Build phase list
    phases = []
    if 'phasefile' in config.options(obsP):
        # From a file
        phases = np.loadtxt(config.get(obsP, 'phasefile'))
    else:
        # From star/end/resolution
        phases = np.arange(config.getfloat(obsP, 'startphase'),
                           config.getfloat(obsP, 'endphase'),
                           config.getfloat(obsP, 'phaserel'))

    # Build velocity axis for profiles
    vaxis = []
    try:
        # From star/end/resolution
        vaxis = np.arange(-config.getfloat(iPar, 'minmax') * 1e5,
                          +config.getfloat(iPar, 'minmax') * 1e5,
                          config.getfloat(iPar, 'vrel') * 1e5)
    except:
        # If information is missing
        errinfo = traceback.format_exc(sys.exc_info()[2]).split('\n')
        for ierr in range(len(errinfo)):
            print '[ERROR] - ', errinfo[ierr]

        print '''
[ERROR] - Missing information in configuration file %s.
[ERROR] - Could not build velocity axis. 
[ERROR] - Check %s section of input file
[ERROR]	- STOPING...
''' % (opt.filename, iPar)
        return -1

    print '	Running...'
    flux = np.zeros(len(phases))
    vv01, pp01 = star.vprof(phases[0], vaxis)

    vprof = np.zeros(len(vv01) * len(phases)).reshape(len(phases), len(vv01))
    vprof[0] = pp01
    print '	|' + '-' * np.min([51, len(phases) + 1]) + '|'
    print '	 =',
    for i in range(1, len(phases)):
        vv01, pp01 = star.vprof(phases[i], vaxis)
        vprof[i] = pp01
        sys.stdout.write('=')
        sys.stdout.flush()
        if not i % 50:
            sys.stdout.write('\n	 ')
    print

    # star.makeSpot(00,+90,20.,0.)
    # star.makeSpot(00,-90,20.,0.)
    # star.makeSpot(90,0,20.,0.)

    if opt.output:
        print '\tSaving results to %s ...' % (opt.output)

        prihdr = pyfits.Header()
        sections = config.sections()
        for sct in sections:
            options = config.options(sct)
            for opts in options:
                prihdr[opts[:8]] = config.get(sct, opts)

        prihdu = pyfits.PrimaryHDU(header=prihdr)

        datacolumns = [pyfits.Column(name='vaxis', format='E', array=vaxis)]

        for i in range(len(vprof)):
            datacolumns.append(pyfits.Column(name='prof%04i' % (i + 1), format='E', array=vprof[i]))

        cols = pyfits.ColDefs(datacolumns)
        datahdu = pyfits.new_table(cols)

        phasecolum = pyfits.ColDefs([pyfits.Column(name='phase', format='E', array=phases),
                                     pyfits.Column(name='flux', format='E', array=flux)])

        phasehdu = pyfits.new_table(phasecolum)
        thdulist = pyfits.HDUList([prihdu, phasehdu, datahdu])

        thdulist.writeto(opt.output)

    if opt.show:
        print '\tShowing results ...'
        H.mollview(star.I, sub=211, rot=(-90, 90))

        # ff = np.loadtxt('/tmp/cl.dat',unpack=True)

        py.subplot(212)

        py.imshow(vprof.T, aspect='auto', interpolation='nearest')

        # py.plot(vv01,-pp01/Flux0,ls='steps-mid')
        # py.plot(vv02,-pp02/Flux0,ls='steps-mid')
        # py.plot(vv03,-pp03/Flux0,ls='steps-mid')

        py.show()

    return 0


################################################################################

if __name__ == "__main__":
    main(sys.argv)

    ################################################################################

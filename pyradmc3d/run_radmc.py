import sys
import os
import subprocess

def run_radmc(options, path='~/bin/radmc3d'):
    print options
    return options
    path = os.path.expanduser(path)

    subprocess.Popen([path]+options)

def make_image(lam=10, incl=0, phi=0, nostar=False, sizeau=None, fluxcons=True, secondorder=False, stokes=False):
    pass

def make_cube(incl=0, phi=0, iline=1, widthkms=10, linenlam=10):
    run_radmc(['image'] + ['%s %i' % (k,v) for k,v in locals().iteritems()])

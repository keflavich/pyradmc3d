import numpy as np
from astropy import constants
import itertools

au = constants.au.to('cm').value
msun = constants.M_sun.to('gram').value
rsun = constants.R_sun.to('cm').value

nx,ny,nz = 16,16,16
sizex,sizey,sizez = itertools.cycle((100*au,))

dvdy     = 1e7/(100*au)
rhogas0  = 1e-16
temp0    = 50.
dusttogas= 0.01
vturb0   = 1e4

mstar = msun
rstar = rsun
tstar = 5780
pstar = [0.,0.,0.]

def make_wavelength():
    lambda1 = 0.1
    lambda2 = 7.0
    lambda3 = 25.
    lambda4 = 1.0e4
    n12     = 20
    n23     = 100
    n34     = 30
    lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12)/(1.*n12))
    lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23)/(1.*n23))
    lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34)/(1.*(n34-1.)))
    lam     = np.hstack([lam12,lam23,lam34])
    return lam

def make_grid(sizex,sizey,sizez,nx,ny,nz):
    xi = -sizex + 2*sizex*np.arange(nx+1)/(1.*nx)
    yi = -sizey + 2*sizey*np.arange(ny+1)/(1.*ny)
    zi = -sizez + 2*sizez*np.arange(nz+1)/(1.*nz)
    return xi,yi,zi
    xc = 0.5 * ( xi[:-1] + xi[1:] )
    yc = 0.5 * ( yi[:-1] + yi[1:] )
    zc = 0.5 * ( zi[:-1] + zi[1:] )
    xx = rebin(xc,nx,ny,nz)
    yy = transpose(rebin(yc,ny,nx,nz),[1,0,2])
    zz = transpose(rebin(zc,nz,ny,nx),[2,1,0])
    rrcyl   = sqrt(xx^2+yy^2)

def build_model(nx,ny,nz,rhogas0,temp0,vturb0):
    rhogas  = np.zeros([nz,ny,nx]) + rhogas0
    tgas    = np.zeros([nz,ny,nx]) + temp0
    vx      = np.zeros([nz,ny,nx])
    vy      = dvdy*yy 
    vz      = np.zeros([nz,ny,nx])
    vturb   = np.zeros([nz,ny,nx]) + vturb0

    return rhogas,tgas,vx,vy,vz,vturb

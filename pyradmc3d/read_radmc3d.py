import numpy as np
import os

def read_dust_temperature(filename='dust_temperature',fileformat='b'):

    if fileformat == 'b':
        with open("dust_temperature.bdat","rb") as f:
            header = np.fromfile(f,count=4,dtype='int64')
            data = np.fromfile(f,dtype='float64')
    elif fileformat == 'ascii':
        pass

    return header,data

    #f.readline()

    #ncells = int(f.readline())
    #nspecies = int(f.readline())

    #temperature = []
    #for i in range(nspecies):
    #    temp = empty((ncells,))

    #    for j in range(ncells):
    #        temp[j] = float(f.readline())

    #    temperature.append(temp)

    #f.close()

    #return temperature


def read_spectrum(filename='spectrum.out'):
    wav,flux = np.loadtxt(filename, skiprows=3).T
    return wav,flux

def image(filename=None,ext=None):
    
    nx=0
    ny=0
    nf=0
    sizepix_x=0.e0
    sizepix_y=0.e0
    
    if (filename == None):
        if (ext == None):
            filename = "image.out"
        else:
            filename = "image_"+str(ext)+".out"

    if (os.path.exists(filename) == False):
        print("Sorry, cannot find {0:s}. Presumably radmc2d exited without success. See above for possible error messages of radmc3d!".format(filename))
        return
    else:
        f = open(filename, "r")

    # Read the image.

    iformat = int(f.readline())

    if (iformat < 1) or (iformat > 4):
        print("ERROR: File format of {0:s} not recognized.".format(filename))
        return

    if (iformat == 1) or (iformat == 3):
        radian = (1 == 0)
    else:
        radian = (1 == 1)

    if (iformat == 1) or (iformat == 2):
        stokes = (1 == 0)
    else:
        stokes = (1 == 1)

    nx, ny = tuple(np.array(f.readline().split(),dtype=int))
    nf = int(f.readline())
    sizepix_x, sizepix_y = tuple(np.array(f.readline().split(),dtype=float))

    lam = np.empty(nf)
    for i in range(nf):
        lam[i] = float(f.readline())
    
    f.readline()

    if stokes:
        image = np.empty((ny,nx,4,nf))
    else:
        image = np.empty((nx,ny,nf))

    for i in range(nf):
        for j in range(ny):
            for k in range(nx):
                if stokes:
                    image[j,k,:,i] = np.array(f.readline().split(), dtype=float)
                else:
                    image[j,k,i] = float(f.readline())

                if (j == ny-1) and (k == nx-1):
                    f.readline()

    f.close()

    # Compute the flux in this image as seen at 1 pc.

    if stokes:
        flux = image[:,:,0,:].sum(axis=0).sum(axis=0)
    else:
        flux = image.sum(axis=0).sum(axis=0)

    flux *= sizepix_x*sizepix_y

    if not radian:
        pc = 3.0857200e18
        flux /= pc**2

    # Compute the x and y coordinates

    #x = np.linspace(-(nx-1)/2.,(nx-1)/2.,nx)*sizepix_x
    #y = np.linspace(-(ny-1)/2.,(ny-1)/2.,ny)*sizepix_y
    
    return flux

class image:
    """RADMC-3D CLASS: Image
    """
    def __init__(self):
        self.radian  = 0
        self.nx      = 0
        self.ny      = 0
        self.nf      = 0
        self.sp_x    = 0
        self.sp_y    = 0
        self.Lambda  = 0
        self.image   = np.array([])
        self.flux    = 0.
        self.size_x  = 0.
        self.size_y  = 0.
        self.x       = np.array([])
        self.y       = np.array([])


def readimage(file = "image.out"):
    """RADMC-3D READING FUNCTION: Read an image
    Usage:
          im = readimage(args)
    Args:
         - file: file to read. Default: image.out
    Return:
         im: an image object"""
    im         = image()
    data       = open(file)
    iformat    = int(data.readline())
    im.radian  = iformat==0
    ndim       = data.readline()
    j          = 0
    nx         = np.zeros(2)
    for i in ndim.split():
        nx[j]  = int(i)
        j     += 1
    im.nx      = nx[0]
    im.ny      = nx[1]
    im.nf      = int(data.readline())
    sizepix    = data.readline()
    j          = 0
    sp         = np.zeros(2)
    for i in sizepix.split():
        sp[j]  = float(i)
        j     += 1
    im.sp_x    = sp[0]
    im.sp_y    = sp[1]
    # read EACH lambda, not just the first...
    im.Lambda  = [float(data.readline()) for ii in xrange(im.nf)]
    col        = np.loadtxt(data, usecols=(0,))
    im.image = col.reshape(im.nx,im.ny,im.nf)
    #flux       = utils.compute_flux(im.image,im.nx,im.ny)
    #if(not(im.radian)):
    #    im.flux = flux/pc.val**2
    #im.size_x  = im.sp_x*im.nx
    #im.size_y  = im.sp_y*im.ny
    #im.x       = np.linspace(-im.size_x/2,im.size_x/2,im.nx)
    #im.y       = np.linspace(-im.size_y/2,im.size_y/2,im.ny)
    data.close()
    return(im)


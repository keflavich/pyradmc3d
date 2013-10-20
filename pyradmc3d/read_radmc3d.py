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

def read_population_to_cube(moleculename, shape, column=0):
    """
    Read in a levelpop_[moleculename].dat file and reshape it to a cube
    """
    fn = 'levelpop_%s.dat' % moleculename
    levels = np.loadtxt(fn, usecols=(column,), skiprows=4)
    return levels.reshape(shape)

def image(filename=None,ext=None):
    """
    This function was taken from Kevin's package:
    https://github.com/Kevtron/radmc3d-py
    """
    
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
    (this module was grabbed from Marc Joos' pyradmc package and modified to work with lines a bit)
    https://bitbucket.org/mjoos/pyradmc/src/581fe44e0bb14a02a6a9303e283d527bae109d89/read.py?at=default
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


def radmcimage_to_fits(imagename,
                       fitsname,
                       dpc=1000,
                       arcsec=False,
                       mas=False,
                       restfreq=None, 
                       vel_offset=None,
                       radeg=None,
                       decdeg=None,):
    """
    This function comes from a post by Eric Jensen on the astropy list:
    http://mail.scipy.org/pipermail/astropy/2013-July/002630.html

    Parameters
    ----------
    dpc: float
        distance in pc?
    """

    import numpy as np
    from astropy.io import fits
    
    cc  = 2.9979245800000e10      # Light speed             [cm/s]
    pc  = 3.08572e18              # Parsec                  [cm]
    im  = readimage(imagename)
    pixdeg_x = 180.0*(im['sizepix_x']/(dpc*pc))/np.pi
    pixdeg_y = 180.0*(im['sizepix_y']/(dpc*pc))/np.pi
    freqhz = 1e4*cc/im['lambda']
    image = im['image']


    n_freq = len(freqhz)
    # Because of the way that Python stores arrays, we need to swap around
    # our data array before writing it, so that it is indexed in the order
    # (freq, y, x), rather than (x, y, freq) - the last index in the Numpy
    # array will become AXIS1 in the FITS image, and so on. 
    if im['stokes']:
        (nx, ny, nstokes, nfreq) = np.shape(image)
        newimage = np.transpose(image, (3, 2, 1, 0))
    else:
        (nx, ny, nfreq) = np.shape(image)
        newimage = np.transpose(image, (2, 1, 0))
        
    #
    # Compute the conversion factor from erg/cm^2/s/Hz/ster to Jy/pixel
    #
    pixsurf_ster = pixdeg_x*pixdeg_y * (np.pi/180.)**2
    factor = 1e+23 * pixsurf_ster
    # And scale the image array accordingly:
    image_in_jypix = factor * newimage
    
    #
    # Make FITS header information:
    #
    header = fits.Header()
    header['BTYPE'] = 'Intensity'
    header['BSCALE'] = 1
    header['BZERO'] = 0
    header['BUNIT'] = 'JY/PIXEL'
    
    if radeg is not None:
        header['EPOCH'] = 2000.
        header['LONPOLE'] = 180.
        header['CTYPE1'] = 'RA---SIN'
        header['CRVAL1'] = radeg

    # Figure out what units we are using per pixel:
    if arcsec:
        unit = 'arcsec'
        multiplier = 3600.
    elif mas:
        unit = 'mas'
        multiplier = 3.6e6
    else:
        unit = 'deg'
        multiplier = 1

    header['CDELT1'] = multiplier*pixdeg_x
    header['CUNIT1'] = unit
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX1'] = 1.0*((nx+1)/2)

    if decdeg is not None:
        header['CTYPE2'] = 'DEC--SIN'
        header['CRVAL2'] = decdeg

    header['CDELT2'] = multiplier*pixdeg_y
    header['CUNIT2'] = unit
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX2'] = 1.0* ((ny+1)/2)

    # If the keyword is set for rest frequency, add that to the header as
    # it gives a reference point for, e.g., the velocity of the frequency
    # steps relative to some line rest frequency.  Note that even though
    # the keyword used here is 'restfreq', the FITS standard for this
    # field is RESTFRQ, i.e. with no second 'E'. 
    if restfreq is not None:
        header['RESTFRQ'] =  (float(restfreq), 'Rest frequency of this transition')
        
    #
    # ...Frequency
    #
    #  If they have specified a velocity offset, we shift all of the
    #  specified frequencies by that amount.  This allows, e.g., mapping a
    #  transition of a particular CO line in radmc3d, but then accounting
    #  for the fact that a real object might have a center-of-mass radial
    #  velocity offset.  Note that we do *not* shift the input 'restfreq'
    #  keyword, as that is presumed to be the real frequency of that
    #  transition.
    if vel_offset is not None:
        # Assume input velocity is in km/s, so
        # specify speed of light in those
        # units:
        c  = 2.9979245800000e5        # Light speed             [km/s]
        freqhz *= (1. - vel_offset/c)
        header['COMMENT'] = "Velocity offset of %0.2f km/s applied to model." % vel_offset

    if n_freq > 1:
        # multiple frequencies - set up the header keywords to define the
        #    third axis as frequency
        header['CTYPE3'] = 'FREQ'
        header['CUNIT3'] = 'Hz'
        header['CRPIX3'] = 1.0
        header['CRVAL3'] = freqhz[0]
        # Calculate the frequency step, assuming equal steps between all:
        delta_freq = freqhz[1] - freqhz[0]
        header['CDELT3'] = delta_freq
    else:                # only one frequency
        header['RESTFREQ'] = freqhz

    # Make a FITS file!
    #

    fits.writeto(fitsname, image_in_jypix, header, output_verify='fix')


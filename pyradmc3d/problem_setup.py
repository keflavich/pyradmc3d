from astropy import units as u
from astropy import constants
import numpy as np

def write_field(filename, 
                suffix,
                data,
                format='b',
                formatnumber=1,
                precision='double',
                nspecies=1,
                is3d=True,
                size=None):
    """
    Create a field file, e.g. dust_density, co_density, dust_temperature 

    Parameters
    ----------
    filename: str
        'dust_density'
        'dust_temperature'
    suffix: str
        'dat' or 'inp'
    format: 'b','u', or 'ascii'
        'b' is binary, 'ascii' is for plain ascii,
        and 'u' is for f77-streaming-record style
        ('u' is not implemented)
    formatnumber: int
        Should be 1 according to the 0.35 docs
    precision: 'single' or 'double'
        Write as float or double? (only matters for binary formats)
    nspecies: int
        Number of dust species.  Right now, only 1 dust species is supported.
        For molecules, write separate files.
    is3d: bool
        Check whether the data is 3d before progressing
    size: int
        Optional - specify the size directly *if* the data is not 3d
    """
    if is3d:
        assert data.ndim == 3
        nz,ny,nx = data.shape
        size = data.size
    elif size is None:
        raise ValueError("If is3d is False, you must specify the number of cells directly with size.")

    if nspecies != 1:
        raise ValueError("Only 1 dust species is implemented.")

    if format=='ascii':
        outfile = filename+"."+suffix
        with open(outfile,'w') as f:
            f.write('1\n')
            f.write('%d\n' % (size))
            f.write('1\n')
            for iz in xrange(nz):
                for iy in xrange(ny):
                    for ix in xrange(nx):
                        f.write('%e\n' % data[iz,iy,ix])
    elif format=='b':
        outfile = filename+".b"+suffix
        header = np.array([formatnumber,
                           4 if precision == 'single' else 8,
                           size,
                           nspecies],
                          dtype='int64')
        dtype = 'float32' if precision == 'single' else 'float64'

        with open(outfile,'wb') as f:
            header.tofile(f)
            # is this right?  Do we need to flatten first?
            data.ravel().astype(dtype).tofile(f)
    elif format=='u':
        raise NotImplementedError("FORTRAN77 style 'records' have not been implemented")

def write_vfield(vx,vy,vz,format='b'):
    """
    """
    assert vx.size == vy.size == vz.size
    data = np.array(zip(vx.flat,vy.flat,vz.flat))
    if format=='b':
        write_field('gas_velocity','inp',data,is3d=False,size=vx.size)
    elif format=='ascii':
        outfile = 'gas_velocity.inp'
        with open(outfile,'w') as f:
            f.write('1\n')
            f.write('%d\n' % (vx.size))
            f.write('1\n')
            for iz in xrange(vz.size):
                f.write('%e %e %e\n' % tuple(data[iz,:].tolist()))


def write_grid(dimensions, xyzdims=None, format='b'):
    """
    """
    # format, gridstyle, coordsystem, gridinfo
    # include x,y,z
    header = np.array([1,0,1,0,
                       1,1,1]+
                      list(dimensions),
                      dtype='int64')

    if xyzdims is None:
        xdims = np.arange(dimensions[2]+1,dtype='float64')
        ydims = np.arange(dimensions[1]+1,dtype='float64')
        zdims = np.arange(dimensions[0]+1,dtype='float64')
    else:
        xdims,ydims,zdims = xyzdims

    if format=='b':
        outfile = 'amr_grid.binp'
        with open(outfile,'wb') as f:
            header.tofile(f)
            xdims.tofile(f)
            ydims.tofile(f)
            zdims.tofile(f)
    elif format=='ascii':
        outfile = 'amr_grid.inp'
        with open(outfile,'w') as f:
            f.write("1\n0\n1\n0\n1 1 1\n")
            f.write(" ".join(['%i' % d for d in dimensions])+"\n")
            f.write(" ".join(["%f" % x for x in xdims])+"\n")
            f.write(" ".join(["%f" % x for x in ydims])+"\n")
            f.write(" ".join(["%f" % x for x in zdims])+"\n")

def write_radmc3d(nphot=1e5,nphot_spec=1e5,out_style='b',do_lines=False,lines_mode=3,
        linesnonlteconvcrit=1e-2,
        linesnonltemaxiter=10):
    with open('radmc3d.inp','w') as f:
        f.write('nphot = %i\n' % nphot)
        #f.write('scattering_mode_max = 1\n')
        f.write('nphot_spec = %i\n' % (nphot_spec))
        #f.write('incl_userdef_srcalp = 1\n')
        #f.write('incl_freefree = 1\n')
        #f.write('incl_dust = 1\n')
        if out_style == 'b':
            f.write('rto_style = 3\n')
        if do_lines:
            f.write('incl_lines = 1\n')
        if lines_mode is not None:
            f.write('lines_mode = %i\n' % lines_mode)

def write_lines(species=['co'],colliders={'co':['h2']}, db={'co':'leiden'}):
    with open('lines.inp','w') as f:
        f.write('2\n')
        f.write('%i\n' % (len(species)))
        for s in species:
            f.write("%s    %s  0 0 %i\n" % (s, db[s], len(colliders[s])))
            for c in colliders[s]:
                f.write("%s\n" % c)

        f.write('\n')

def write_dustopac(species_names=['silicate'],formatnumber=2,readmethods=[1],quantum=[0]):
    """

    readmethod: int
        * 1: dustkappa_<name>.inp
        * 10: dustkapscatmat_<name>.inp
        * -1: dustopac_<name>.inp
    quantum: 0
        Must be zero, but left flexible in case quantum grain heating is
        implemented in the future.
    """
    header = """
{formatnumber}               Format number of this file
{nspecies}               Nr of dust species
============================================================================
    """.strip().lstrip()
    body = """
{0}               Way in which this dust species is read
{1}               0=Thermal grain
{2}        Extension of name of dustkappa_***.inp file
----------------------------------------------------------------------------
""".strip().lstrip()

    with open('dustopac.inp','w') as f:
        output = header.format(formatnumber=formatnumber, nspecies=len(species_names))
        output += "\n"
        output += "\n".join([body.format(*x) for x in zip(readmethods, quantum, species_names)])
        f.write(output)
        f.write("\n")

def make_default_wavelength():
    lambda1=0.1
    lambda2=7.0
    lambda3=25.0
    lambda4=0.8e4
    n12=20
    n23=100
    n34=30
    lam12 = [lambda1*(lambda2/lambda1)**(i/(1.0*n12)) for i in range(n12)]
    lam23 = [lambda2*(lambda3/lambda2)**(i/(1.0*n23)) for i in range(n23)]
    lam34 = [lambda3*(lambda4/lambda3)**(i/(1.0*n34)) for i in range(n34)]
    [lam12.append(lam23[i]) for i in range(n23)]
    [lam12.append(lam34[i]) for i in range(n34)]
    return lam12

def write_wavelength():
    with open('wavelength_micron.inp','w') as f:
        lam = make_default_wavelength()
        f.write('%i\n' % len(lam))
        f.write("\n".join(["%f" % l for l in lam]))
        f.write("\n")

def write_external_bb(lam, tem):
    intensity = (2*constants.h*constants.c**2/(lam*u.micron)**5 *
                (np.exp(constants.c*constants.h/(lam*u.micron * constants.k_B * tem * u.K)).to(1) - 1)**(-1)).value
    with open('external_source.inp','w') as f:
        f.write('2\n')
        f.write('%i\n' % len(lam))
        f.write("\n".join(["%f" % l for l in lam]))
        f.write("\n")
        f.write("\n".join(["%g" % i for i in intensity]))
        f.write("\n")

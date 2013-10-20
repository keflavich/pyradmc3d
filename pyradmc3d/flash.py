"""
Tools to turn a fixed-grid FLASH simulation into the format appropriate
for RADMC-3D 0.35 (October 19, 2013)
"""
import h5py

def load_field(fn, kw='dens'):
    """
    Load a field from an HDF5 dataset into an array
    
    Parameters
    ----------
    fn: str
        HDF5 dataset filename
    kw: str
        Field name

    Examples
    --------
    >>> data = load_field('sim_dens','dens')
    """
    f = h5py.File(fn)
    data = f[kw]
    return data[:]


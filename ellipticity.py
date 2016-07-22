#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table

from lsst.validate.drp import validate, util
import lsst.daf.persistence as dafPersist


def calculate_ellipticity(data):
    I_xx = data['base_SdssShape_xx']
    I_xy = data['base_SdssShape_xy']
    I_yy = data['base_SdssShape_yy']
    e = (I_xx - I_yy + 2j*I_xy) / (I_xx+I_yy + 2*np.sqrt(I_xx*I_yy - I_xy*2))
    e1 = np.imag(e)
    e2 = np.real(e)
    return e, e1, e2


def plot_ellipticity(e1, e2):
    plt.scatter(e1, e2)
    plt.xlabel('e1')
    plt.ylabel('e2')
#    plt.show()


def plot_spatial_ellipticity(x, y, e1, e2):
    """Plot a quiver plot of e1/e2 vs. x, y"""
    plt.quiver(x, y, e1, e2)
    plt.xlabel('x')
    plt.ylabel('y')
#    plt.show()


def plot_angular_ellipticity(x, y, e1, e2):
    """Plot a histogram of the e1/e2 vs. x, y"""
    plt.hist(np.arctan(e1/e2), bins=np.linspace(-np.pi/2,+np.pi/2,50), range=(-np.pi/2,+np.pi/2),
             histtype='stepfilled')
    plt.xlabel('atan(e1/e2)')
#    plt.show()


def get_psf(filename):
    """Get the PSF for an exposure, x, y"""
    pass


def subtract_psf(data, psf):
    """Subtract the PSF e1, e2 from the measure e1, e2

    Generate new m_e1, m_e2
    """
    data['m_e1'] = np.array([d['e1'] - psf(d) for d in data])
    return data


if __name__ == "__main__":
    repo = "/Users/wmwv/local/lsst/validation_data_cfht/data/"
    filename = "/Users/wmwv/local/lsst/validation_data_cfht/data/src/06AL01/D3/2006-05-20/r/SRC-849375-21.fits"
    visits = [176837, 176846]
    dId = {'visit': 849375, 'filter': 'r', 'ccd': 12}
    butler = dafPersist.Butler(repo)

    data = butler.get("src", dId)
#    psf = butler.get("psf", dId)

#    data = Table.read(filename)
    ra, dec = data['coord_ra'], data['coord_dec']
    e, e1, e2 = calculate_ellipticity(data)
    ellipticity_norm = np.real(e * np.conj(e))
    w, = np.where(data['base_PsfFlux_flux']/data['base_PsfFlux_fluxSigma'] > 100)
    print(np.median(ellipticity_norm[w]))

    plt.subplot(2,2,1)
    plot_ellipticity(e1, e2)
    plt.xlim(-1,+1)
    plt.ylim(-1,+1)
    plt.text(-0.9, +0.9, "{:.5f}".format(np.median(ellipticity_norm[w])))
    plt.subplot(2,2,2)
    plot_spatial_ellipticity(data['base_SdssCentroid_x'], data['base_SdssCentroid_y'], e1, e2)
    plt.subplot(2,2,3)
    plot_angular_ellipticity(data['base_SdssCentroid_x'], data['base_SdssCentroid_y'], e1, e2)

"""Variable definitions and metadata for the GRASP-kernel netCDF files (V5 onward).

Used by convertkernels.py (when writing a file) and updatekernelmeta.py (when updating the
metadata of an existing file), so that both produce identical descriptions.

Conventions (see also the 'efficiency_convention' global attribute):
- x = 2*pi*r_v/lambda, r_v the volume-equivalent sphere radius.
- ext, abs, sca = C/V (1/um at the reference wavelength).
- qext, qsca, qabs, qb are relative to pi*r_v^2 (Mishchenko; Dubovik et al. 2006; GRASP),
  not to the particle's orientation-averaged projected area A; C/A = q*sphericity.
- scama is dimensionless: 0.5*integral(F11 sin(theta) dtheta) = 1, i.e. the integral over all
  solid angles is 4*pi; the phase function per steradian is F11/(4*pi).
"""
import numpy as np

VERSION = 'V5'

EFFICIENCY_CONVENTION = (
    'qext, qsca, qabs and qb are cross sections divided by pi*r_v^2, the geometric cross section of the '
    'volume-equivalent sphere (Mishchenko; Dubovik et al. 2006; GRASP). For non-spherical particles this is not '
    'C/A with A the orientation-averaged projected area (the convention of TAMUdust2020, Saito et al. 2021): '
    'qext tends to 2/sphericity at large x, and C/A = q*sphericity.')
SIZE_CONVENTION = 'x = 2*pi*r_v/lambda, with r_v the radius of the sphere of equal volume.'
WAVELENGTH_SCALING = (
    'ext, abs, sca (1/um) and cext, csca, cabs (um^2) are given at reference_wavelength_um. At another wavelength '
    'lambda and the same x, multiply ext, abs and sca by reference_wavelength_um/lambda and cext, csca and cabs by '
    '(lambda/reference_wavelength_um)^2. All other variables depend on x only.')


def var_attrs(lamb):
    """{variable: {attribute: value}} for every variable written by convertkernels.py."""
    L = '%.3f um' % lamb
    a = {}
    a['ratio'] = dict(units='1', long_name='shape label', description=(
        'Spheroid aspect ratio (rotational/equatorial semi-axis; < 1 oblate, > 1 prolate) for single spheroid shapes, '
        '1 for spheres. For pre-integrated or non-spheroidal shapes (the GRASP spheroid mixture, the TAMUdust2020 '
        'hexahedra) it is the ratio label of the GRASP kernel files (2.986), not a property of the particles. '
        'See shape_description and sphericity.'))
    a['shape_description'] = dict(long_name='description of each shape (ratio) entry')
    a['sphericity'] = dict(units='1', long_name='sphericity', description=(
        'Psi = pi^(1/3)*(6V)^(2/3)/S_surface = pi*r_v^2/A for convex particles (A orientation-averaged projected '
        'area = S_surface/4); 1/Psi = A/(pi*r_v^2). For an ensemble, A and V are ensemble averages.'))
    a['mr'] = dict(units='1', long_name='real part of the refractive index')
    a['mi'] = dict(units='1', long_name='imaginary part of the refractive index',
                   description='Stored negative (GRASP convention); the absorbing part is |mi|.')
    a['x'] = dict(units='1', long_name='size parameter of the volume-equivalent sphere', description=SIZE_CONVENTION)
    a['angle'] = dict(units='degree', long_name='scattering angle')
    for se in ('scattering_element_index', 'scattering_element'):
        a[se] = dict(units='1', long_name='scattering matrix element',
                     description='Order of the elements in scama: 11, 22, 33, 44, 12, 34.')
    for short, full in (('ext', 'extinction'), ('abs', 'absorption'), ('sca', 'scattering')):
        a[short] = dict(units='1/um', long_name='volume %s coefficient' % full, description=(
            '%s cross section per unit particle volume, C_%s/V, at %s; equivalently the %s coefficient for unit volume '
            'concentration. At another wavelength lambda and the same x multiply by reference_wavelength_um/lambda. '
            'Mass %s efficiency = %s/rho (numerically, m^2/g = %s[1/um] / rho[g/cm^3]).'
            % (full.capitalize(), short, L, full, full, short, short)))
        a['q' + short] = dict(units='1', long_name='%s efficiency relative to the volume-equivalent sphere' % full,
                              description=(
            'C_%s/(pi*r_v^2) = %s*4*r_v/3 (Mishchenko/Dubovik/GRASP convention). For non-spherical particles this is not '
            'C/A with A the orientation-averaged projected area: qext tends to 2/sphericity at large x, and '
            'C_%s/A = q%s*sphericity.' % (short, short, short, short)))
        a['c' + short] = dict(units='um^2', long_name='%s cross section' % full, description=(
            'C_%s at %s = %s*(4/3)*pi*r_v^3 = q%s*pi*r_v^2. At another wavelength lambda and the same x multiply by '
            '(lambda/reference_wavelength_um)^2.' % (short, L, short, short)))
    a['qb'] = dict(units='1', long_name='backscattering efficiency relative to the volume-equivalent sphere', description=(
        'qsca*F11(180 deg), with F11 normalized so that its integral over all solid angles is 4*pi (see scama); '
        'dimensionless. The backscattering efficiency per steradian is qb/(4*pi) [1/sr], and the backscattering cross '
        'section per steradian is qb*pi*r_v^2/(4*pi).'))
    a['lidar_ratio'] = dict(units='sr', long_name='lidar ratio', description='Extinction-to-backscatter ratio, 4*pi*qext/qb.')
    a['g'] = dict(units='1', long_name='asymmetry parameter',
                  description='Mean cosine of the scattering angle, 0.5*integral(F11 cos(theta) sin(theta) dtheta).')
    a['scama'] = dict(units='1', long_name='normalized scattering matrix elements', description=(
        'Elements F_ij in the order of scattering_element_index (11, 22, 33, 44, 12, 34), normalized so that '
        '0.5*integral(F11 sin(theta) dtheta) from 0 to pi = 1, i.e. the integral over all solid angles is 4*pi; '
        'dimensionless (F11/(4*pi) is the phase function per steradian). F12 is negative in the Rayleigh limit; F34 '
        'has the sign of the GRASP kernels. Below x = renormUpperX (global attribute) all elements are renormalized '
        'so that F11 meets this normalization.'))
    return a


def spheroid_sphericity(e):
    """Sphericity of a spheroid with axis ratio e (rotational/equatorial semi-axis)."""
    if abs(e - 1) < 1e-9:
        return 1.0
    a = e ** (-1 / 3)
    c = e * a  # volume-equivalent radius 1
    if e > 1:
        ecc = np.sqrt(1 - a * a / (c * c))
        S = 2 * np.pi * a * a * (1 + c / (a * ecc) * np.arcsin(ecc))
    else:
        ecc = np.sqrt(1 - c * c / (a * a))
        S = 2 * np.pi * a * a * (1 + (1 - ecc ** 2) / ecc * np.arctanh(ecc))
    return 4 * np.pi / S


def mixture_sphericity(shape_dist_file):
    """Sphericity of a spheroid mixture at fixed r_v, from a GEOSmie kernel_shape_dist file
    (columns: fraction, aspect ratio): 1 / sum(w_k / Psi_k)."""
    d = np.loadtxt(shape_dist_file)
    w = d[:, 0] / d[:, 0].sum()
    return 1.0 / np.sum([wk / spheroid_sphericity(ek) for wk, ek in zip(w, d[:, 1])])


def set_attrs(var, attrs):
    for k, v in attrs.items():
        if k in var.ncattrs():
            var.delncattr(k)
        var.setncattr(k, v)


def apply(ncdf, lamb, ratio_values, sphericity, shape_description, global_attrs):
    """Write the V5 metadata into an open netCDF4 Dataset: per-variable attributes, the ratio values,
    the sphericity and shape_description variables, and global attributes."""
    ncdf.variables['ratio'][:] = np.asarray(ratio_values, float)
    if 'sphericity' not in ncdf.variables:
        ncdf.createVariable('sphericity', 'f8', ('ratio',))
    ncdf.variables['sphericity'][:] = np.asarray(sphericity, float)
    if 'shape_description' not in ncdf.variables:
        ncdf.createVariable('shape_description', str, ('ratio',))
    ncdf.variables['shape_description'][:] = np.array(shape_description, dtype=object)
    for name, attrs in var_attrs(lamb).items():
        if name in ncdf.variables:
            set_attrs(ncdf.variables[name], attrs)
    g = {'title': 'GRASP-format single-scattering kernels', 'version': VERSION,
         'reference_wavelength_um': float(lamb), 'size_convention': SIZE_CONVENTION,
         'efficiency_convention': EFFICIENCY_CONVENTION, 'wavelength_scaling': WAVELENGTH_SCALING}
    g.update(global_attrs)
    for k, v in g.items():
        if k in ncdf.ncattrs():
            ncdf.delncattr(k)
        ncdf.setncattr(k, v)

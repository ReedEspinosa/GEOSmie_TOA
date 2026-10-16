"""Write a V5 copy of an existing GRASP-kernel netCDF file with updated metadata only.

  python updatekernelmeta.py --filename input_files/kernelconversion/<spec>.json

The spec JSON has:
  src, dest        input (e.g. V4) and output netCDF paths; dest must not exist
  ratioValues      values for the 'ratio' variable (e.g. exact axis ratios from the kernel headers)
  sphericity       one entry per ratio: a number, "spheroid" (analytic from the ratio value), or
                   "mixture:<kernel_shape_dist file>" (spheroid mixture at fixed r_v)
  shapeDescription one string per ratio, or "spheroid" to generate "spheroid, axis ratio <value>"
  globalAttrs      extra global attributes (e.g. history)

The data are copied unchanged, except that 'ratio' takes ratioValues and 'lidar_ratio'
(never written by earlier converter versions) is filled with 4*pi*qext/qb. Variable
attributes, the new 'sphericity' and 'shape_description' variables and the global attributes
come from kernelmeta.py, the same definitions convertkernels.py uses.
"""
import datetime
import json
import os
import re
import shutil
from optparse import OptionParser

import netCDF4 as nc
import numpy as np

import kernelmeta
from convertkernels import gitCommit


def sphericity_of(entry, ratio_value):
    if isinstance(entry, str):
        if entry == 'spheroid':
            return kernelmeta.spheroid_sphericity(ratio_value)
        if entry.startswith('mixture:'):
            return kernelmeta.mixture_sphericity(entry.split(':', 1)[1])
        raise ValueError('unknown sphericity entry %r' % entry)
    return float(entry)


def fun(ifn):
    with open(ifn) as fp:
        spec = json.load(fp)
    src, dest = spec['src'], spec['dest']
    if os.path.exists(dest):
        raise FileExistsError(dest)
    rv = [float(v) for v in spec['ratioValues']]
    sph = [sphericity_of(e, v) for e, v in zip(spec['sphericity'], rv)]
    desc = spec['shapeDescription']
    if desc == 'spheroid':
        desc = ['spheroid, axis ratio %g (rotational/equatorial semi-axis)' % v for v in rv]
    print('copying %s -> %s' % (src, dest))
    shutil.copyfile(src, dest)
    with nc.Dataset(dest, 'a') as ds:
        if len(ds.dimensions['ratio']) != len(rv):
            raise ValueError('ratio dimension has %d entries, spec has %d' % (len(ds.dimensions['ratio']), len(rv)))
        m = re.search(r'λ=([0-9.]+)', ds['cext'].getncattr('long_name'))  # V4 files record λ_ref only here
        lamb = float(spec.get('referenceWavelength', m.group(1) if m else 0.34))
        qext, qb = ds['qext'][:].filled(np.nan), ds['qb'][:].filled(np.nan)
        ds['lidar_ratio'][:] = 4 * np.pi * qext / np.where(qb > 0, qb, np.nan)
        g = {'source_file': os.path.abspath(src), 'geosmie_commit': gitCommit(),
             'created': datetime.datetime.now().isoformat(timespec='seconds')}
        g.update(spec.get('globalAttrs', {}))
        kernelmeta.apply(ds, lamb, rv, sph, desc, g)
    return dest


if __name__ == '__main__':
    parser = OptionParser(usage='Usage: %prog --filename <spec.json>')
    parser.add_option('--filename', dest='filename', default='', help='metadata update spec (JSON)')
    (options, args) = parser.parse_args()
    if not os.path.exists(options.filename):
        parser.error('Input file path (--filename) does not exist')
    fun(options.filename)
    print('Done!')

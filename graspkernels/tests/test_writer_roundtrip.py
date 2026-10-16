"""Writer check: re-write a set of GRASP text kernels from the values read by GEOSmie's
convertkernels.py and compare with the originals.

  python graspkernels/tests/test_writer_roundtrip.py <kernel dir> <scratch dir> [ratio]

Passes if every rewritten file equals the original line by line (ignoring trailing
blanks), or, where layouts differ, if all numbers agree exactly.
"""
import os
import re
import sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
import convertkernels as ck
from graspkernels import graspkerneltext as gkt

src, tmp = sys.argv[1], sys.argv[2]
ratio = sys.argv[3] if len(sys.argv) > 3 else '100'
r, wl, ang = gkt.read_grid(os.path.join(src, 'grid1.dat'))
mr, mi = gkt.read_ri_nodes(os.path.join(src, 'kernels_%s_00.txt' % ratio))
nr, nang = len(r), len(ang)
cj = gkt.conversion_json('x', src, 'kernels', ratio, nr, len(mr), len(mi), nang, [])
e, a, _ = ck.read00(src, ratio, 'kernels', cj['contlen'])
shape = (len(mr), len(mi), nr)
ext, absb = np.array(e).reshape(shape), np.array(a).reshape(shape)
pm = {}
for el in ['11', '12', '22', '33', '34', '44']:
    s = ck.readScatEle(src, ratio, el, 'kernels', cj['scathdrlen'], cj['scatcontlen'], cj['scatelemlen'])
    pm[el] = np.array(s).reshape(shape + (nang,))
ratio_value = float(open(os.path.join(src, 'kernels_%s_00.txt' % ratio)).readline().split()[2])
gkt.write_kernels(tmp, 'kernels', ratio, ratio_value, wl, r, mr, mi, ext, absb, ang, pm)
gkt.write_grid(tmp, r, wl, ang)

num = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+)(?:[DE][-+]?\d+)?')
ok = True
for fn in ['grid1.dat'] + ['kernels_%s_%s.txt' % (ratio, el) for el in ['00', '11', '12', '22', '33', '34', '44']]:
    A = [l.rstrip() for l in open(os.path.join(src, fn))]
    B = [l.rstrip() for l in open(os.path.join(tmp, fn))]
    while A and not A[-1]:
        A.pop()
    while B and not B[-1]:
        B.pop()
    if A == B:
        print('%-22s identical (%d lines)' % (fn, len(A)))
        continue
    na = [float(t) for t in num.findall('\n'.join(A))]
    nb = [float(t) for t in num.findall('\n'.join(B))]
    same = len(na) == len(nb) and np.array_equal(na, nb)
    ndiff = sum(x != y for x, y in zip(A, B)) + abs(len(A) - len(B))
    print('%-22s layout differs in %d lines; numbers %s' % (fn, ndiff, 'identical' if same else 'DIFFER'))
    if not same:
        ok = False
        for x, y in zip(A, B):
            if x != y:
                print('   orig: %r\n   new : %r' % (x, y))
                break
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)

"""Regression test: spher_main (upstream demo main) + library spher.f must reproduce the
reference output that Mishchenko ships at the top of spher.f, number by number.

  python graspkernels/tests/test_spher_upstream.py <scratch dir>
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, '..', 'spher_main')
tmp = sys.argv[1]
os.makedirs(tmp, exist_ok=True)
subprocess.run([EXE], cwd=tmp, check=True)          # writes spher.print and spher.write
num = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[DE][-+]?\d+)?')
toks = lambda fn: [float(t.replace('D', 'E')) for t in num.findall(open(fn).read())]
a = toks(os.path.join(tmp, 'spher.print'))
b = toks(os.path.join(HERE, 'spher_upstream_test_output.txt'))
bad = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
print('%d numbers compared, %d differ' % (len(b), len(bad) + abs(len(a) - len(b))))
ok = len(a) == len(b) and not bad
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)

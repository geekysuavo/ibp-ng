#!/usr/bin/env python

# enable imports from the ../../expts directory.
from subprocess import call
import os, stat, sys, time
sys.path.insert(0, '../../expts')

# import the utility and pdb modules.
from util import *
import pdb

# if supplied, get the dihedral expansion value.
eps = 0
if len(sys.argv) >= 2:
  eps = float(sys.argv[1])

# if supplied, get the distance restraint limit.
dmax = 10
if len(sys.argv) >= 3:
  dmax = float(sys.argv[2])

# build the output filenames.
fres = 'loops-{}-{}.res'.format(eps, dmax)
frun = 'loops-{}-{}.run'.format(eps, dmax)

# parse the target pdb file (from ../talos/talos.dcd).
s = pdb.pdb('loops.pdb')

# get the residue index extents.
(rmin, rmax) = minmax([atom['resSeq'] for atom in s.models[0]])

# initialize the lines of output.
lines = ['', '{{* loops.res: eps = {}, dmax = {} *}}'.format(eps, dmax)]

# add dihedral restraints.
for resid in range(rmin, rmax + 1):
  h = resid - rmin
  i = resid - rmin + 1
  j = resid - rmin + 2

  lines.append('')
  lines.append('{{* resid {} *}}'.format(i))

  Ch  = s.select(resid - 1, 'C')
  Ni  = s.select(resid,     'N')
  CAi = s.select(resid,     'CA')
  Ci  = s.select(resid,     'C')
  Nj  = s.select(resid + 1, 'N')
  CAj = s.select(resid + 1, 'CA')

  delta = 0.1
  if i in [11,12,13, 30,31,32,33,34, 46,47,48,49, 61,62,63,65]:
    delta = delta + eps

  if Ch and Ni and CAi and Ci:
    phi = pdb.dihed(Ch, Ni, CAi, Ci)[0]
    line = dihed.format(h, 'C', i, 'N', i, 'CA', i, 'C', phi, delta)
    lines.append(line)

  if Ni and CAi and Ci and Nj:
    psi = pdb.dihed(Ni, CAi, Ci, Nj)[0]
    line = dihed.format(i, 'N', i, 'CA', i, 'C', j, 'N', psi, delta)
    lines.append(line)

  if CAi and Ci and Nj and CAj:
    omega = pdb.dihed(CAi, Ci, Nj, CAj)[0]
    line = dihed.format(i, 'CA', i, 'C', j, 'N', j, 'CA', omega, 0.1)
    lines.append(line)

# add distance restraints.
numDist = 0
lines.append('')
lines.append('{* distance restraints *}')
for resid1 in range(rmin, rmax + 1):
  i1 = resid1 - rmin + 1
  H1 = s.select(resid1, 'H')
  if not H1:
    continue

  for resid2 in range(resid1 + 5, rmax + 1):
    i2 = resid2 - rmin + 1
    H2 = s.select(resid2, 'H')
    if not H2:
      continue

    d = pdb.dist(H1, H2)[0]
    if d <= dmax:
      numDist = numDist + 1
      line = dist.format(i1, 'H1', i2, 'H1', d, 0.25, 0.25)
      lines.append(line)

# write the restraints.
lines.append('')
f = open(fres, 'w')
f.write('\n'.join(lines))
f.close()

# build the runfile string.
runstr = ("""#!/bin/bash

../../bin/ibp-ng -v -v                     \\
  --input      loops.fa                    \\
  --psf        loops-{0}-{1}.psf           \\
  --dmdgp      loops-{0}-{1}.dmdgp         \\
  --output     loops-{0}-{1}.dcd           \\
  --format     dcd                         \\
  --restraints loops-{0}-{1}.res           \\
  --params     ../../lib/ibp-protein.par   \\
  --topology   ../../lib/ibp-protein.top   \\
  --reorder    ../../lib/ibp-protein.ord   \\
  --threads 1 --method dist,impr --limit 1 \\
  --branch-eps 0.01 --branch-max 16        \\
  --vdw-scale 0.5 --ddf-tol 0.1            \\
    1>loops-{0}-{1}.out 2>loops-{0}-{1}.err

""".format(eps, dmax))

# write the runfile.
f = open(frun, 'w')
f.write(runstr)
f.close()

# set the runfile as executable.
st = os.stat(frun)
os.chmod(frun, st.st_mode | stat.S_IEXEC)

# execute the runfile.
startTime = time.time()
call(['./loops-{}-{}.run'.format(eps, dmax)])
endTime = time.time()
timing = endTime - startTime

# output the final run results.
print('{} {} {} {}'.format(eps, dmax, numDist, timing))


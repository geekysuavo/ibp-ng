
# oneletter: conversion to single-letter codes.
# (switches all PRO -> ALA)
oneletter = {
  'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
  'TYR': 'Y', 'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T',
  'MET': 'M', 'ALA': 'A', 'GLY': 'G', 'PRO': 'A', 'CYS': 'C'}

# select: atom selection format substring.
select = '(resid {:>2d} and name {:>2s})'

# dist: distance restraint format string.
dist = ('assign ' + select + ' ' + select +
        ' {:>8.3f} {:>8.3f} {:>8.3f}')

# dihed: dihedral restraint format string.
dihed = ('assign ' + select + ' ' + select + '\n' +
         '       ' + select + ' ' + select +
         ' 1.0 {:>8.3f} {:>8.3f} 1')

# runstr: ibp-ng experiment script format string.
runstr = ("""#!/bin/bash

if [ "${1}" == "clean" ]; then
  rm -f output.*
  exit
fi

../../bin/ibp-ng -v -v        \\
  --input      input.fa       \\
  --psf        output.psf     \\
  --dmdgp      output.dmdgp   \\
  --output     output.dcd     \\
  --format     dcd            \\
  --restraints input.res      \\
  --params     input.par      \\
  --topology   ../protein.top \\
  --reorder    ../protein.ord \\
  --threads 1                 \\
  --method dist,impr          \\
  --limit 10000 --branch-eps 0.01 --branch-max 16 \\
  --vdw-scale 0.25 --ddf-tol 0.1

""")


# minmax: get the minimum and maximum values in a list.
#
def minmax(L):
  return (min(L), max(L))


# midrange: compute the center and range of the values in a list.
#
def midrange(L):
  (lmin, lmax) = minmax(L)
  return ((lmin + lmax) / 2, lmax - lmin)


# meanbounds: get the mean, upper and lower bounds of values in a list.
#
def meanbounds(L):
  (lmin, lmax) = minmax(L)
  lmean = sum(L) / float(len(L))
  return (lmean, lmean - lmin, lmax - lmean)


# fixomega: correct for circular shifts in the omega angle.
#
def fixomega(L):
  if type(L) is list:
    return [fixomega(l) for l in L]
  if L < -90:
    return L + 360
  return L


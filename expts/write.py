
# import the required modules.
import os, stat, pdb
from util import *

# runfile: write the run script string to a file.
#
def runfile(filename):
  # write the output file.
  f = open(filename, 'w')
  f.write(runstr)
  f.close()

  # set the file as executable.
  st = os.stat(filename)
  os.chmod(filename, st.st_mode | stat.S_IEXEC)


# fasta: write the sequence of a structure to a fasta file.
#
def fasta(struct, filename):
  # build a dictionary mapping residue number to name.
  smap = {atom['resSeq']: atom['resName'] for atom in struct.models[0]}
  seq = [oneletter[name] for name in smap.values()]

  # include the header.
  lines = ['> {}'.format(struct.filename)]

  # include the residues with regular line breaks.
  i = 0
  while i < len(seq):
    line = ''.join(seq[i : i + 50])
    lines.append(line)
    i = i + 50

  # include a couple blank lines.
  lines.append('')
  lines.append('')

  # build the output string and write it to file.
  f = open(filename, 'w')
  f.write('\n'.join(lines))
  f.close()


# restraints: write a restraint file from a pdb structure.
#
def restraints(struct, filename):
  f = open(filename, 'w')
  rbase(struct, f)
  rdihed(struct, f)
  rdist(struct, f)
  f.close()


# rbase: write a set of base restraints from a pdb structure.
#
def rbase(s, f):
  # get the start and end residue numbers.
  (rmin, rmax) = minmax([atom['resSeq'] for atom in s.models[0]])

  # include a heading.
  lines = ['', '{* base restraints *}']

  # loop over the residues.
  for resid in range(rmin, rmax + 1):
    # get the residue index.
    i = resid - rmin + 1

    # write the N-CA-C-HA dihedral restraint.
    line = dihed.format(i, 'N', i, 'CA', i, 'C', i, 'HA', -119, 0)
    lines.append(line)

  # write the CTER CA-O-C-O2 dihderal restraint.
  i = rmax - rmin + 1
  line = dihed.format(i, 'CA', i, 'O', i, 'C', i, 'O2', 171.5, 0)
  lines.append(line)

  # build the output string and write it to the file.
  lines.append('')
  f.write('\n'.join(lines))


# rdihed: write a set of dihedral restraints from a pdb structure.
#
def rdihed(s, f):
  # get the start and end residue numbers.
  (rmin, rmax) = minmax([atom['resSeq'] for atom in s.models[0]])

  # loop over the residues.
  lines = []
  for resid in range(rmin, rmax + 1):
    # get the residue indices.
    h = resid - rmin
    i = resid - rmin + 1
    j = resid - rmin + 2

    # print a heading.
    lines.append('')
    lines.append('{{* resid {} *}}'.format(i))

    # get the relevant atom positions.
    Ch  = s.select(resid - 1, 'C')
    Ni  = s.select(resid,     'N')
    CAi = s.select(resid,     'CA')
    Ci  = s.select(resid,     'C')
    Nj  = s.select(resid + 1, 'N')
    CAj = s.select(resid + 1, 'CA')

    # compute the phi dihedral.
    if Ch and Ni and CAi and Ci:
      phi = pdb.dihed(Ch, Ni, CAi, Ci)
      (phi0, dphi) = midrange(phi)
      line = dihed.format(h, 'C', i, 'N', i, 'CA', i, 'C', phi0, dphi)
      lines.append(line)

    # compute the psi dihedral.
    if Ni and CAi and Ci and Nj:
      psi = pdb.dihed(Ni, CAi, Ci, Nj)
      (psi0, dpsi) = midrange(psi)
      line = dihed.format(i, 'N', i, 'CA', i, 'C', j, 'N', psi0, dpsi)
      lines.append(line)

    # compute the omega dihedral.
    if CAi and Ci and Nj and CAj:
      omega = pdb.dihed(CAi, Ci, Nj, CAj)
      (omega0, domega) = midrange(fixomega(omega))
      line = dihed.format(i, 'CA', i, 'C', j, 'N', j, 'CA', omega0, domega)
      lines.append(line)

  # build the output string and write it to the file.
  lines.append('')
  f.write('\n'.join(lines))


# rdist: write a set of distance restraints from a pdb structure.
#
def rdist(s, f):
  # get the start and end residue numbers.
  (rmin, rmax) = minmax([atom['resSeq'] for atom in s.models[0]])

  # print a heading.
  lines = ['', '{* distance restraints *}']

  # loop over the residues.
  for resid1 in range(rmin, rmax + 1):
    # get the first residue index and amide proton.
    i1 = resid1 - rmin + 1
    H1 = s.select(resid1, 'H')
    if not H1:
      continue

    # loop again over the residues.
    for resid2 in range(resid1 + 5, rmax + 1):
      # get the second residue index and amide proton.
      i2 = resid2 - rmin + 1
      H2 = s.select(resid2, 'H')
      if not H2:
        continue

      # compute the distance.
      d = pdb.dist(H1, H2)
      if len([True for di in d if di <= 6]):
        (d0, da, db) = meanbounds(d)
        line = dist.format(i1, 'H1', i2, 'H1', d0, da, db)
        lines.append(line)

  # build the output string and write it to the file.
  lines.append('')
  f.write('\n'.join(lines))


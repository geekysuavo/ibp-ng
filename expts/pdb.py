
# import the required modules.
from util import median
import math

# pdb: class for parsing protein data bank files.
#
class pdb:
  # default values:
  ok = ['N', 'CA', 'C', 'O', 'H', 'HA']
  filename = 'input.pdb'
  models = []


  # __init__: construct a pdb object.
  #
  def __init__(self, filename = None):
    # if supplied, set the input filename.
    if filename:
      self.filename = filename

    # open the input file.
    f = open(self.filename)
    lines = [line.strip() for line in f.readlines()]
    nblines = zip(range(len(lines)), lines)
    f.close()

    # get the MODEL and ENDMDL line numbers.
    model = [line[0] for line in nblines if line[1].startswith('MODEL ')]
    endmdl = [line[0] for line in nblines if line[1].startswith('ENDMDL')]

    # get the list of ATOM blocks.
    ablks = []
    for i in range(len(model)):
      blk = lines[model[i] + 1 : endmdl[i]]
      blk = [line for line in blk if line.startswith('ATOM  ')]
      blk = [self.__atom(line) for line in blk]
      blk = [atom for atom in blk if self.__accept(atom)]
      ablks.append(blk)

    # store the parsed lists.
    self.models = ablks


  # select: get the coordinates of selected atoms in a pdb file.
  #
  def select(self, resid = 1, name = 'H'):
    # get the initial selection.
    sel = [[atom['pos'] for atom in mdl
            if atom['name'] == name
            and atom['resSeq'] == resid]
           for mdl in self.models]

    # check for empty results.
    sel = [atom[0] for atom in sel if len(atom)]
    if sel == []:
      return None

    # return the final selection.
    return sel


  # bond: get bond length statistics from a pdb file.
  #
  def bond(self, ni, nj, ij = 0):
    # get the residue index extents.
    rlist = [atom['resSeq'] for atom in self.models[0]]
    rmin = min(rlist)
    rmax = max(rlist)

    # build the list of distances.
    d = []
    for i in range(rmin, rmax + 1):
      j = i + ij
      ai = self.select(i, ni)
      aj = self.select(j, nj)
      if ai and aj:
        d = d + dist(ai, aj)

    # compute and return statistics on the list.
    return (median(d), min(d), max(d))


  # angle: get angle statistics from a pdb file.
  #
  def angle(self, ni, nj, nk, ij = 0, ik = 0):
    # get the residue index extents.
    rlist = [atom['resSeq'] for atom in self.models[0]]
    rmin = min(rlist)
    rmax = max(rlist)

    # build the list of angles.
    theta = []
    for i in range(rmin, rmax + 1):
      j = i + ij
      k = i + ik
      ai = self.select(i, ni)
      aj = self.select(j, nj)
      ak = self.select(k, nk)
      if ai and aj and ak:
        theta = theta + angle(ai, aj, ak)

    # compute and return statistics on the list.
    return (median(theta), min(theta), max(theta))


  # __atom: parse a single ATOM line in a pdb file.
  #
  def __atom(self, line):
    # create the field dictionary.
    fields = {}

    # parse the fields.
    fields['serial'] = int(line[6:11].strip())
    fields['name'] = line[12:16].strip()
    fields['altLoc'] = line[16]
    fields['resName'] = line[17:20].strip()
    fields['chainID'] = line[21]
    fields['resSeq'] = int(line[22:26].strip())
    fields['iCode'] = line[26]
    fields['x'] = float(line[30:38].strip())
    fields['y'] = float(line[38:46].strip())
    fields['z'] = float(line[46:54].strip())

    # finally, add a new field for convenience.
    fields['pos'] = (fields['x'], fields['y'], fields['z'])

    # return the final field dictionary.
    return fields


  # __accept: determine whether an atom should be accepted or not.
  #
  def __accept(self, atom):
    # check for acceptance.
    if atom['chainID'] == 'A' and atom['name'] in self.ok:
      return True
    else:
      return False

# ===

# sub3: subtract 3-vectors.
#
def sub3(a, b):
  return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


# dot3: compute 3-vector dot products.
#
def dot3(a, b):
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# cross3: compute 3-vector cross products.
#
def cross3(a, b):
  return (a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0])


# len3: compute the length of a 3-vector.
#
def len3(v):
  return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


# norm3: normalize a 3-vector to unit length.
#
def norm3(v):
  vlen = len3(v)
  return (v[0] / vlen, v[1] / vlen, v[2] / vlen)


# dist: compute the distance between two points.
#
def dist(a, b):
  # handle coordinate lists.
  if type(a) is list:
    return [dist(*z) for z in zip(a, b)]

  # compute and return the distance.
  return len3(sub3(a, b))


# angle: compute the angle formed by three points.
#
def angle(a, b, c):
  # handle coordinate lists.
  if type(a) is list:
    return [angle(*z) for z in zip(a, b, c)]

  # compute and return the angle.
  theta = math.acos(dot3(norm3(sub3(a, b)), norm3(sub3(c, b))))
  return math.degrees(theta)


# dihed: compute the dihedral formed by four points.
#
def dihed(a, b, c, d):
  # handle coordinate lists.
  if type(a) is list:
    return [dihed(*z) for z in zip(a, b, c, d)]

  # set up the three vector system.
  b1 = sub3(a, b)
  b2 = sub3(b, c)
  b3 = sub3(c, d)

  # compute the normal vectors.
  n1 = norm3(cross3(b1, b2))
  n2 = norm3(cross3(b2, b3))
  m = cross3(n1, norm3(b2))

  # compute and return the dihedral.
  omega = math.atan2(dot3(m, n2), dot3(n1, n2))
  return math.degrees(omega)


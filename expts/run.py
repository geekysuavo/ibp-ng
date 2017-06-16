
# import the required modules.
import os, pdb, write, urllib2

# ids: list of pdb identifiers to download and process for experiments.
ids = [
  '2d9g', '2k4x', '1u3o', '2db6', '2lgn', '2ee9', '1z66', '2ofq',
  '2fft', '1lv3', '1iw4', '2jrj', '1ri9', '2row', '2yrt', '1fyj',
  '2kkt', '1x47', '1je3', '1q8l', '2jr7', '1jjr', '1vdi', '2cpm',
  '1fsh', '2joi', '1cfa', '2keo', '2kbq', '2eap'
]

# loop over the identifiers.
for ident in ids:
  # skip existing directories.
  if os.path.exists(ident):
    continue

  # create the directory.
  os.mkdir(ident)
  uc = ident.upper()

  # download the pdb file.
  url = 'https://files.rcsb.org/download/' + uc + '.pdb'
  response = urllib2.urlopen(url)
  data = response.read()

  # write the result to disk.
  fpdb = os.path.join(ident, 'input.pdb')
  f = open(fpdb, 'w')
  f.write(data)
  f.close()

  # open and parse the pdb file.
  struct = pdb.pdb(fpdb)

  # write the required input files.
  write.runfile(os.path.join(ident, 'run'))
  write.fasta(struct, os.path.join(ident, 'input.fa'))
  write.restraints(struct, os.path.join(ident, 'input.res'))


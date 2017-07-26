
# import the required modules.
from struct import unpack

# dcd: class for parsing binary trajectories.
#
class dcd:
  # default values.
  filename = 'input.dcd'
  frames = []

  # __init__: construct a dcd object.
  #
  def __init__(self, filename = None):
    # if supplied, set the input filename.
    if filename:
      self.filename = filename

    # open the input file.
    f = open(self.filename, 'rb')
    data = f.read()
    f.close()

    # unpack the first data record.
    header = unpack('=i4s9id10i', data[0:92])

    # unpack the second data record.
    comment = unpack('=2i160si', data[92:264])[2]

    # unpack the third data record.
    atoms = unpack('=3i', data[264:276])[1]

    # determine the block format.
    block_format = '=' + 'i{}fi'.format(atoms) * 3
    block_size = 24 + 12 * atoms
    block_start = 276

    # loop to read each frame.
    while block_start + block_size - 1 <= len(data):
      # unpack the current data block.
      dat = data[block_start : block_start + block_size]
      block = unpack(block_format, dat)

      # extract the frame coordinates.
      x = block[1 : atoms + 1]
      y = block[atoms + 3 : 2 * atoms + 3]
      z = block[2 * atoms + 5 : 3 * atoms + 5]

      # store the frame coordinates.
      xyz = zip(x, y, z)
      self.frames.append(xyz)

      # move to the next block.
      block_start += block_size


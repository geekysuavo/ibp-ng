
# enable imports from the ../../expts directory.
sys.path.insert(0, '../../expts')

# import the numpy and dcd modules.
import numpy as np
import dcd

# define a header string for the output latex document.
header = '''
\\documentclass{article}
\\begin{document}

\\title{Results of ``re-folding'' the TALOS instance}
\\maketitle

These tables summarize the results of my analysis of how iBP-ng handles
flexibility of loops regions in protein structures. Procedure:
\\begin{enumerate}
 \\item{The TALOS instance was solved using its original restraints.}
 \\item{The dihedrals of the resulting solution were computed and used to
  generate new restraints, where all $(\\phi,\\psi)$ angles in flexible
  ``loop'' regions were expanded by the amount $\\Delta\\phi$,$\\Delta\\psi$.}
 \\item{All $H_N$--$H_N$ distances less than $d_{max}$ from the solution
  were also added to the new restraint set, yielding a set of $N_{dist}$
  distance restraints.\\footnote{A 0.5\\r{A} interval was added to each
  distance.}}
 \\item{The time required to obtain a solution using the new restraints
  was recorded. Where the time is marked as ``--'', iBP-ng ran for its
  maximum time of 5 hours without yielding a solution.}
\\end{enumerate}

\\begin{center}
'''

# define a format string for printing table rows.
out = ('  {0:2d} & {1:2d} & {2:4d} & {3:s} & {4:s} \\\\')

# fields(): format a line of statistics text as latex, with additional
# post-processing results included.
#
def fields(line):
  s = line.strip().split()

  eps = int(round(float(s[0])))
  dmax = int(round(float(s[1])))
  rest = int(s[2])
  done = (s[3] == 'True')
  time = float(s[4])

  if done:
    final = '{0:.1f}'.format(time)
    P = dcd.dcd('../talos/talos.dcd')
    Q = dcd.dcd('loops-{0}-{1}.dcd'.format(*s))
    rmsd = float(rmsd.strip())
    rmsd = '{0:.3f}'.format(rmsd)

  else:
    final = '--'
    rmsd = '--'

  return (eps, dmax, rest, final, rmsd)

# opener(): open a new table.
#
def opener():
  print('\\begin{tabular}{llllr}')
  print('  $\Delta\phi$,$\Delta\psi$ &')
  print('  $d_{max}$ (\\r{A}) &')
  print('  $N_{dist}$ &')
  print('  $t$ (s) &')
  print('  rmsd (\\r{A}) \\\\')

# closer(): close an open table.
#
def closer():
  print('\\end{tabular}')
  print


# the file 'stats.all' is the concatenation of all 'stats-*' files
# produced by running './run' in the current directory.
L = [fields(l) for l in open('stats.all').readlines()]
L = sorted(L, key = lambda s: 1000 * s[0] + s[1])

# begin the latex document.
print(header)

# loop over every row produced by processing the 'stats.all' file.
for i in range(len(L)):
  if i > 0 and i % 27 == 0:
    closer()

  if i % 27 == 0:
    opener()

  if i % 9 == 0:
    print(' \\hline')

  print(out.format(*L[i]))

  if i == len(L) - 1:
    closer()

# end the latex document.
print
print('\\end{center}')
print('\\end{document}')
print


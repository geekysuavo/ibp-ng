
# import the necessary modules.
import numpy as np
import math

# rmsd: compute the *aligned* rmsd between two structures.
#
def rmsd(X, Y):
  D = len(X[0])
  N = len(X)

  x0 = sum(X) / N
  y0 = sum(Y) / N

  Xc = []
  Yc = []

  for i in range(N):
    Xc.append(X[i] - x0)
    Yc.append(Y[i] - y0)

  P = np.array(Xc)
  Q = np.array(Yc)

  C = np.dot(P.T, Q)
  (V, S, W) = np.linalg.svd(C)

  d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
  if d:
    S[-1] = -S[-1]
    V[:,-1] = -V[:,-1]

  U = np.dot(V, W)
  P = np.dot(P, U)

  dev = 0.0
  for p, q in zip(P, Q):
    dev += sum([(p[i] - q[i])**2.0 for i in range(D)])

  return np.sqrt(dev / N)


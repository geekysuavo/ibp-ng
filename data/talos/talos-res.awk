/^REMARK/ {
  next
}

NF == 11 && $NF != "None" {
  I = $1
  H = I - 1
  J = I + 1
  resname = $2

  phi = $3
  psi = $4

  dphi = $5
  dpsi = $6

  printf("{* %s%s: phi, psi, omega *}\n", resname, I)

  printf("assign (resid %2d and name C)  (resid %2d and name N)\n", H, I)
  printf("       (resid %2d and name CA) (resid %2d and name C)", I, I)
  printf("  1.0 %8.3f %8.3f 1\n", phi, dphi)

  printf("assign (resid %2d and name N)  (resid %2d and name CA)\n", I, I)
  printf("       (resid %2d and name C)  (resid %2d and name N)", I, J)
  printf("  1.0 %8.3f %8.3f 1\n", psi, dpsi)

  printf("assign (resid %2d and name CA) (resid %2d and name C)\n", I, I)
  printf("       (resid %2d and name N)  (resid %2d and name CA)", J, J)
  printf(" 1.0 %8.3f %8.3f 1\n", 180.0, 0.1)

  printf("\n")
}

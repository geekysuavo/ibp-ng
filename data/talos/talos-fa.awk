BEGIN {
  w = 0

  printf("> talos.fa\n")
}

/^REMARK/ {
  next
}

NF == 11 {
  resname = $2
  w++

  printf("%s", resname)
  if (w >= 50) {
    printf("\n")
    w = 0
  }
}

END {
  printf("\n")
}

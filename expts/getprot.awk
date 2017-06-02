/^ATOM/ {
  name = substr($0, 13, 4);
  chn = substr($0, 22, 1);
  gsub(/ /, "", name);

  if (name != "H" || chn != "A") { next; }

  resid = substr($0, 23, 4);
  x = substr($0, 31, 8);
  y = substr($0, 39, 8);
  z = substr($0, 47, 8);

  print resid, x, y, z;
}

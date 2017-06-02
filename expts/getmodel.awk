BEGIN { p = 0 }
$0 ~ /^MODEL/ && $2 == id { p = 1 }
{
  if (p) {
    print $0
  }
}
$0 ~ /^ENDMDL/ { p = 0 }

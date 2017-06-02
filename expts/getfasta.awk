/^CHN/ { print ">", $2 }
/^SEQ/ { gsub(/P/, "A", $3); print $3 }

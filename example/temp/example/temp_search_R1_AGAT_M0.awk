match($0, /AGAT/) {s=$0; m=0; while((n=match(s, /AGAT/))>0){m+=n; printf "%s,", m; m+=3; s=substr(s, n+4)}}

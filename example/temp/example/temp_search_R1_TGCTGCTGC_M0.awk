match($0, /TGCTGCTGC/) {s=$0; m=0; while((n=match(s, /TGCTGCTGC/))>0){m+=n; printf "%s,", m; m+=8; s=substr(s, n+9)}}

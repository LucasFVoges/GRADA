match($0, /ACTTCTGGACT/) {s=$0; m=0; while((n=match(s, /ACTTCTGGACT/))>0){m+=n; printf "%s,", m; m+=10; s=substr(s, n+11)}}

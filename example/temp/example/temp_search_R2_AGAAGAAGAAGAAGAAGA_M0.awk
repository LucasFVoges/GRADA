match($0, /AGAAGAAGAAGAAGAAGA/) {s=$0; m=0; while((n=match(s, /AGAAGAAGAAGAAGAAGA/))>0){m+=n; printf "%s,", m; m+=17; s=substr(s, n+18)}}

match($0, /(.{0,1}|A.)GAAGAAGAAGAAGAAGA|A(.{0,1}|G.)AAGAAGAAGAAGAAGA|AG(.{0,1}|A.)AGAAGAAGAAGAAGA|AGA(.{0,1}|A.)GAAGAAGAAGAAGA|AGAA(.{0,1}|G.)AAGAAGAAGAAGA|AGAAG(.{0,1}|A.)AGAAGAAGAAGA|AGAAGA(.{0,1}|A.)GAAGAAGAAGA|AGAAGAA(.{0,1}|G.)AAGAAGAAGA|AGAAGAAG(.{0,1}|A.)AGAAGAAGA|AGAAGAAGA(.{0,1}|A.)GAAGAAGA|AGAAGAAGAA(.{0,1}|G.)AAGAAGA|AGAAGAAGAAG(.{0,1}|A.)AGAAGA|AGAAGAAGAAGA(.{0,1}|A.)GAAGA|AGAAGAAGAAGAA(.{0,1}|G.)AAGA|AGAAGAAGAAGAAG(.{0,1}|A.)AGA|AGAAGAAGAAGAAGA(.{0,1}|A.)GA|AGAAGAAGAAGAAGAA(.{0,1}|G.)A|AGAAGAAGAAGAAGAAG.{0,1}/) {s=$0; m=0; while((n=match(s, /(.{0,1}|A.)GAAGAAGAAGAAGAAGA|A(.{0,1}|G.)AAGAAGAAGAAGAAGA|AG(.{0,1}|A.)AGAAGAAGAAGAAGA|AGA(.{0,1}|A.)GAAGAAGAAGAAGA|AGAA(.{0,1}|G.)AAGAAGAAGAAGA|AGAAG(.{0,1}|A.)AGAAGAAGAAGA|AGAAGA(.{0,1}|A.)GAAGAAGAAGA|AGAAGAA(.{0,1}|G.)AAGAAGAAGA|AGAAGAAG(.{0,1}|A.)AGAAGAAGA|AGAAGAAGA(.{0,1}|A.)GAAGAAGA|AGAAGAAGAA(.{0,1}|G.)AAGAAGA|AGAAGAAGAAG(.{0,1}|A.)AGAAGA|AGAAGAAGAAGA(.{0,1}|A.)GAAGA|AGAAGAAGAAGAA(.{0,1}|G.)AAGA|AGAAGAAGAAGAAG(.{0,1}|A.)AGA|AGAAGAAGAAGAAGA(.{0,1}|A.)GA|AGAAGAAGAAGAAGAA(.{0,1}|G.)A|AGAAGAAGAAGAAGAAG.{0,1}/))>0){m+=n; printf "%s,", m; m+=17; s=substr(s, n+18)}}
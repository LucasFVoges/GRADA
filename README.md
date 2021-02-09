# GRADA
simple GRep ADapter Analyser

> R-Script utilazing the unix bash powers for adapter (sequence) analysis in a read file.


> Many programms like fastp, fastQC, PRINSEQ are great for analyzing and preprocessing NGS read files. Though they are relativ complex programms and I wanted to see on a easy to understand level the contamination of a specific sequence (eg the adapter) in a read (fastq) file. This is possible by grep / agrep / wc commands. This scipt allows to get an overview of the contamination of a sequence in your read files.

System requirements:

- UNIX System (developed on Linux Mint 20.1 Cinnamon)
- R-Studio
  - packages:
  - DT
  - parallel

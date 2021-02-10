# GRADA - R-Package
simple GRep ADapter Analyser

> R-Script utilazing the unix bash powers for adapter (sequence) analysis in a read file.

> Many programms like fastp, fastQC, PRINSEQ are great for analyzing and preprocessing NGS read files. Though they are relativ complex programms and I wanted to see on a easy to understand level the contamination of a specific sequence (eg the adapter) in a read (fastq) file. This is possible by grep / agrep / wc commands. This scipt allows to get an overview of the contamination of a sequence in your read files.

## System requirements:

- UNIX System (developed on Linux Mint 20.1 Cinnamon)
- R-Studio
  - packages are suggested:
  - DT
  - parallel
  
## Installation:
just run:
`devtools::install_github("LucasFVoges/GRADA", build_vignettes = TRUE)`

## BEFORE USING THE SCRIPT
Please note, that GRADA will create an temp/ folder in your working directory. It will save the results here but also the .txt files wich will have the corresponding reads inside.

These files can be very big and can be deleted afterwards!

## Usage:
```
library(GRADA)

# recommended at the moment:
library(parallel) 
library(DT)
```

You can load some example data:
```
read1 <- system.file("extdata", "grada_R1.fastq", package = "GRADA")
read2 <- system.file("extdata", "grada_R2.fastq", package = "GRADA")
seq <- system.file("extdata", "adapter_list.txt", package = "GRADA")
```

Then you can call the table and plot function (this will change a little in future releases)
```
grada_table(seq = seq, read1 = read1, read2 = read2)
grada_plot()
```
There are additional options to these functions.

## Example:
GRADA comes with an example vignette and example data (very basic).

`browseVignettes("GRADA")`

or: 

```
library(GRADA)   
vignette("example")
```

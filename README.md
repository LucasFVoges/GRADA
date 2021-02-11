# GRADA - R-Package
simple GRep ADapter Analyser

> R-Script utilazing the unix bash powers for adapter (sequence) analysis in a read file.

> Many programms like fastp, fastQC, PRINSEQ are great for analyzing and preprocessing NGS read files. Though they are relativ complex programms and I wanted to see on a easy to understand level the contamination of a specific sequence (eg the adapter) in a read (fastq) file. This is possible by grep / agrep / wc commands. This scipt allows to get an overview of the contamination of a sequence in your read files.

A complete introduction: [see in Vignettes](http://htmlpreview.github.io/?https://github.com/LucasFVoges/GRADA/blob/master/vignettes/introduction.html)

## System requirements:

- UNIX System (developed on Linux Mint 20.1 Cinnamon)
- R-Studio
- R-packages are suggested:
  - DT
  - parallel
  - knitr
  - rmarkdown
  
## Installation:
just run:  
`devtools::install_github("LucasFVoges/GRADA", build_vignettes = TRUE)`
or for latest development version:  
`devtools::install_github("LucasFVoges/GRADA", branch = "dev", build_vignettes = TRUE)`

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

Then you can call the analyze functions. The table and plot function will render the results. For rendering only the "grada_table.txt" and "adapter_positions.Rdata" are needed.
```
grada_analyze(PE = TRUE, seq = seq, read1 = read1, read2 = read2)
grada_analyze_positions(PE = TRUE, readlength = 150)
grada_table()
grada_plot()
```
There are additional options to these functions.

### Table

For the table there is:
```
# For a kable-table:
grada_table_simple()
# For a rmarkdown-table (requires "rmarkdown" package):
grada_table_md()
# For a DT interactive table (requires "DT" package):
grada_table_DT() = grada_table()
```
But you could use your own table-script. you can load the data with: `load("temp/Adapter_Positions.Rdata")`

### Plot

For the plots there is:
```
# For a standard barplot:
grada_plot_bar() = grada_plot()
```

## Example:
GRADA comes with an example vignette and example data (very basic).

[See in Vignettes](http://htmlpreview.github.io/?https://github.com/LucasFVoges/GRADA/blob/master/vignettes/example.html)

`browseVignettes("GRADA")`

or: 

```
library(GRADA)   
vignette("example")
```

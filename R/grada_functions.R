#' GRADA - analyze
#'
#' This function will perform a unix "agrep" and "wc" to look through the read.fastq files. This will be iterativly done for each mismatch allowed (note: mismatches are only other characters - if the adapter sequence is overlaping it will not be found). The Mismatches can't be bigger then the shortest sequence!
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE)
#' @param seq sequences to search for (adapters) (A text file containing:">Name
#' Sequence
#' IlumniaUniversalAdapter
#' AGATCGGAAGAGC")
#' @param read1 Path to R1 read file (std. NULL)
#' @param read2 Path to R2 read file (if paired data, std. NULL)
#' @param M_min minimal mismatches allowed (std. 0)
#' @param M_max maximal mismatches allowed (std. 2)
#' @param output the folder where all data will be created. (std. "temp/")
#' @param numCores Number of cores to use. If numCores=1 then the normal lapply function is used and the parallel package is not neccessary! (std. detectCores()/2)
#' @return A Table as .txt of the found sequences (adapters) and .txt files containig reads per sequence and mistake.
#' @export
grada_analyze <- function(PE=TRUE, seq=NULL, read1=NULL, read2=NULL, M_min=0, M_max=2, output= "temp/", numCores=detectCores()/2){
  ####### Testing function calling #####
  writeLines(paste0("#### GRADA v.1.2 ####\nfun: grada_table()\nPaired data: ", PE, "\nR1: ", read1, ", R2: ", read2,"\nSequences: ", seq, "\nFrom ", M_min, " to ", M_max, " mismatches\n\nSave to: ", output, "\n#####################\n"))
  if (is.null(read1)){stop("No data file (read1 = ?)")}
  if (M_min < 0){stop("You have entered a negative number...")}
  if (M_min > M_max){stop("minimum mismatches can't be higher then maximum!")}
  if (PE == TRUE && is.null(read2) == TRUE){stop("No paired data supplied?")}
  if (M_max - M_min > 5){print("WARNING: that will take time and space!  many mismatches...")}

  ####### USER INPUT ##########

  # Input for the Adapters
  input1 <- file(seq, "r")
  # Output folder:
  if (!(dir.exists(output))){dir.create(output)}

  #############################

  ##### INIT #####
  adapters <- matrix(ncol = 3, nrow = 0)
  colnames(adapters) <- c("Adapter", "Sequence", "Length")

  ### Read the adapter file ###
  while (TRUE) {
    line = trimws(readLines(input1, n = 1))       # read one line of the file
    if (length(line) == 0 ) {                   # End of file
      break
    }
    if (substr(line, 1, 1) == ">") {            # get the Adapter Name
      adap_name <- substring(line, 2)
      next
    } else if (line == "") {                    # Skip any Blank Lines
      next
    }
    if (adap_name == ""){
      print("Could not read complete adapter file. There is a formating problem. Please check your adapters!")
      break
    }
    adapters <- rbind(adapters, c(adap_name, line, nchar(line)))
    adap_name <- ""
  }
  close(input1)

  ##### Functions #####

  # writing reads :
  find_adap_read <- function(Rnum, adapter, Mnum){
    if (Mnum == M_max){
      if (Rnum == 1){read <- read1} else {read <- read2}
    } else {
      read <- paste0(output, "temp_R", Rnum, "_", adapter, "_M", Mnum+1, ".txt")
    }
    # This is the AGREP Unix function for finding the sequences:
    system(paste0("agrep -", Mnum, " " , adapter, " ", read, " > ", output, "temp_R", Rnum, "_", adapter, "_M", Mnum, ".txt"), intern = FALSE, wait = TRUE)
    # This is the counting Unix wc function:
    system(paste0("wc -l ", output, "temp_R", Rnum, "_", adapter,"_M", Mnum, ".txt | cut -f1 -d' ' > ", output, "counts_temp_R", Rnum, "_", adapter,"_M", Mnum, ".txt"), intern = FALSE, wait = TRUE)
  }

  write_adap_read <- function(adapter){
    if (PE){Rnum <- 1:2} else {Rnum <- 1:1}                  # check if paired data is TRUE.
    for (M in M_max:M_min) {
      if (numCores == 1){
        lapply(Rnum, find_adap_read, adapter, M)  
      } else {
        mclapply(Rnum, find_adap_read, adapter, M, mc.silent = TRUE, mc.cores = tail(Rnum, n=1))        # and give to the find function
      }
    }
    return(paste0(adapter, " seems to has finished"))
  }

  ##### writing temp adapter_files #####
  if (numCores == 1){
    lapply(adapters[,"Sequence"], write_adap_read)
    } else {
    mclapply(adapters[,"Sequence"], write_adap_read, mc.silent = TRUE, mc.cores = numCores)
    }
  ##### saving the count files to the list #####
  if (PE) {Rnum <- 2} else {Rnum <- 1}
  for (R in 1:Rnum) {
    for (M in M_max:M_min){
      count_col <- c()
      for (adapter in adapters[,"Sequence"]) {
        counts <- read.delim(paste0(output,"counts_temp_R", R, "_", adapter,"_M", M, ".txt"), header = FALSE)
        count_col <- c(count_col, as.integer(counts))
      }
      adapters <- cbind(adapters, count_col)
      colnames(adapters)[colnames(adapters) == "count_col"] <- paste0("R",R,"M",M)
    }
  }

  ##### Results - writing table #####
  write.table(adapters, file = paste0(output, "grada_table.txt"), row.names = FALSE)
  writeLines("GRADA has made a table!\n")

  #### Delete count files ####
  system(paste0("rm ", output, "counts_*"), intern = FALSE)
}


#' GRADA - Analyze Positions
#'
#' This function will count the 1. positions of the sequences found in grada_table. Therefore grada_table needs to be run prior. It will only search for mismatches = 0.
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE)
#' @param readlength longest read! for X-axis. (std. 150)
#' @param input input folder where the "grada_table.txt" is (output from grada_table())
#' @param numCores Number of cores to use. If numCores=1 then the normal lapply function is used and the parallel package is not neccessary! (std. detectCores()/2)
#' @return Rdata matrix and can be used for your own plots as well.
#' @export
grada_analyze_positions <- function(PE = TRUE, readlength = 150, input="temp/", numCores=detectCores()/2){
  missM <- 0 # No Effect until now... (unix awk command has to change)
  # Size of the rads must be set here!
  writeLines(paste0("#### GRADA v.1.2 ####\nfun: grada_plot()\nPaired data: ", PE, "\nRead length: ", readlength, "\nInput: ", input, "grada_table.txt\nMismatches", missM, " (fixed)\n#####################\n"))
  adapter_positions <- matrix(ncol = readlength, nrow = 0)
  colnames(adapter_positions) <- c(1:readlength)
  if (!(file.exists(paste0(input, "grada_table.txt")))){stop("could not find grada_table.txt - pleas run grada_table() first.")}
  adapters <- read.table(paste0(input, "grada_table.txt"), header = TRUE)
  if (!(file.exists(paste0(input, "temp_R1_", adapters[1,"Sequence"], "_M0.txt")))){stop("File with 0 mismatches not available!?")}
  if (PE){
    if (!(file.exists(paste0(input, "temp_R2_", adapters[1,"Sequence"], "_M0.txt")))){stop("File with 0 mismatches not available for R2 reads!?")}
  }

  # Function:
  find_positions <- function(adapter){
    if (PE) {Rnum <- 1:2} else {Rnum <- 1:1}
    if (numCores == 1){
      lapply(Rnum, find_positions_sec, adapter)
    } else {
      mclapply(Rnum, find_positions_sec, adapter, mc.cores = tail(Rnum, n=1))
    }  
    return(paste0(adapter, " has been prepared for plotting"))
  }

  find_positions_sec <- function(R, adapter){
    poslist <- c(1:(readlength))
    # was before: poslist <- c(0:(readlength-1)) This is wrong. And was used wrongly! :(
    awk_positions <- system(intern = TRUE, paste0("awk -v p=", adapter," 'index($0,p) {s=$0; m=0; while((n=index(s,p))>0){m+=n; printf \"%s,\", m; s=substr(s, n+1)}}' ", input, "temp_R", R, "_", adapter, "_M", missM,".txt"))
    # build vector
    positions <- as.vector(as.numeric(unlist(strsplit(awk_positions, ","))))
    poslist <- append(poslist, positions)
    # change points to frequency over the readlength nt.
    poscount <- table(poslist)
    # correction loop due to first creation of poslist.
    for (i in 1:readlength) {
      poscount[i] <- poscount[i] -1
    }
    save(poscount, file = paste0(input, paste0("Adapter_Positions_R", R, "_", adapter, "_M", missM, ".Rdata")))
  }

  # call function
  if (numCores == 1){
    lapply(adapters[,"Sequence"], find_positions)
  } else {  
    mclapply(adapters[,"Sequence"], find_positions, mc.silent = TRUE,  mc.cores = numCores)
  }
  
  # load all data in a table
  for (adapter_i in adapters[,"Sequence"]) {
    if (PE){Rnum <- 1:2} else {Rnum <- 1:1}
    for (R in Rnum){
      load(paste0(input, paste0("Adapter_Positions_R", R ,"_", adapter_i, "_M", missM, ".Rdata")))
      adapter_positions <- rbind(adapter_positions, poscount)
      # Change Rowname to Adapter
      row.names(adapter_positions)[nrow(adapter_positions)] <- paste0(as.character(adapters[, "Adapter"][adapters[,"Sequence"] == adapter_i]), " R_", R)
    }
  }

  # Save Data File for plotting.
  save(adapter_positions, file = paste0(input, "Adapter_Positions.Rdata"))

  #### Delete temp files ####
  system(paste0("rm ", input, "Adapter_Positions_*"), intern = FALSE)

 writeLines(paste0("GRADA has analyzed the positions!\n"))
}

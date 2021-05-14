#' GRADA - analyze
#'
#' This function will perform a unix "agrep" and "wc" (and "zcat") to look through the read.fastq(.gz) files.
#' It will take more time if used on zipped (.gz) files. Consider to unzip...
#' This search will be iteratively done for each mismatch allowed. The number of mismatches can't be equal or higher then the shortest adapter sequence!
#' 
#' ! This works only case sensitive !
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE)
#' @param seq "sequence text file" to search for (adapters) (A text file containing 2 lines per Sequence: 1:>AdapterName 2:Sequence for example: 1:>IlumniaUniversalAdapter 2:AGATCGGAAGAGC") - look in examples for an, well, example.
#' @param read1 "R1 read file.fastq(.gz)" (std. NULL)
#' @param read2 "R2 read file.fastq(.gz)" (if paired data, std. NULL)
#' @param M_min minimal mismatches allowed (std. 0)
#' @param M_max maximal mismatches allowed (std. 2)
#' @param output the folder where all data will be created. (std. "temp/")
#' @param omitMeta if TRUE: only sequences will be analyzed and the rest of the read data omitted. This can be beneficial if sequences are found in the title line. Although fewer lines need to be searched, this option will be slower. (std. FALSE)
#' @param keepfiles if FALSE: will delete temp files after counting. Note that analyze_positions is depending on these files! But this could save disk space in case of big files or adapter lists. (std. TRUE)
#' @param numCores Number of cores to use. If numCores=1 then the normal lapply function is used and the parallel package is not necessary! (std. detectCores()%/%2)
#' @return A Table as .txt of the found sequences (adapters) and .txt files containing reads per sequence and mistake.
#' @export
grada_analyze <- function(PE=TRUE, seq=NULL, read1=NULL, read2=NULL, M_min=0, M_max=2, output= "temp/", omitMeta=FALSE, keepfiles=TRUE, numCores=detectCores()%/%2){
  ####### Testing function calling #####
  writeLines(paste0("#### GRADA v.1.2 ####\nfun: grada_analyze()\nPaired data: ", PE, "\nR1: ", read1, ", R2: ", read2,"\nSequences: ", seq, "\nFrom ", M_min, " to ", M_max, " mismatches\n\nSave to: ", output, "\n#####################\n"))
  # Test if Unix commands are available:
  if(!length(system("which agrep", intern = TRUE)) >= 1){stop("X - 'agrep' seems to be missing on system...")}
  if(!length(system("which wc", intern = TRUE)) >= 1){stop("X - 'wc' seems to be missing on system...")}
  if(!length(system("which zcat", intern = TRUE)) >= 1 && omitMeta == TRUE){stop("X - 'zcat' seems to be missing on system...")}
  # Test if input is acceptable:
  if (is.null(read1)){stop("X - No data file (read1 = ?)")}
  if (M_min < 0){stop("X - You have entered a negative number...")}
  if (M_max > 8){
    M_max <- 8
    writeLines("! - M_max was corrected to 8, which is the maximum supported by agrep.\n")
  }
  if (M_min > M_max){stop("X - Minimum mismatches can't be higher then maximum!")}
 
  if (PE == TRUE && is.null(read2) == TRUE){stop("X - No paired data supplied?")}
  if (M_max - M_min > 5){print("! - This will take time and space! (many mismatches...)")}

  ####### USER INPUT ##########

  # Input for the Adapters
  input1 <- file(seq, "r")
  # Output folder (create if not exist):
  if (!(dir.exists(output))){
    dir.create(output, recursive = TRUE)
  }

  #############################
  
  ### CHECK NAMES and INPUT###
  read1 <- trimws(read1)
  if(PE){read2 <- trimws(read2)}
  
  if(numCores == 0){numCores <- 1}
  
  ##### INIT #####
  adapters <- matrix(ncol = 3, nrow = 0)
  colnames(adapters) <- c("Adapter", "Sequence", "Length")

  ### Read the adapter file ###
  search_name <- TRUE
  while (TRUE) {
    line = trimws(readLines(input1, n = 1))     # read one line of the file
    if (length(line) == 0 ) {                   # End of file
      break
    }
    if (substr(line, 1, 1) == ">") {            # get the Adapter Name
      adap_name <- trimws(line, whitespace = "[>,<, ,!]")
      next
    } else if (line == "") {                    # Skip any Blank Lines
      next
    }
    if (line %in% adapters[,"Sequence"]){
      writeLines(paste0("! - A doubled sequence: ", line, " is skipped...!\n"))
      next
    } 
    if (nchar(line) <= M_max){
      writeLines("! - An adapter is shorter or equal the M_max mismatches. This adapter will be skipped!\n")
      next
    }
    if (adap_name == ""){
      writeLines("! - Could not read complete adapter file. There is a formating problem. Please check your adapters!\n")
      break
    }
    # Check for too small sequenzes

    adapters <- rbind(adapters, c(adap_name, line, nchar(line)))
    adap_name <- ""
  }
  close(input1)

  ##### Functions #####

  # writing reads :
  find_adap_read <- function(Rnum, adapter, Mnum){
    ### This part handles if the read file or the temp files are used:
    if (Mnum == M_max){
      if (Rnum == 1){read <- read1} else {read <- read2}
      # check if file is compressed:
      if(substr(read, nchar(read)-2, nchar(read)) == ".gz"){
        compressed <- TRUE
      } else {
        compressed <- FALSE
      }
      # check if omitMeta:
      if(omitMeta){
        omit <- TRUE
      } else{
        omit <- FALSE
      }
    } else {
      read <- paste0(output, "temp_R", Rnum, "_", adapter, "_M", Mnum+1, ".txt")
      compressed <- FALSE
      omit <- FALSE
    }
      
    ### the AGREP Unix function ###
    if(!compressed && !omit){
      system(paste0("agrep -", Mnum, " " , adapter, " '", read, "' > ", output, "temp_R", Rnum, "_", adapter, "_M", Mnum, ".txt"), intern = FALSE, wait = TRUE)
      } else {
        ### Make the omit/compressed string:
        if(compressed && omit){
          pipestringA <- paste0("zcat '", read, "' | sed '/^$/d' | awk 'NR % 4 == 2' | ")
        } else if (compressed && !omit){
          pipestringA <- paste0("zcat '", read, "' | ")
        } else if (!compressed && omit){
          pipestringA <- paste0("sed '/^$/d' '", read, "' | awk 'NR % 4 == 2' | ")
        }
        system(paste0(pipestringA, "agrep -", Mnum, " " , adapter, " > ", output, "temp_R", Rnum, "_", adapter, "_M", Mnum, ".txt"), intern = FALSE, wait = TRUE)
      }
    
    # This is the unix wc counting function:
    system(paste0("wc -l ", output, "temp_R", Rnum, "_", adapter,"_M", Mnum, ".txt | cut -f1 -d' ' > ", output, "counts_temp_R", Rnum, "_", adapter,"_M", Mnum, ".txt"), intern = FALSE, wait = TRUE)
    # Remove orig temp_ file if keepfiles=FALSE
    if (!keepfiles) {
      # Delete previous file?
      if (!Mnum == M_max) {
        Mnum_del <- Mnum + 1
        system(paste0("rm ", output, "temp_R", Rnum, "_", adapter, "_M", Mnum_del, ".txt"), intern = FALSE)
      }
      # Delete actual file?
      if (Mnum == M_min) {
        system(paste0("rm ", output, "temp_R", Rnum, "_", adapter, "_M", Mnum, ".txt"), intern = FALSE)
      }
    }
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
#' This function will count the 1. positions of the sequences found in grada_table. Therefore grada_analyze() needs to be run prior, with appropriate mismatches. Although up to 8 mismatches can be analyzed, depending on the size and findings this could take a long time. GRADA uses the utils aregexec() function to find the results. Another possible way is to search with unix "awk". This takes a regular expression (which can be awful long for many mismatches) and therefore will be killed (out of memory?) by the system sometimes!
#' 
#' !method="awk" will search for multiple results in one read (line) method="R" only finds the first occurrence per read (line)
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE)
#' @param readlength longest read! for x-axis. (std. 150)
#' @param input input folder where the "grada_table.txt" is (output from grada_table())
#' @param M_min minimal mismatches (std. 0) 
#' @param M_max maximal mismatches (std. 2)
#' @param keepfiles if FALSE: delete the temp files from grada_analyze() (std. FALSE)
#' @param method different implementation of position finding. "R" (1) uses R utils package. "awk" (2) uses unix awk command. (std. 1)
#' @param numCores Number of cores to use. If numCores=1 then the normal lapply function is used and the parallel package is not necessary! (std. detectCores()%/%2)
#' @return Rdata matrix and can be used for your own plots as well.
#' @export
grada_analyze_positions <- function(PE = TRUE, readlength = 150, input="temp/", M_min=0, M_max=2, keepfiles=FALSE, method=1, numCores=detectCores()%/%2){
  #### Size of the rads must be set here!
  writeLines(paste0("#### GRADA v.1.2 ####\nfun: grada_analyze_positions()\nPaired data: ", PE, "\nRead length: ", readlength, "\nInput: ", input, "grada_table.txt\nMismatches ", M_min, " to ", M_min,  "\n#####################\n"))
  

  
  #### CHECK INPUT ####
  if(numCores == 0){numCores <- 1}
  if (M_min < 0){stop("You have entered a negative number...")}
  if (M_min > M_max){stop("minimum mismatches can't be higher then maximum!")}
  if (M_max - M_min > 5){print("WARNING: this will take time and space! (many mismatches...)")}
  
  if (!(file.exists(paste0(input, "grada_table.txt")))){stop("could not find grada_table.txt - pleas run grada_table() first.")}
  adapters <- read.table(paste0(input, "grada_table.txt"), header = TRUE)
  
  for (MMs in M_max:M_min) {
    if (!(file.exists(paste0(input, "temp_R1_", adapters[1,"Sequence"], "_M", MMs ,".txt")))){
      stop(paste0("R1-File with ", MMs, " mismatches not available!?"))
    }
    if (PE){
      if (!(file.exists(paste0(input, "temp_R2_", adapters[1,"Sequence"], "_M", MMs ,".txt")))){
        stop(paste0("R2-File with ", MMs, " mismatches not available!?"))
      }
    }
  }
  
  # Method:
  if (method == 1 || method == "R" || method == "r"){
    search_method <- 1
  } else if (method == 2 || method == "awk" || method == "AWK"){
    search_method <- 2
  } else {  
    writeLines("! - Method not detectable. Using 'R' instead.\n###\n")
    search_method <- 1
  }
  
  #### Test if Unix commands are available:
  if (search_method == 2){
    if(!length(system("which awk", intern = TRUE)) >= 1){
      stop("awk seems to be missing on system...")
    }
  }
  
  # Function:
  find_positions <- function(adapter, missM){
    if (PE) {Rnum <- 1:2} else {Rnum <- 1:1}
    if (search_method == 2){
      if (numCores == 1){
        lapply(Rnum, find_positions_awk, adapter, missM)
      } else {
        mclapply(Rnum, find_positions_awk, adapter, missM, mc.cores = tail(Rnum, n=1))
      }  
    } else {  
      if (numCores == 1){
        lapply(Rnum, find_positions_R, adapter, missM)
      } else {
        mclapply(Rnum, find_positions_R, adapter, missM, mc.cores = tail(Rnum, n=1))
      }  
    }  
    return(paste0(adapter, " has been prepared for plotting"))
  }
  
  ### AWK VERSION of POSITIONS with regex!
  find_positions_awk <- function(R, adapter, missM){
    poslist <- c(1:(readlength))                  # "1:" because awk works with "1"
    adapter <- as.character(adapter)
    if (missM > 0) {
      regex_string <- ""
      length_adapter <- nchar(adapter)
      char_adapter <- strsplit(adapter, "")[[1]]
      # This will look in the pascal triangle for max number of possible combinations:
      # Old: pt_possibils <- lapply(length_adapter, function(i) choose(i, missM))[[1]]
      pt_possibils <- choose(length_adapter, missM)

      # Here some info for the regex
      writeLines(paste0("Regex for: ", missM, " mismatches. Possibilities: ", pt_possibils))
      # the combinations are available via combn()
      pt_combinations <- combn(1:length_adapter, missM, simplify = TRUE)    # error if missM higher than length of adapter!
      # generate the adapter regex:
      for (possibles in 1:pt_possibils) {
        this_char_adapter <- char_adapter
        for (this_pos in pt_combinations[,possibles]){
          if (this_pos == length_adapter) {
            # because insertion at the end is pointless.
            this_char_adapter[this_pos] <- ".{0,1}"
          } else {
            # This allown Indels as well. Without it would be just "."
            this_char_adapter[this_pos] <- paste0("(.{0,1}|", char_adapter[this_pos], ".)")
          }
        }
        regex_string <- paste0(regex_string, paste0(this_char_adapter, collapse = ""), "|")
        if (possibles == pt_possibils) {          # so in the end is no "|"
          regex_string <- substr(regex_string, 1, nchar(regex_string)-1)
        }
      }
    } else if (missM == 0){
      regex_string <- adapter
    } else {
      stop("This should never happen. Error in nr of mismatches. it is not 0 and not higher than 0...?")
    }
    
    ### AWK search for positions: (a ".awk" is needed because bash does not allow that long expressions.)
    awk_string <- paste0(input, "temp_search_R", R, "_", adapter, "_M", missM, ".awk")
    write(x = paste0("match($0, /", regex_string, "/) {s=$0; m=0; while((n=match(s, /", regex_string, "/))>0){m+=n; printf \"%s,\", m; m+=", nchar(adapter)-1, "; s=substr(s, n+", nchar(adapter),")}}"), file = awk_string)
    awk_positions <- system(intern = TRUE, paste0("awk -f ", awk_string, " ", input, "temp_R", R, "_", adapter, "_M", missM,".txt"))
    
    # BACKUP: awk_positions <- system(intern = TRUE, paste0("awk 'match($0,", regex, ") {s=$0; m=0; while((n=match(s,", regex, "))>0){m+=n; printf \"%s,\", m; m+=", nchar(adapter)-1, "; s=substr(s, n+", nchar(adapter),")}}' ", input, "temp_R", R, "_", adapter, "_M", missM,".txt"))
    
    # BACKUP OLD: awk_positions <- system(intern = TRUE, paste0("awk -v p=", adapter," 'index($0,p) {s=$0; m=0; while((n=index(s,p))>0){m+=n; printf \"%s,\", m; s=substr(s, n+1)}}' ", input, "temp_R", R, "_", adapter, "_M", missM,".txt"))
    
    ###DOCUMENTATION of the awk Command:
    # builds a string: "1,5,..." which are the positions inside the lines!
    # awk -v p=ACTTCTGGACT    # -v allows the variable   | This is not neccesary as seen in example 2!
    # 'index($0,p)            # per line with findings it will do the following:
    #  s=$0; m=0;             # start param
    #  while(                 # loop printing all results
    #  (n=index(s,p))>0){     # n is the index wherep is found.
    #    m+=n;                # sets m to this position (+= because remaining string starts again at 1)
    #    printf "%s,", m;     # the printing of m
    #    m+= nchar(adapter)-1 # sets M adequate to remaining string.
    #    s=substr(s, n+nchar(adapter))     # set s to remaining string.
    #

    ### build vector
    positions <- as.vector(as.numeric(unlist(strsplit(awk_positions, ","))))
    poslist <- append(poslist, positions)

    ### change points to frequency over the readlength nt.
    poscount <- table(poslist)

    ### correction loop due to first creation of poslist.
    for (i in 1:readlength) {
      poscount[i] <- poscount[i] -1
    }
    save(poscount, file = paste0(input, paste0("Adapter_Positions_R", R, "_", adapter, "_M", missM, ".Rdata")))
  }

  
  # # NEW AGREP VERSION    # PROBLEM agrep behaves a little strange... see issues
  # # 
  # # This is abandoned!
  # #
  # find_positions_agrep <- function(R, adapter){
  #   poslist <- c(1:(readlength))    
  #   adapter <- as.character(adapter)
  #   length_adapter <- nchar(adapter)
  #   this_file <- paste0(input, "temp_R", R, "_", adapter, "_M", missM,".txt")
  #   
  #   # OLD: per line search:
  #   # agrep_position <- c()
  #   # lines_of_file <- as.integer(system(intern=TRUE, paste0("wc -l ", input, "temp_R", R, "_", adapter, "_M", missM,".txt | cut -f1 -d' '")))
  #   # if (lines_of_file > 0) {
  #   #   # The approach with agrep per line:
  #   #   for (readline in 1:lines_of_file) {
  #   #     agrep_position <- c(agrep_position, as.integer(system(intern = TRUE, paste0("sed -n ", readline, "p ", input, "temp_R", R, "_", adapter, "_M", missM,".txt | agrep -", missM, " -b ", adapter, " | cut -f1 -d'='"))) - (length_adapter-1))
  #   #   }
  #   #     
  #   # } else { # end lines_of_file > 0
  #   #   agrep_position <- c()
  #   # }
  #   
  #   # New approach, as the line by line is too slow! PROBLEM: agrep seems to run in bug -.-
  #   # if (file.size(this_file) != 0){
  #   #   # make agrep vector for the positions  # PROBLEM: it will find more than one per line:. how?
  #   #   agrep_adapter_positions <- as.integer(system(intern = TRUE, paste0("agrep -", missM, " -b '", adapter, "' ", this_file,  " | cut -f1 -d'='")))
  #   #   agrep_newline_positions <- c(-1, as.integer(system(intern = TRUE, paste0("agrep -b '$' ", this_file, " | cut -f1 -d'='"))))
  #   #   # The last endofline is not interesting. But we need to account for the first line, which is done in the command above.
  #   #   agrep_newline_positions <- agrep_newline_positions[-length(agrep_newline_positions)]
  #   #   if (length(agrep_adapter_positions) != length(agrep_newline_positions)) {
  #   #     stop("There is a problem that should not happen. Are the settings identical to grada_analyze()? Are the files untempered?")
  #   #   }
  #   #   # Get the real position (starting at 1!):
  #   #   agrep_adapter_positions <- agrep_adapter_positions - agrep_newline_positions - (length_adapter - 1)
  #   #   
  #   #   ### add vector to poslist
  #   #   poslist <- append(poslist, agrep_adapter_positions)
  #   # }  
  #   
  #   ### change points to frequency over the readlength nt.
  #   poscount <- table(poslist)
  #   
  #   ### correction loop due to first creation of poslist.
  #   for (i in 1:readlength) {
  #     poscount[i] <- poscount[i] - 1
  #   }
  #   save(poscount, file = paste0(input, paste0("Adapter_Positions_R", R, "_", adapter, "_M", missM, ".Rdata")))
  # }
  
  # As the files are smaller now, I think it is a possible solution to use it in R or python.
  find_positions_R <- function(R, adapter, missM){
    poslist <- c(1:(readlength))    
    adapter <- as.character(adapter)
    length_adapter <- nchar(adapter)
    this_file <- paste0(input, "temp_R", R, "_", adapter, "_M", missM,".txt")
    open_file <- file(this_file, "r")
    content <- readLines(open_file)
    close(open_file)
    findings <- unlist(aregexec(adapter, content, fixed = TRUE, max.distance = missM))
    # add vector to poslist:
    poslist <- append(poslist, findings)
    ### change points to frequency over the readlength nt.
    poscount <- table(poslist)
    ### correction loop due to first creation of poslist.
    for (i in 1:readlength) {
      poscount[i] <- poscount[i] - 1
    }
    save(poscount, file = paste0(input, paste0("Adapter_Positions_R", R, "_", adapter, "_M", missM, ".Rdata")))
  }  

  # call function   # This will make 1 file for each mismatch!
  for (m in M_min:M_max) {
    adapter_positions <- matrix(ncol = readlength, nrow = 0)
    colnames(adapter_positions) <- c(1:readlength)
    if (numCores == 1){
      lapply(adapters[,"Sequence"], find_positions, m)
    } else {  
      mclapply(adapters[,"Sequence"], find_positions, m, mc.silent = TRUE, mc.cores = numCores)
    }
    # load all data in a table
    for (adapter_i in adapters[,"Sequence"]) {
      if (PE){Rnum <- 1:2} else {Rnum <- 1:1}
      for (R in Rnum){
        load(paste0(input, paste0("Adapter_Positions_R", R ,"_", adapter_i, "_M", m, ".Rdata")))
        adapter_positions <- rbind(adapter_positions, poscount[1:readlength])
        # Change Rowname to Adapter
        row.names(adapter_positions)[nrow(adapter_positions)] <- paste0(as.character(adapters[, "Adapter"][adapters[,"Sequence"] == adapter_i]), " R_", R)
      }
    }
    # Save Data File for plotting.
    save(adapter_positions, file = paste0(input, "Adapter_Positions_M", m, ".Rdata"))
  }
  #### Delete temp files ####
  system(paste0("rm ", input, "Adapter_Positions_R*"), intern = FALSE)
  if (!keepfiles) {
    system(paste0("rm ", input, "temp_R*"), intern = FALSE)
  }

 writeLines(paste0("GRADA has analyzed the positions!\n"))
}

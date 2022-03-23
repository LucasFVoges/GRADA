#' GRADA - Python - RC Adapters 
#'
#' This function will reverse complement all adapters
#' 
#' Will work with DNA bases only (A,T,C,G,N)
#' 
#' @param file Path to the adapter file
#' @export
grada_py_make_reverse_complement <- function(file){
  path <- paste(system.file(package="GRADA"), "make_reverse_complement.py", sep="/")
  command <- paste("python", path, "-f", file, sep = " ")
  try(suppressWarnings(response <- system(command, intern=TRUE)), silent = TRUE)
  
  # if(!is.null(attr(response,"status"))){
  #   if(attr(response,"status") == 1){
  #     response <- ""
  #     cat("To Field Empty") } } 
  # 
  
}

#' GRADA - Python - average read length
#'
#' This function will return the average read length for a standard .fastq file as a .txt file.
#' 
#' The first 100 reads will be evaluated (if not set to another parameter)
#' 
#' @param file Path to the fastq file
#' @param sample number of reads to look at (std. 100)
#' @export
grada_py_AVG_LengthOfReads <- function(file, sample = 100){
  path <- paste(system.file(package="GRADA"), "AVG_LengthOfReads.py", sep="/")
  command <- paste("python", path, "-f", file, "-nb", sample, sep = " ")
  try(suppressWarnings(response <- system(command, intern=TRUE)), silent = TRUE)
}
#' GRADA - Python
#'
#' This function will reverse complement all adapters (if they have standard bases)
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
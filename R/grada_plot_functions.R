#' GRADA - Plot
#'
#' This function will refer to the stdandard plot function. 
#' 
#' It will only accept standard settings. You cannot define anything here. Please use another _plot function, if you want to change the in-/output.
#' 
#' Standard plot function: grada_plot_bar()
#'
#' @return R barplots of the adapter positions.
#' @export
grada_plot <- function(){
  grada_plot_bar()
  }  

#' GRADA - Barplot
#'
#' This function will plot the results of grada_analyze_positions(). The results are stored in "Adapter_Positions.Rdata"
#' 
#' The first position of the found Sequence will be plottet! If the read is 100 bases long and the adapter 10, the maximal position will therefor be 91!
#'
#' plots will be skipped if no adapters are found. If there are any difficulties with this function, may set skip to FALSE.
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE)
#' @param input input folder where the "dapter_Positions.Rdata" is (output from grada_table() std. "temp/")
#' @param skip if TRUE, it will skip plots for empty Data (if adapter is not found) so you will not get empty plots (sometimes it happens that the first plost is present anyway. this is due to a matrix R-command). (std. TRUE)
#' @param plot_row is the par(mfrow=c(plot_row,plot_col)) for arrenging the plots (std. 2)
#' @param plot_col is the par(mfrow=c(plot_row,plot_col)) for arrenging the plots (std. 2)
#' @return R barplots of the adapter positions.
#' @export
grada_plot_bar <- function(PE = TRUE, input="temp/", skip=TRUE, plot_row=2, plot_col=2){
  missM <- 0 # No Effect until now... (unix awk command has to change)
  writeLines(paste0("#### GRADA v.1.0 ####\nfun: grada_plot()\nPaired data: ", PE, "\nInput: ", input, "Adapter_Positions.Rdata\nSkip empty plots: ", skip, "\nMismatches", missM, " (fixed)\n#####################\n"))

  #### Load Data ####
  if (substr(input, nchar(input) - 1 + 1, nchar(input)) == "/"){
    load(paste0(input, "Adapter_Positions.Rdata"))  
  } else {
    stop("input folder problem. (Please add '/' to the end)")
  }
    
  #### PLOT ####
  # delete the data with no content!
  if (skip){
    if (PE){
      kills <- c()
      for (i in row.names(adapter_positions)){
        a <- which(row.names(adapter_positions) == i)
        if (a %% 2 == 0) {
          next
        }
        if (sum(adapter_positions[a,]) == 0 && sum(adapter_positions[a+1,]) == 0){
          kills <- c(kills, a, a+1)
        }
      }
    } else {
      kills <- c()
      for (i in row.names(adapter_positions)){
        a <- which(row.names(adapter_positions) == i)
        if (sum(adapter_positions[a,]) == 0){
          kills <- c(kills, a)
        }
      }
    }

    # follow line is just for fixing a bug. I need to fix!
    while ((nrow(adapter_positions) - length(kills)) < 2){kills <- head(kills, -1)}
    adapter_positions <- adapter_positions[-kills,]
  }
  # Graph generation: (barplot simple one after another.)
  par(mfrow=c(plot_row,plot_col))
  for (i in row.names(adapter_positions)){
    barplot(adapter_positions[i,],
            col = "#1b98e0",
            main = sprintf("%s (1. pos!)", i),
            #ylim = c(0, 90000),
            las = 2,
            cex.names = .8,
            cex.main = .8,
            border = NA)
  }
 writeLines(paste0("GRADA has plotted in " , plot_row, "x", plot_col, " arrangement!\n"))
}

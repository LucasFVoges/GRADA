#' GRADA - Plot
#'
#' This function will refer to the standard plot function. 
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
#' This function will plot the results of grada_analyze_positions() which are stored in "Adapter_Positions.Rdata"
#' 
#' standard: The first position of the found Sequence will be plotted! If the read is 100 bases long and the adapter 10, the maximal position will therefore be 91!
#' full_length = TRUE will set this plot to use all adapter positions. In the above mentioned example: 91,92,93,94,95,96,97,98,99,100 will be counted. Be aware, that this will not take into account any indel mismatches.
#'
#' plots will be skipped if no adapters are found. If there are any difficulties with this function, may set skip to FALSE.
#'
#' @param PE paired data? TRUE / FALSE (std. TRUE) (obsolete option. is no longer used)
#' @param input input folder where the "adapter_Positions.Rdata" is (output from grada_analyze_positions() std. "temp/")
#' @param M_min min mismatches. (std. 0)
#' @param M_max max mismatches. (std. 2)
#' @param full_length if TRUE it will plot all positions of the adapter (instead of the first pos only) (std. FALSE)
#' @param skip if TRUE, it will skip plots for empty Data, so you will not get empty plots. (std. FALSE)
#' @param plot_row is the par(mfrow=c(plot_row,plot_col)) for arrenging the plots (std. 2)
#' @param plot_col is the par(mfrow=c(plot_row,plot_col)) for arrenging the plots (std. 2)
#' @param colour colour scheme for plot. (1 = std. / 2 = differen / 3 = black and white) or you can use your own in the format: colour = c("red", "blue", #FFFFFF ...)
#' @return R barplots of the adapter positions.
#' @export
grada_plot_bar <- function(PE = TRUE, input="temp/", M_min=0, M_max=2, full_length=FALSE, skip=FALSE, plot_row=2, plot_col=2, colour=1){
  #### CHECK INPUT ####
  if (M_min < 0){stop("You have entered a negative number...")}
  if (M_min > M_max){stop("minimum mismatches can't be higher then maximum!")}

  # doubled variable if needed to change inside code.
  missM <- M_max
  
  ### Colour Scheme
  if (colour[1] == 1) {
  color_chem <- c("#000929", "#000e42", "#00186e", "#05269c", "#0d33b8", "#143cc9", "#234bdb", "#3761fa", "#577bff")
  } else if (colour[1] == 2) {
  color_chem <- c("#E8985E", "#A9714B", "#54442B", "#262A10", "#141204", "#000000", "#000000", "#000000", "#000000")
  } else if(colour[1] == 3) {
    color_chem <- c("#400100", "#960200", "#CE6C47", "#FFD046", "#EADAA2", "#F4EBCD", "#FFF5D6", "#FFEBEB", "#FBF2EF")
  } else if (colour[1] == 4) {
    color_chem <- c("#000000", "#1B1B1B", "#303030", "#474747", "#5E5E5E", "#777777", "#919191", "#ABABAB", "#C6C6C6")
  } else {
  color_chem <- colour
  }
    
  #### Load Data ####
  plotlist <- list()
  if (substr(input, nchar(input) - 1 + 1, nchar(input)) == "/"){
    for (MMs in M_max:M_min) {
      load(paste0(input, "Adapter_Positions_M", MMs, ".Rdata"))
      plotlist[[MMs+1]] <- adapter_positions
    }
    adapter_content <- read.table(paste0(input, "grada_table.txt"), header = TRUE, row.names = 1)
  } else {
    stop("input folder problem. (Please add '/' to the end)")
  }
  
  
  #### Full Plot advance
  if (full_length){
    # This screams for lapply. ... 
    for (adapter_name in row.names(plotlist[[missM+1]])){
      adapter <- substr(adapter_name, 1, nchar(adapter_name)-4)
      adapter_length <- adapter_content[adapter,"Length"]
      # with one adapter go over all tables:
      for (MMs in M_max:M_min) {
        for(zpos in ncol(plotlist[[MMs+1]]):1){
          if (plotlist[[MMs+1]][adapter_name,zpos] > 0){
            for(j in 1:adapter_length-1){
              if (!zpos+j > ncol(plotlist[[MMs+1]])){ #adapter problematc with indels
                plotlist[[MMs+1]][adapter_name,zpos+j] <- plotlist[[MMs+1]][adapter_name,zpos+j] + plotlist[[MMs+1]][adapter_name,zpos]
              }
            }  
          }  
        }
      }  
    }  
  }  
  
  
    
  #### CLEAR UP ####        OUTDATED - New Data handling and other method used. (Tough this one had R1/R2 skip dependent output)
  # delete the data with no content! 
  # if (skip == TRUE && M_max == M_min){
  #   if (PE){
  #     kills <- c()
  #     for (i in row.names(adapter_positions)){
  #       a <- which(row.names(adapter_positions) == i)
  #       if (a %% 2 == 0) {
  #         next
  #       }
  #       if (sum(adapter_positions[a,]) == 0 && sum(adapter_positions[a+1,]) == 0){
  #         kills <- c(kills, a, a+1)
  #       }
  #     }
  #   } else {
  #     kills <- c()
  #     for (i in row.names(adapter_positions)){
  #       a <- which(row.names(adapter_positions) == i)
  #       if (sum(adapter_positions[a,]) == 0){
  #         kills <- c(kills, a)
  #       }
  #     }
  #   }
  # 
  #   # follow line is just for fixing a bug, so that at leas one plot will remain.
  #   while ((nrow(adapter_positions) - length(kills)) < 2){kills <- head(kills, -1)}
  #   if(!is.null(kills)){
  #     adapter_positions <- adapter_positions[-kills,]
  #   }
  # }
  
  #### PLOT ####
  # Graph generation: (barplot simple one after another.)
  par(mfrow=c(plot_row,plot_col))
  lab <- c(1, rep(NA, (ncol(plotlist[[missM+1]])-1)))
  lab[seq(10, ncol(plotlist[[missM+1]]), 10)] <- seq(10, ncol(plotlist[[missM+1]]), 10)
  for (i in row.names(plotlist[[missM+1]])){
    plotdata <- c()
    sum <- 0
    thiscolor_chem <- color_chem # color for plotting
    for (MMs in M_max:M_min) {
      plotdata <- rbind(plotdata, plotlist[[MMs+1]][i,])
      thissum <- sum(plotlist[[MMs+1]][i,]) 
      # This will change the color for plotting, if first mismatch is empty plot... (if not used, the colors and legend will get desynchronised)
      if(thissum <= 0){
        thiscolor_chem <- thiscolor_chem[-1]
      }
      sum <- sum + thissum
    }
    # barplot will be stacked, so the number need to be subtracted
    if (nrow(plotdata) > 1) {
      for(row in 1:(nrow(plotdata)-1)){
        plotdata[row, ] <- plotdata[row, ] - plotdata[row+1, ]
        # ! - Problem can be, that the mismatch is not on the same position anymore. (indel) So negative number is possible.
        plotdata[row, ] <- replace(plotdata[row, ], plotdata[row, ] < 0, 0)
      }
    }
    # for skipping empty plots. (this should consider R1/R2 later)
    if (skip && sum <= 0) {
     next
    }
    limit <- max(plotdata)
    
    barplot(plotdata[nrow(plotdata):1,],
            col = thiscolor_chem,
            main = sprintf("%s (1.pos)", i),
            xlab="read position",
            ylab="counts", 
            axis.lty = 1,
            names.arg=lab,
            ylim = c(0, limit),
            space = 0,
            las = 2,
            cex.names = .8,
            cex.axis = .8,
            cex.main = .8,
            border = NA)
    legend("topright",
           legend = c(M_min:M_max), 
           lwd=4, 
           bty = "n", 
           col = color_chem,
           title = "Mismatches:",
           cex = .8)
  }
}


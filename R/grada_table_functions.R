#' GRADA - table
#'
#' This function will refer to the stdandard table function. 
#' 
#' It will only accept standard settings. You cannot define anything here. Please use another _table function, if you want to change the in-/output.
#' 
#' Standard table function: grada_table_DT()
#' 
#' @return a table of the found adapters (forward R1 and reverse R2) with mismatches (M)
#' @export
grada_table <- function(){
  grada_table_DT()
  }

#' GRADA - simple table
#'
#' This function will render a table of the results from grada_analyze()
#'
#' @param input the folder where "grada_table.txt" is saved. (std. "temp/")
#' @return a table of the found adapters (forward R1 and reverse R2) with mismatches (M)
#' @export
grada_table_simple <- function(input= "temp/"){
  # nothing here jet
  }


#' GRADA - Datatable (DT) 
#'
#' requires: library("DT")!
#' 
#' This function will render a table of the results from grada_analyze()
#' 
#' @param input the folder where "grada_table.txt" is saved. (std. "temp/")
#' @return a table of the found adapters (forward R1 and reverse R2) with mismatches (M)
#' @export
grada_table_DT <- function(input= "temp/"){
  ####### Testing function calling #####
  writeLines(paste0("#### GRADA v.1.1 ####\nfun: grada_table_DT()\n#####################\n"))

  #### return DT table ####
  adapter_content <- read.table(paste0(input, "grada_table.txt"), header = TRUE)
  content_l <- length(adapter_content)-1
  datatable(adapter_content, rownames = FALSE,
            caption = "Table 1: Sequence content in the read files",
            class = 'cell-border stripe',
            extensions = c('Buttons', 'FixedColumns'),
            options = list(dom = 'Bfrtip',
                           pageLength = 20,
                           lengthMenu = c(5, 10, 100),
                           buttons = list('copy', 'excel', 'pdf', 'pageLength', list(extend = 'colvis', columns = c(1:content_l))),
                           autoWidth = FALSE,
                           # seems not to effect anything:
                           # columnDefs = list(list(width = '100px', targets = c(1,2))),
                           fixedColumns = TRUE
            )) %>%  formatRound(c(3:content_l+1), 0)
}

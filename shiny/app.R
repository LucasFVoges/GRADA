library(shiny)
library(shinythemes)
library(data.table)
library(GRADA)

# Define UI for application
ui <- fluidPage(theme = shinytheme("superhero"),
    navbarPage(
    "GRADA v.1.2",
    tabPanel("Analyze Data",
             sidebarPanel(
                 tags$h3("Input:"),
                 textInput("txt1", "Given Name:", ""),
                 textInput("txt2", "Surname", ""),

                 
             ), # SidebarPanel
             mainPanel(
                 h1("GRADA"),
                 h3("Read Files"),
                 fileInput("file_R1", accept = ".fastq", label = h6("Read 1 file input")),
                 fileInput("file_R2", accept = ".fastq", label = h6("Read 2 file input")),
                 actionButton("analyzebutton", "GRADA analyze", class = "btn btn-primary"),
                 h3("Summary"),
                 verbatimTextOutput("txtout"),
                 
             ) # Mainpanel
             
    ), # tabPanel
    
    tabPanel("Results", 
            sidebarPanel(
                 tags$h3("Plot options:"),
                 sliderInput(inputId = "bins",
                             label = "Area of reads:",
                             min = 0,
                             max = 150,
                             value = c(0,150),
                             step = 5)
             ), # SidebarPanel
             
            mainPanel(
                 h1("GRADA"),
                 h3("Plots of Sequence - 1. Position"),
                 "Oh boy, how can I plot multiple times? I do not know how many plots there will be?",
                 plotOutput(outputId = "barplot"),
                 h3("Plots of Sequence - complete"),
                 "Okay, the plot is working, but to get interactive results, another plot function of grada is neccessary. Do not just copy the Code? ... well maybe its okay",
                 plotOutput(outputId = "barplot_full")
             ) # mainPanel
    ), # tabPanel
    
    
    tabPanel("About", "This page will contain some information")
    
    ) # navbarPage
) # fluidPage                


# Define server logic
server <- function(input, output) {
    
    output$txtout <- renderText({
        paste(input$txt1, input$txt2, sep = " ")
    })
    
    # Plots:
    output$barplot <- renderPlot({
        load(file = "Adapter_Positions.Rdata")
        plot_data <- adapter_positions[1,]
        min_data <- input$bins[1]
        max_data <- input$bins[2]
        if (min_data == max_data){max_data <- max_data + 1}
       
        barplot(plot_data[min_data:max_data], col = "#1b98e0", 
             border = NA,
             space = 0,
             las = 2,
             cex.names = .8,
             xlab = "read nt",
             ylab = "frequency",
             main = paste0(row.names(adapter_positions)[1]," - 1. Position"))          
    })
    output$barplot_full <- renderPlot({
        grada_plot_bar_full(PE = FALSE, input = "../grada/")
    })

} # server

# Run the application 
shinyApp(ui = ui, server = server)

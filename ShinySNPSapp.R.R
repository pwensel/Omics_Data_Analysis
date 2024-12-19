#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



library(shiny)
library(data.table) ##read txt files


SNPS<- data.frame(fread("SNPdataset.txt")) ##your path to the file
str(SNPS)



## INDICATIONS
## fill the gaps in the code below, by following the example ("04_mpg") from the link https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/
## the main idea is to plot the distribution of the age of individuals comparing by sex and case/controls

##1st. step
##change characters to factors

##2nd. step
##similar to the example, fill the gaps

##you can also try with other ideas



# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("SNPS dataset"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for variable to plot against mpg ----
      selectInput("variable", "Variable:",
                  c("XXXXXXX" = "XXXXXXX",
                    "XXXXXXX" = "XXXXXX")),
      
      # Input: Checkbox for whether outliers should be included ----
      checkboxInput("outliers", "Show outliers", TRUE)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      h3(textOutput("caption")),
      
      # Output: Plot of the requested variable against age ----
      plotOutput("Age_Plot")
      
    )
  )
)

# Define server logic to plot various variables against age ----
server <- function(input, output) {
  
  # Compute the formula text ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$Age_Plot functions
  
  formulaText <- reactive({
    paste("XXXXXXXX", input$variable) ### fill
  })
  
  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })
  
  # Generate a plot of the requested variable against mpg ----
  # and only exclude outliers if requested
  
  
  output$#####___fill___### <- renderPlot({
    boxplot(as.formula(formulaText()),
            data = SNPS,
            outline = input$outliers,
            col = "#75AADB", pch = 19)
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)
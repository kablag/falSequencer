library(shiny)
library(plotly)

shinyUI(fluidPage(

  # Application title
  titlePanel("Old Faithful Geyser Data"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      width = 2,
      numericInput("peakWidth", "Peak Width", 40, 5, 100),
      actionButton("generateChromatogram", "Generate")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      textInput("sequence", "Sequence", 
                "atagcccgatcyatatrggatgchagtmtancatgg",
                placeholder = "Input IUPAC sequence",
                width = "100%"),
      plotlyOutput("chromatogramPlot")
    )
  )
))

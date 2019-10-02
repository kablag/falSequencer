library(shiny)
library(plotly)

shinyUI(fluidPage(

  # Application title
  titlePanel("falSequencer"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      width = 2,
      sliderInput("peakWidth", "Points Per Peak", 40, 5, 100),
      sliderInput("subSeqFreq", "Subpeaks Frequency", 0, 1, 0.15),
      sliderInput("hScale", "Peak Width Compression", 0, 1, 0.3),
      sliderInput("vScale", "Min Peak Height", 0.5, 1, 0.7),
      sliderInput("subSeqVScale", "Max Subpeak Height", 0.1, 0.5, 0.2),
      sliderInput("backgroundVScale", "Max Background Height", 0, 0.05, 0.005, step = 0.001),
      actionButton("generateChromatogram", "Generate Chromatogram")
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

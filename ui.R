library(shiny)
# library(plotly)

shinyUI(fluidPage(

  # Application title
  titlePanel("falSequencer"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      width = 2,
      numericInput("samplesPerNuc", "Samples per nucleotide",
                   20, 5),
      numericInput("ppositionVariety", "Peak position variety",
                   0.1, 0),
      numericInput("pwidth", "Peak width (fraction of Samples per nucleotide)",
                   0.1, 0.01),
      numericInput("pwidthVariety", "Peak width variety",
                   0.1, 0),
      numericInput("pheight", "Peak height",
                   500, 10),
      numericInput("pheightVariety", "Peak height variety",
                   0.1, 0),
      numericInput("subpProb", "Subpeak probability",
                   0.1, 0),
      numericInput("subpHeight", "Subpeak height (fraction of Peak height)",
                   0.3, 0),
      numericInput("noisepProb", "Noise peak probability",
                   0.1, 0),
      numericInput("noisepHeight", "Noise peak height (fraction of Peak height)",
                   0.1, 0),
      textInput("comments", "Comments", "")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      textInput("sequence", "Sequence", 
                "GGCTGGAGAAGCTGCATCGCTCACCMGGGGCTGGTGGTCACTTTTTGATCTATTTCAGTGCATTGCTAAGGAACTGATTCCAGAAGCCACT",
                placeholder = "Input IUPAC sequence",
                width = "100%"),
      textInput("subsequence", "Optional subsequence", 
                "",
                placeholder = "Input IUPAC sequence",
                width = "100%"),
      fluidRow(
        column(2, actionButton("generateChromatogram", "Generate Chromatogram")),
        column(2, downloadButton("downloadSCF", "Download chromatogram as SCF"))
        ),
      plotOutput("chromatogramPlot")
    )
  )
))

library(shiny)
# library(scales)
# library(tidyverse)

source("generic.R")

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  scf <- eventReactive(
    input$generateChromatogram,
    {
      # isolate({
      req(input$sequence)
      subsequence <- if (input$subsequence == "") NULL else input$subsequence
      create.scf(gsub(" ", "", input$sequence),
                 samplesPerNuc = input$samplesPerNuc,
                 ppositionVariety = input$ppositionVariety,
                 pwidth = input$pwidth, pwidthVariety = input$pwidthVariety,
                 pheight = input$pheight, pheightVariety = input$pheightVariety,
                 subpSeq = subsequence, 
                 subpProb = input$subpProb, subpHeight = input$subpHeight,
                 noisepProb = input$noisepProb, noisepHeight = input$noisepHeight,
                 comments = input$comments)
  })
  
  output$chromatogramPlot <- renderPlot({
    req(scf())
    chromatogram(sangerseq(scf()))
  })

})

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
  
  output$chromatogramPlot <- renderPlot(
    width = 1200,
    # width = function() {
    #   input$nucsInRow * input$nucWidthPixels
    # },
    height = function() {
      height <- ceiling(scf()@header@bases / 80) * 120
      if (height < 200)
        height <- 200
      height
    },
    {
    req(scf())
    chromatogram(sangerseq(scf()))
  })
  
  output$downloadSCF <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), ".scf", sep = "")
    },
    content = function(file) {
      write.scf(scf(), file)
    }
  )

})

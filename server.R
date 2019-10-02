library(shiny)
library(scales)
library(tidyverse)

source("generic.R")

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  sequence <- reactive({
    req(input$sequence)
    strsplit(input$sequence %>% tolower(), NULL)[[1]]
  })
  
  chromatogram <- eventReactive(
    input$generateChromatogram,
    {
      # isolate({
      req(sequence())
      
      nucleotides <- list("a" = c("a", "r", "m", "w", "d", "h", "v", "n"),
                          "t" = c("t", "y", "k", "w", "b", "d", "h", "n"),
                          "g" = c("g", "r", "k", "s", "b", "d", "v", "n"),
                          "c" = c("c", "y", "m", "s", "b", "h", "v", "n"))
      
      
      genCurveElement <- function(currentNucleotide, nucleotide,
                                  subSeqFreq = 0.15,
                                  baseNucleotideWidth,
                                  hScale = 0.3,
                                  vScale = 0.5,
                                  subSeqVScale = 0.2,
                                  backgroundVScale = 0.005) {
        width <- round(baseNucleotideWidth * runif(1, 1 - hScale, 1) * 0.8)
        # createCurve <- (currentNucleotide %in% nucleotides[[nucleotide]]) ||
        #   (runif(1) < subSeqFreq)
        subSeqNum <- runif(1)
        subSeq <- subSeqNum <= subSeqFreq && subSeqNum > backgroundVScale 
        height <- { 
          if (currentNucleotide %in% nucleotides[[nucleotide]]) {
            runif(1, vScale, 1)
          } else {
            if (subSeq)
              runif(1, 0, subSeqVScale)
              # 0
            else
              runif(1, 0, backgroundVScale)
              # 0
          }
        }
        preWidth <- runif(1)
        curve <- cos(seq(pi, pi*3, pi/width)) * height
        # prenoise <- rnorm(100) %>% density
        # noise <- prenoise$y[seq(1, length(prenoise$y), length(prenoise$y) / baseNucleotideWidth)] / 100
        curve <- c(rep(0, round((baseNucleotideWidth - width) * preWidth)),
                   (curve - first(curve)),# / noise,
                   rep(0, round((baseNucleotideWidth - width) * (1 - preWidth))))
        # c(curve, rep(0, baseNucleotideWidth - width))
      }
      
      genCurve <- function(sequence, nucleotide, baseNucleotideWidth) {
        preCurve <- rep(0, (length(sequence) + 4) * baseNucleotideWidth)
        for (i in 1:length(sequence)) {
          position <- ((i + 1) * baseNucleotideWidth)
          curveElement <- genCurveElement(sequence[i], 
                                          nucleotide,
                                          input$subSeqFreq,
                                          baseNucleotideWidth,
                                          input$hScale,
                                          input$vScale,
                                          input$subSeqVScale,
                                          input$backgroundVScale)
          preCurve[position:(position + length(curveElement) - 1)] <- 
            preCurve[position:(position + length(curveElement) - 1)] + 
            curveElement
        }
        preCurve * 0.5
      }
      
      map(c("A" = "a", "T" = "t", "G" = "g", "C" = "c"),
          ~ genCurve(sequence(), ., input$peakWidth)
      ) %>% 
        as_tibble(.) %>% 
        mutate(seqTime = 1:((length(sequence()) + 4) * input$peakWidth)) %>% 
        gather("Nucleotide", "Signal", 1:4)
      # })
  })
  
  output$chromatogramPlot <- renderPlotly({
    req(chromatogram())
    ggplot(chromatogram(),
           aes(x = seqTime, y = Signal, color = Nucleotide)) +
      geom_line(size = 1) + 
      theme_bw()
  })

})

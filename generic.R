library(sangerseqR)

wInt <- function(obj, f, size) 
  writeBin(as.integer(obj),
           f, endian = "big", size = size)
wInt32 <- function(obj, f) wInt(obj, f, 4)
wInt16 <- function(obj, f) wInt(obj, f, 2)
wInt8 <-  function(obj, f) wInt(obj, f, 1)

wChar <- function(obj, f, nchars = nchar(obj)) 
  writeChar(obj, f, nchars = nchars, eos = NULL)


deltaIt <- function(samples) {
  p_delta  <-  0;
  for (i in 1:length(samples)) {
    p_sample <-  samples[i];
    samples[i] <-  samples[i] - p_delta;
    p_delta  <-  p_sample;
  }
  p_delta  <-  0;
  for (i in 1:length(samples)) {
    p_sample <-  samples[i];
    samples[i] <-  samples[i] - p_delta;
    p_delta  <-  p_sample;
  }
  samples
}

write.scf <- function(scf, filename) {
  fc <- file(filename, open = "wb")
  
  tryCatch({
    # write header
    wChar(".scf", fc, nchars = 4)
    wInt32(scf@header@samples, fc)
    wInt32(scf@header@samples_offset, fc)
    wInt32(scf@header@bases, fc)
    wInt32(scf@header@bases_left_clip, fc)
    wInt32(scf@header@bases_right_clip, fc)
    wInt32(scf@header@bases_offset, fc)
    wInt32(scf@header@comments_size, fc)
    wInt32(scf@header@comments_offset, fc)
    wChar(sprintf("%.2f", scf@header@version), fc, nchars = 4)
    wInt32(scf@header@sample_size, fc)
    wInt32(scf@header@code_set, fc)
    wInt32(scf@header@private_size, fc)
    wInt32(scf@header@private_offset, fc)
    wInt32(rep(0, 18), fc)
    
    # write samples vectors
    wInt16(deltaIt(scf@sample_points[, 1]), fc) # A
    wInt16(deltaIt(scf@sample_points[, 2]), fc) # C
    wInt16(deltaIt(scf@sample_points[, 3]), fc) # G
    wInt16(deltaIt(scf@sample_points[, 4]), fc) # T
    
    # write basecall positions
    wInt32(scf@basecall_positions, fc)
    
    # write basecall probabilities
    wInt8(scf@sequence_probs[, 1], fc) # A
    wInt8(scf@sequence_probs[, 2], fc) # C
    wInt8(scf@sequence_probs[, 3], fc) # G
    wInt8(scf@sequence_probs[, 4], fc) # T
    
    # write basecalls
    wChar(scf@basecalls, fc)
    
    # write reserved field
    wInt8(rep(0, nrow(scf@sequence_probs) * 3), fc)
    
    # write comments
    wChar(scf@comments, fc)
    
    # write private
    wInt8(scf@private, fc)
  },
  finally = close(fc)
  )
}

generatePeak <- function(pregion, 
                         pwidth, pwidthVariety,
                         pheight, pheightVariety#,
                         # noise
                         ) {
  peak <- dnorm(pregion, 
                mean = mean(pregion), 
                sd = pwidth * length(pregion) * 
                  abs(rnorm(1, mean = 1, sd = pwidthVariety)))
  # # add noise
  # peak <- peak + rnorm(pregion, sd = 
  #                        abs(rnorm(1, sd = noise)))
  # rescale peak
  peak <- punif(peak, min = min(peak), max = max(peak)) * pheight *
    abs(rnorm(1, mean = 1, sd = pheightVariety))
  peak
}

addNuc <- function(samplePoints,
                   nuc,
                   pregion,
                   pwidth, pwidthVariety,
                   pheight, pheightVariety) {
  # A, C, G, T only
  addRealNuc <- function(realNuc, pheight) {
    samplePoints[pregion, realNuc] <<- 
      samplePoints[pregion, realNuc] +
      generatePeak(pregion = pregion,
                   pwidth = pwidth, pwidthVariety = pwidthVariety,
                   pheight = pheight,
                   pheightVariety = pheightVariety#,
                   # noise = noise
      )
  }
  switch(nuc,
    "A" = { addRealNuc("A", pheight) },
    "C" = { addRealNuc("C", pheight) },
    "G" = { addRealNuc("G", pheight) },
    "T" = { addRealNuc("T", pheight) },
    "R" = { addRealNuc("A", pheight * 0.7); addRealNuc("G", pheight * 0.7)},
    "Y" = { addRealNuc("C", pheight * 0.7); addRealNuc("T", pheight * 0.7)},
    "S" = { addRealNuc("G", pheight * 0.7); addRealNuc("C", pheight * 0.7)},
    "W" = { addRealNuc("A", pheight * 0.7); addRealNuc("T", pheight * 0.7)},
    "K" = { addRealNuc("G", pheight * 0.7); addRealNuc("T", pheight * 0.7)},
    "M" = { addRealNuc("A", pheight * 0.7); addRealNuc("C", pheight * 0.7)},
    "B" = { addRealNuc("C", pheight * 0.45); addRealNuc("G", pheight * 0.45); addRealNuc("T", pheight * 0.45);},
    "D" = { addRealNuc("A", pheight * 0.45); addRealNuc("G", pheight * 0.45); addRealNuc("T", pheight * 0.45);},
    "H" = { addRealNuc("A", pheight * 0.45); addRealNuc("C", pheight * 0.45); addRealNuc("T", pheight * 0.45);},
    "V" = { addRealNuc("A", pheight * 0.45); addRealNuc("C", pheight * 0.45); addRealNuc("G", pheight * 0.45);},
    "N" = { addRealNuc("A", pheight * 0.3); addRealNuc("C", pheight * 0.3); addRealNuc("G", pheight * 0.3); addRealNuc("T", pheight * 0.3);},
    "*"
  )
  samplePoints
}

genPeakRegion <- function(peakPosition, nsamples) {
  if (peakPosition[1] < 1) peakPosition[1] <- 1
  if (peakPosition[2] > nsamples) peakPosition[2] <- nsamples
  peakPosition[1]:peakPosition[2]
}

generateSamplePoints <- 
  function(nucSeq, samplesPerNuc,
           ppositionVariety,
           pwidth, pwidthVariety,
           pheight, pheightVariety,
           # noise,
           subpSeq, subpProb, subpHeight,
           noisepProb, noisepHeight
           # k
           ) {
  nucs <- strsplit(toupper(nucSeq), "")[[1]]
  nsamples <- length(nucs) * samplesPerNuc
  
  samplePoints <- 
    matrix(0L,
           nrow = nsamples,
           ncol = 4L,
           dimnames = list(1:nsamples,
                           c("A", "C", "G", "T"))
    )
  
  positions <- data.frame(nuc = nucs,
                          pos = 0L,
                          A = 0L,
                          C = 0L,
                          G = 0L,
                          T = 0L,
                          stringsAsFactors = FALSE)
  if (!is.null(subpSeq)) {
    subpProb <- 1
    subpSeq <- strsplit(toupper(subpSeq), "")[[1]]
    subpSeq <- c(subpSeq, rep("_", length(nucs) - length(subpSeq)))
    subpSeqGenerator <- function(i) {
        subpSeq[i]
    }
  } else {
    subpSeqGenerator <- function(i) {
      sample(c("A", "C", "G", "T"), 1)
    }
  }
  
  for (i in seq_along(nucs)) {
    peakPositionDelta <- rnorm(1, sd = ppositionVariety)
    peakPosition <- c(
      as.integer((i - 1.5 + peakPositionDelta) * samplesPerNuc),
      as.integer((i + 1.5 + peakPositionDelta) * samplesPerNuc)
    )
    peakRegion <- genPeakRegion(peakPosition, nsamples)
    samplePoints <- 
      addNuc(samplePoints = samplePoints,
             nuc = nucs[i],
             pregion = peakRegion,
             pwidth = pwidth, pwidthVariety = pwidthVariety,
             pheight = pheight, pheightVariety = pheightVariety
      )
    # samplePoints[peakRegion, nucs[i]] <- 
    # samplePoints[peakRegion, nucs[i]] +
      # generatePeak(pregion = peakRegion,
      #              pwidth = pwidth, pwidthVariety = pwidthVariety,
      #              pheight = pheight, pheightVariety = pheightVariety#,
      #              # noise = noise
      #              )
    positions[i, "pos"] <- as.integer(mean(peakRegion))
    
    # add subpeak
    if (runif(1) <= subpProb) {
      peakPosition <- c(
        # more delta for random subpeaks
        as.integer((i - 1.5 + peakPositionDelta * 2) * samplesPerNuc),
        as.integer((i + 1.5 + peakPositionDelta * 2) * samplesPerNuc)
      )
      peakRegion <- genPeakRegion(peakPosition, nsamples)
      randomNuc <- subpSeqGenerator(i)
      
      samplePoints <- 
        addNuc(samplePoints = samplePoints,
               nuc = randomNuc,
               pregion = peakRegion,
               pwidth = pwidth, pwidthVariety = pwidthVariety,
               pheight = pheight * subpHeight, 
               pheightVariety = pheightVariety
        )
    }
    
    # add random position noise peak
    if (runif(1) <= noisepProb) {
      position <- runif(1) * nsamples
      peakPosition <- c(
        # more delta for random subpeaks
        as.integer(position),
        as.integer(position + abs(rnorm(1, mean = 1) * samplesPerNuc) + samplesPerNuc)
      )
      peakRegion <- genPeakRegion(peakPosition, nsamples)
      randomNuc <- sample(c("A", "C", "G", "T"), 1)
      
      samplePoints <- 
        addNuc(samplePoints = samplePoints,
               nuc = randomNuc,
               pregion = peakRegion,
               pwidth = pwidth * 4, pwidthVariety = pwidthVariety,
               pheight = pheight * noisepHeight, 
               pheightVariety = pheightVariety * 4
        )
      
      # samplePoints[peakRegion, randomNuc] <- 
      #   samplePoints[peakRegion, randomNuc] +
      #   generatePeak(pregion = peakRegion,
      #                pwidth = pwidth * 4, pwidthVariety = pwidthVariety,
      #                pheight = pheight * noisepHeight,
      #                pheightVariety = pheightVariety * 4)
    }
  }
  
  # # smooth curves and make integer
  # samplePoints <- 
  #   apply(samplePoints, 2, zoo::rollmean, k = k, na.pad = TRUE)
  samplePoints <- 
    apply(samplePoints, 2, as.integer)
  
  samplePoints[is.na(samplePoints) | samplePoints < 0] <- 0L
  # Old : Q = -10 * log10(nuc_area_under_peak / total_nucs_area_under_peak)
  # Q = -10 * log10(nuc_peak_height / sum_nucs_peaks_height)
  for (i in seq_along(nucs)) {
    # peakPosition <- c(positions[i, "pos"] - 0.5 * samplesPerNuc,
    #   positions[i, "pos"] + 0.5 * samplesPerNuc)
    # peakRegion <- genPeakRegion(peakPosition, nsamples)
    # positions[i, nucs[i]] <- 
    #   as.integer(-10 *
    #                log10(sum(samplePoints[peakRegion, nucs[i]]) / 
    #                        sum(samplePoints[peakRegion,])))
    for (nuc in c("A", "C", "G", "T")) {
      positions[i, nuc] <- 
        tryCatch(
          as.integer(-10 *
                       log10(1 - (samplePoints[positions[i, "pos"], nuc] / 
                                    sum(samplePoints[positions[i, "pos"],])))),
          warning = function(w) 40) # 0.0001 error
    }
  }
  
  list(samplePoints = samplePoints, 
       positions = positions)
}



create.scf <- function(nucSeq, 
                       samplesPerNuc = 20L,
                       ppositionVariety = 0.1,
                       pwidth = 0.1, pwidthVariety = 0.1,
                       pheight = 500, pheightVariety = 0.1,
                       # noise = 0.0005, 
                       subpSeq = NULL, subpProb = 0.2, subpHeight = 0.3,
                       noisepProb = 0.2, noisepHeight = 0.15,
                       # k = samplesPerNuc / 2,
                       comments = "Generated by falSequencer",
                       private = raw(1)) {
  scf <- new("scf")
  
  nucSeql <- nchar(nucSeq)
  
  peaks <- generateSamplePoints(
    nucSeq, samplesPerNuc = samplesPerNuc,
    ppositionVariety = ppositionVariety,
    pwidth = pwidth, pwidthVariety = pwidthVariety,
    pheight = pheight, pheightVariety = pheightVariety,
    # noise = noise, 
    subpSeq = subpSeq, subpProb = subpProb, subpHeight = subpHeight,
    noisepProb = noisepProb, noisepHeight = noisepHeight
    # k = k
    )
  
  scf@header@scf <- ".scf"
  scf@header@samples <- nucSeql * samplesPerNuc 
  scf@header@samples_offset <- 128L
  scf@header@bases <- nucSeql
  scf@header@bases_left_clip <- 0L
  scf@header@bases_right_clip <- 0L
  
  scf@header@bases_offset <- scf@header@samples_offset +
    scf@header@samples * 4L * # ACGT
    2L # 16 bits = 2 bytes
  
  scf@header@comments_size <- nchar(comments)
  scf@header@comments_offset <- scf@header@bases_offset +
    nucSeql * 4L +   # 32 bits = 4 bytes peak_index
    nucSeql * 4L +   # probabilities
    nucSeql +        # basecalls
    nucSeql * 3L     # spare
    
    
  scf@header@version <- 3L
  scf@header@sample_size <- 2L # Size of samples in bytes 1=8bits, 2=16bits
  scf@header@code_set <- 0L
  scf@header@private_size <- length(private)
  scf@header@private_offset <- scf@header@comments_offset +
    scf@header@comments_size
  
  scf@sample_points <- peaks$samplePoints
  scf@basecall_positions <- peaks$positions$pos
  scf@sequence_probs <- as.matrix(peaks$positions[, c("A", "C", "G", "T")])
  scf@basecalls <- nucSeq
  scf@comments <- comments
  scf@private <- private
  
  scf
}




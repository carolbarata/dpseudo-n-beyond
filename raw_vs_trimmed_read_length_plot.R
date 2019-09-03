#!/usr/bin/env Rscript

########################################################
#       Raw vs trimmed read length distribution        #
########################################################

args <- commandArgs(TRUE)
vargs <- strsplit(args, ",")
lengthdatafile <- as.character(vargs[[1]])
outfilename <- as.character(vargs[[2]])

plot_read_len <-
  function(lengthdatafile, outfilename) {
    lengthdata <- read.table(lengthdatafile, sep = ",")
    samplename <- head(strsplit(tail(strsplit(lengthdatafile, split = "/")[[1]], n = 1), split = "[.]")[[1]], n = 1)
    allreads <- data.frame(length = c(lengthdata$V1, lengthdata$V2), 
                           status = c(rep("Raw reads", each = length(lengthdata$V1)), rep("Trimmed reads", each = length(lengthdata$V2))))
    rm(lengthdata)
    gc()
    library(ggplot2)
    length_distn <- ggplot(allreads, aes(x = length, colour = status)) + 
      geom_freqpoly(binwidth = 1) + 
      ylab("Counts") + xlab("Read length") + 
      labs(title = paste0("Read length distribution for ", samplename))
    ggsave(filename = outfilename, plot = length_distn)
    rm(list = ls())
    gc()
  }

plot_read_len(lengthdatafile, outfilename)
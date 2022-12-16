#!/usr/bin/env Rscript

########################################################
#     Raw vs trimmed mapping quality distribution      #
########################################################

args <- commandArgs(TRUE)
vargs <- strsplit(args, ",")
rawmapqfile <- as.character(vargs[[1]])
trimmedmapqfile <- as.character(vargs[[2]])
outfilename <- as.character(vargs[[3]])


plot_map_quality <-
  function(rawmapqfile, trimmedmapqfile, outfilename) {
    rawmapq <- unlist(read.table(rawmapqfile))
    trimmedmapq <- unlist(read.table(trimmedmapqfile))
    samplename <- head(strsplit(tail(strsplit(rawmapqfile, split = "/")[[1]], n = 1), split = "-")[[1]], n = 1)
    allmapq <- data.frame(mapq = c(rawmapq, trimmedmapq), 
                           status = c(rep("Raw data", each = length(rawmapq)), rep("Trimmed data", each = length(trimmedmapq))))
    rm(rawmapq)
    rm(trimmedmapq)
    gc()
    library(ggplot2)
    mapq_distn <- ggplot(allmapq, aes(x = mapq, colour = status)) + 
      geom_freqpoly(binwidth = 1) + 
      ylab("Counts") + xlab("Mapping quality") + 
      labs(title = paste0("Mapping quality distribution for ", samplename))
    ggsave(filename = outfilename, plot = mapq_distn)
    rm(list = ls())
    gc()
  }

plot_map_quality(rawmapqfile, trimmedmapqfile, outfilename)

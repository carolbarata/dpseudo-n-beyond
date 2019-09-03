#!/usr/bin/env Rscript

########################################################
#             Coverage distribution plot               #
########################################################

args <- commandArgs(TRUE)
vargs <- strsplit(args, ",")
covdatafile <- as.character(vargs[[1]])
outfilename <- as.character(vargs[[2]])

plot_coverage <-
  function(covdatafile, outfilename) {
    library(ggplot2)
    covdata <- read.table(covdatafile)
    samplename <- head(strsplit(tail(strsplit(covdatafile, split = "/")[[1]], n = 1), split = "[.]")[[1]], n = 1)
    nosing_covdata <- covdata[which(!startsWith(as.character(covdata$V1), "Unknown_singleton") &
                                      covdata$V3 < 151), ]
    nosing_covdata$V1 <- droplevels(nosing_covdata$V1)
    rm(covdata)
    gc()
    cov_distn <- ggplot(data = nosing_covdata) + 
      geom_freqpoly(mapping = aes(V3), 
                    binwidth = 2) + 
      ylab("Counts") + xlab("Coverage") +
      xlim(0, 150) +
      labs(title = paste0("Genome-wide coverage distribution for ", samplename))
    ggsave(filename = outfilename, plot = cov_distn)
    rm(list = ls())
    gc()
  }

plot_coverage(covdatafile, outfilename)
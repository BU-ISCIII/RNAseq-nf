#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

require(ballgown)
my.data = ballgown(dataDir="ballgown/", samplePattern='*_ballgown', meas='all')

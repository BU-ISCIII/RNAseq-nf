#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

require(ballgown)
my.data = ballgown(dataDir="./", samplePattern='*_ballgown', meas='all')

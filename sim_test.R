#!/usr/bin/env Rscript

library(optparse) # to create a command line interface
library(sys) # to infer system environment variables
library(compiler) # to accelerate using a just in time compiler

options(echo=TRUE)

#we can create a default output directory path, here:
default.out = paste("/home/", Sys.getenv("USER"), sep="", collapse=NULL)

option_list <- list(
  make_option(c("-N","--OtherItems"), type="integer",
              help="Number of OtherItems",action = "store",
              metavar="number"),
  make_option(c("-K","--NPL"), type="integer",
              help="Number of NPLs",action = "store",
              metavar="number"),
  make_option(c("-F","--nFreetime"), type="integer",
              help="Number of Freetime Conditions",action = "store",
              metavar="number"),
  make_option(c("-R","--nRetrievals"), type="integer", action="store",
              help="Number of Retrievals",
              metavar="number")

)

args <- parse_args(OptionParser(option_list=option_list))

nFT <- args$nFreetime
nRet <- args$nRetrievals

N <- args$OtherItems
K <- args$NPL


print(nFT)
print(nRet)
print(N)
print(K)


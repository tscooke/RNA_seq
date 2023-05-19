# This code walks through the tutorial from Law et al's 2018 paper:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/

library(limma)
library(edgeR)

files <- c(
  "GSM1545535_10_6_5_11.txt",
  "GSM1545536_9_6_5_11.txt",
  "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt",
  "GSM1545540_JMS8-3.txt",
  "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt",
  "GSM1545544_JMS9-P7c.txt",
  "GSM1545545_JMS9-P8c.txt"
)

setwd(paste0(getwd(), "/00. Raw Data"))

read.delim(files[1], nrow = 5)
  # Prove to myself that data is in the right place and can read it into R

x <- readDGE(files, columns = c(1,3))
  # readDGE() can conveniently read in multiple sample files at the same time,
  # and will create a separate column for each sample
class(x)
dim(x)

# For convenience sake, we will remove the "GSM" text from sample IDs
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

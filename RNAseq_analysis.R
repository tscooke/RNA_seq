# This code walks through the tutorial from Law et al's 2018 paper:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/

library(limma)
library(edgeR)
library(Mus.musculus)

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
colnames(x) <- samplenames

group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
x$samples$group <- group

lane <- as.factor(rep(c("L004", "L006", "L008"), c(3,4,2)))
x$samples$lane <- lane
  # We added sample metadata to the 'samples' dataframe inside of DGEList x

x$samples

# Next, we're going to set up metadata for the genes included in x
geneid <- rownames(x)
genes <- select(
  Mus.musculus,
  keys = geneid,
  columns = c("SYMBOL", "TXCHROM"),
  keytype = "ENTREZID"
)
head(genes)

genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes
x

cpm <- cpm(x) # Calculate counts per million
lcpm <- cpm(x, log = TRUE) # Calculate log2-counts per million

# If we're curious about what the offset is for lcpm, can determine calculated input L
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

# We're not interested in genes that have low expression across all samples, so we'll filter those out
# Note that filtering for low expression happens BEFORE normalization

table(rowSums(x$counts == 0) == 9)
  # rowSums(x$counts == 0) is the number of columns where count == 0; the whole expression
  # tells us how many rows have 0 for all 9 columns

keep.exprs <- filterByExpr(x, group = group)
  # filterByExpr produces a logical vector, which we then use to extract a subset from x
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
dim(x)  

# To visually inspect the distribution of the data, and get a sense of where our cut-off should be,
# we can plot a density function of log-cpm
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples,"Paired")
par(mfrow = c(1,2)) # Setting graphical parameters for output to display as 1x2 array
plot(
  density(lcpm[,1]), # This is the 'lcpm' that we calculated BEFORE excluding low expression values
  col = col[1],
  lwd = 2,
  ylim = c(0, 0.26),
  las = 2,
  main = "",
  xlab = ""
)
title(main = "A. Raw Data", xlab = "Log-cpm")
abline(v = lcpm.cutoff, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", samplenames, text.col = col, bty = "n")

lcpm <- cpm(x, log = TRUE) # We recalculate 'lcpm' AFTER excluding low expression
plot(
  density(lcpm[,1]), # The second graph uses the NEW 'lcpm' to show filtered data only
  col = col[1],
  lwd = 2,
  ylim = c(0, 0.26),
  las = 2,
  main = "",
  xlab = ""
)
title(main = "B. Filtered Data", xlab = "Log-cpm")
abline(v = lcpm.cutoff, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", samplenames, text.col = col, bty = "n")


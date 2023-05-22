# This code walks through the tutorial from Law et al's 2018 paper:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/

library(limma)
library(edgeR)
library(Glimma)
library(gplots)
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

# Log-cpm isn't enough to account for full normalization of data, so we'll use TMM instead
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# Creating an exploratory multi-dimensional scaling plot via plotMDS() is the recommended 
# way to begin looking for differential gene expression between samples / conditions.
# MDS is basically a PCA under the hood
lcpm <- cpm(x, log = TRUE) # I think we're now incorporating normalization factors into lcpm calculation
par(mfrow = c(1,2))

col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

col.lane <- lane
levels(col.lane) <- brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)

plotMDS(lcpm, labels = group, col = col.group)
title(main = "A. Sample groups")
plotMDS(lcpm, labels = lane, col = col.lane, dim = c(3,4))
title(main = "B. Sequencing lanes")

# The Glimma package allows for interactive MDS and other kinds of plots
# Can visit the authors' interactive plot at:
  # http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MDS-Plot.html
glMDSPlot(
  lcpm,
  labels = paste(group, lane, sep = "_"),
  groups = x$samples[,c(2,5)],
  launch = FALSE
)

# Next we're actually going to start looking at the differential gene expression 
# To tell the algorithms what we actually want to compare, we make a design matrix
  # that delineates the different groups that we want to compare between, 
  # in this case, group and lane
design <- model.matrix(~0 + group + lane)
colnames(design) <- gsub("group", "", colnames(design))
design

# Since we have three groups, we will probably be interested in looking at pairwise
  # comparisons between the groups
contr.matrix <- makeContrasts(
  BasalvsLP = Basal - LP,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design)
)
contr.matrix

# Next we use voom to create precision-weight normalization factors to facilitate 
  # inter-sample comparisons and differential gene expression analysis
v <- voom(x, design, plot = TRUE)
v

# voom transforms the data so that we can safely use linear models on it
# The next step is to prepare those linear models
# We also use an empirical Bayes transformation "to obtain more precise estimates of gene-wise variability"
  # My understanding of empirical Bayes methods are that they're used to account for 
  # multiple levels of variation within data, but I'm not entirely sure what we're targeting here...
  # Maybe the fact that longer genes tend to have higher counts?
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

# Finally, we're ready to actually start looking at which genes are differentially expressed
summary(decideTests(efit))
  
#decideTests() uses an adjusted p-value cut-off of 5%, but we can use stricter criteria if we want
tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)

# We might want to ask which genes are differentially expressed in Basal vs LP AND Basal vs ML
de.common <- which(dt[,1] != 0 & dt[,2] != 0)
length(de.common)
head(tfit$genes$SYMBOL[de.common])
vennDiagram(dt[,1:2], circle.col = c("turquoise", "salmon"))
write.fit(tfit, dt, file = "results.txt")

# So far we've seen how to get a list of genes that are differentially expressed in either
# direction, but what if we want to see them in order
basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)
basal.vs.ml <- topTreat(tfit, coef = 2, n = Inf)
head(basal.vs.lp)

plotMD(
  tfit,
  column = 1,
  status = dt[,1],
  main = colnames(tfit)[1],
  xlim = c(-8, 13)
)

# We can use the Glimma package to make this graph ~interactive~

glMDPlot(
  tfit,
  coef = 1,
  status = dt,
  main = colnames(tfit)[1],
  side.main = "ENTREZID",
  counts = lcpm,
  groups = group,
  launch = TRUE
)

# Heatmaps can also be a helpful way to visualize RNA-seq data

basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000, "blue", "white", "red")

heatmap.2(
  lcpm[i,], 
  scale = "row",
  labRow = v$genes$SYMBOL[i],
  cexRow = 0.5,
  labCol = group,
  col = mycol,
  trace = "none",
  density.info = "none",
  margins = c(8,6),
  lhei = c(2,8.5),
  dendrogram = "column"
)

# Next, can perform gene set enrichment analysis from inside of R, using CAMERA instead of GSEA
  # CAMERA is better for smaller sample sets, uses linear modeling approach to establishing null hypothesis baseline

load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
idx <- ids2indices(Mm.c2, id = rownames(v))

cam.BasalvsLP <- camera(v, idx, design, contrast = contr.matrix[,1])
head(cam.BasalvsLP,5)

cam.BasalvsML <- camera(v,idx, design, contrast = contr.matrix[,2])
head(cam.BasalvsML, 5)

cam.LPvsML <- camera(v, idx, design, contrast = contr.matrix[,3])
head(cam.LPvsML)

barcodeplot(
  efit$t[,3],
  index = idx$LIM_MAMMARY_LUMINAL_MATURE_UP,
  index2 = idx$LIM_MAMMARY_LUMINAL_MATURE_DN,
  main = "LPvsML"
)

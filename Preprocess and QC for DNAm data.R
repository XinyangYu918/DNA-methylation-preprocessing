#### Preprocess and QC methylation data arrayed using the EPIC V2 platform ####
library(RColorBrewer)
library(minfi)
library(limma)
library(FlowSorted.Blood.EPIC)
library(dplyr)
library(SummarizedExperiment)
library(corpcor)

# Replace with your local working directory
setwd(".") 
datadir <- "."

#### Read raw methylation data ####
targets=read.metharray.sheet(datadir, "Methylation_Samplesheet_IMAGEN_FU3.csv", recursive=T)
RGset <- read.metharray.exp(targets = targets, verbose=T, force = T)

# Update annotation
RGset@annotation[["annotation"]] <- "20a1.hg38" 
RGset@annotation[["array"]] <- "IlluminaHumanMethylationEPICv2"
#save(RGset, file="./RGset_raw.rda")

# Extracte phenotype data from the sample sheet
pd <- pData(RGset) 

# Generate QC report
qcReport(RGset, 
         sampNames = pd$Sample_ID, 
         pdf = "./qcReport_raw.pdf")

#### Run preprocessRaw to get beta values ####
object <- preprocessRaw(RGset)

# Assign probes with their physical location in the genome
object <- mapToGenome(object)

# Convert raw methylation data
object <- ratioConvert(object, type="Illumina")

# Get beta values and phenotype data
beta <- getBeta(object) 
dat <- object
pd <- pData(dat) 

#### Remove data from the X and Y chromosomes ####
keepIndex <- which(!seqnames(dat) %in% c("chrX", "chrY"))
beta <- beta[keepIndex, ]
 
### Use MDS to check bad samples ####
mydist <- dist(t(beta[sample(1:nrow(beta),10000),]))
mds <- cmdscale(mydist)

# Plot batch effects
# For sentrix barcodes
pdf("./Batch_effect_mds.pdf", width=11, height=4)
plot(mds[,1], col = as.numeric(as.factor(pd[,"SentrixBarcode_A"]))+1,
     xlab="",
     ylab="First Principal component", 
     xaxt="n")
dev.off()

#### Use PCA to check batch effects ####
b <- beta - rowMeans(beta)

# Performe singular value decomposition (SVD), equivalent to principal component analysis with covariance matrix
ss <- fast.svd(b, tol = 0)

# Check variances explained by each component and plot PCA results
percvar <- ss$d^2/sum(ss$d^2)
pdf(file = "./PCA_distribution.pdf")
plot(percvar,
     xlab="Principal Components",
     ylab="Variance explained")
dev.off()
 
# Plot batch effects over methylation PCs (here I use the first 4 components)
for (i in c("SentrixBarcode_A", "SentrixPosition_A"))
{ 
  pdf(file = paste("./batch_", i, ".pdf", sep = ""))
  variabletocolor = pd[, i] 
  bob = levels(factor(variabletocolor))
  colors = match(variabletocolor, bob)
  pairs(ss$v[, 1:4], col = colors, labels = c("PC1", "PC2", "PC3", "PC4"))
  par(mfrow = c(2,2))
  boxplot(ss$v[, 1] ~ pd[,i], ylab = "PC1", xlab = i)
  boxplot(ss$v[, 2] ~ pd[,i], ylab = "PC2", xlab = i)
  boxplot(ss$v[, 3] ~ pd[,i], ylab = "PC3", xlab = i)
  boxplot(ss$v[, 4] ~ pd[,i], ylab = "PC4", xlab = i)
  dev.off()
}

# Save PCA results for downstream analysis
save(ss, file="./fast_svd.rda")

# Mark individuals out of normal range (i.e. median+3SD or median-3SD)
RM <- rep(FALSE, nrow(ss$v))
for (i in 1:4) 
{
  # Calculate median and standard deviation of each component
  Median <- median(ss$v[, i]) 
  SD <- sd(ss$v[, i])
         
  # Mark individuals outside 3SD range as TRUE
  RM <- RM|(abs(ss$v[, i] - Median) > 3*SD) 
}

# Get the total number of participants to be removed due to out of 3SD
sum(RM)


### Predict sex using methylation data ####
predictedSex <- getSex(dat, cutoff = -2) 
predictedSex <- as.numeric(as.factor(predictedSex[, 3]))

# Check number of each sex (predicted sex)
tabulate(predictedSex)

# Check number of each sex (self-reported sex)
tabulate(as.factor(pd$Gender))

#### Plot predicted and self-reported sex ####
pdf(file="./Sex_Plot.pdf")
Jitter1 <- jitter(as.numeric(as.factor(pd$Gender)))
Jitter2 <- jitter(as.numeric(as.factor(predictedSex)))
plot(Jitter1, 
     Jitter2,
     xlab = "Sex (self reported)", 
     ylab="Sex(predicted)", 
     xaxt = "n", 
     yaxt = "n")
axis(2, c(1, 2), c("Female", "Male"))
axis(1, c(1, 2), c("Female", "Male"))

# Recode sex
predictedSexRecode <- ifelse(predictedSex == 1, "Female", ifelse(predictedSex == 2, "Male", NA))
predictedSexRecode <- as.factor(predictedSexRecode)

# Match predicted sex and self-reported sex
idx <- pd[, "Gender"] == predictedSexRecode
#idx <- as.factor(pd[, "Gender"])== as.factor(predictedSex)
#Outliers <- pd[!idx, "File"] 
# To show IDs of mismatched samples
# text(Jitter1[!idx], Jitter2[!idx]-0.05, Outliers, cex=0.5) 
dev.off()
 
# Get the total number of participants to be removed due to sex discrepancy
sum(idx == "FALSE")
 
# Mark and remove all outliers
RM <- RM|(pd[,"Sample_ID"] %in% Outliers)
pd <- pd[-which(RM),]


#### Reload data without outliers ####
RGset <- read.metharray.exp(targets = pd, verbose = TRUE) 

# Save raw data with outliers removed
save(RGset, file="./RGset_QCed.rda")
 
# Run preprocessQuantile, samples of low intensity will be removed
RGset@annotation[["annotation"]] <- "20a1.hg38" 
RGset@annotation[["array"]] <- "IlluminaHumanMethylationEPICv2"
object <- preprocessQuantile(RGset, 
                          fixOutliers = TRUE, 
                          removeBadSamples = TRUE, 
                          badSampleCutoff = 10.5, 
                          quantileNormalize = TRUE, 
                          stratified = TRUE, 
                          mergeManifest = FALSE, 
                          sex = NULL, 
                          verbose = TRUE)
save(object, file="./Quan-norm.rda")
 
# Get beta values of normalized data
beta <- getBeta(object)

# Mark and removed probes on X and Y chromosomes
keepIndex <- which(!seqnames(object) %in% c("chrX", "chrY")) 
beta <- beta[keepIndex, ]

# Save beta values of the quantile normalised data
save(beta, file="./Beta_Quantile.rda")

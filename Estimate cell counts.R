# The following was downloaded from: https://github.com/immunomethylomics/FlowSorted.Blood.EPIC/issues/10
#### To estimate the counts of cell types for methylation data arrayed with EPIC V2 platform ####
library(minfi)
library(sesame)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(FlowSorted.Blood.EPIC)

# Load RGset and annotate the methylation array data
RGset = read.metharray.exp(workdir,recursive = TRUE)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGset)["annotation"] = "20a1.hg38"

# Preprocess using preprocessNoob
MSet <- preprocessNoob(RGset)

# Get bata values
Betas <- getBeta(MSet)
Betas <- sesame::betasCollapseToPfx(Betas) #you can also use ENmix::rm.cgsuffix(Betas) or other function to remove replicates 

# Load and filter optimised CpGs:
IDOLOptimizedCpGsBloodv2<- IDOLOptimizedCpGs[which(IDOLOptimizedCpGs%in%rownames(Betas))]
identical(rownames(IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,]), IDOLOptimizedCpGsBloodv2)

# Project cell types
propEPIC <- projectCellType_CP(
    Betas[IDOLOptimizedCpGsBloodv2, ],
    IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,],
    contrastWBC = NULL, nonnegative = TRUE,
    lessThanOne = FALSE
)

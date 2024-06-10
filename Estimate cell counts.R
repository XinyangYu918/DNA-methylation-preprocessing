# The following was downloaded from: https://github.com/immunomethylomics/FlowSorted.Blood.EPIC/issues/10
#### To estimate the counts of cell types for methylation data arrayed with EPIC V2 platform ####
library(minfi)
library(sesame)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(FlowSorted.Blood.EPIC)
library(reshape2)
library(ggplot2)

# Load RGset and annotate the methylation array data
#RGset = read.metharray.exp(workdir,recursive = TRUE)
load("./RGset_QCed.rda")
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
    contrastWBC = NULL, 
    nonnegative = TRUE,
    lessThanOne = FALSE
)

# Plot the proportion of each cell type
cellcount <- data.frame(cbind(rownames(propEPIC), propEPIC))
cellcount[c(2:7)] <- lapply(cellcount[c(2:7)], as.numeric)
names(cellcount) <- c("ID", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
data_melt <- melt(cellcount, id.vars = "ID")

p <- ggplot(data_melt, aes(x = variable, y = value)) +
geom_boxplot(fill = "gray") +
labs(x = "Cell Types", 
     y = "Proportion") +
  theme_minimal()
print(p)

# Perform PCA of cell count proportion
pca_result <- prcomp(cellcount[2:7])  
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create a data frame for PCA plotting
scree_data <- data.frame(
  Principal_Component = 1:length(variance_explained),
  Variance_Explained = variance_explained)

# Generate the scree plot
scree_plot <- ggplot(scree_data, aes(x = Principal_Component, y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_line(aes(x = Principal_Component, y = Variance_Explained), color = "black") +
  geom_point(aes(x = Principal_Component, y = Variance_Explained), color = "red") +
  labs(title = "", x = "Principal Component", y = "Proportion of Variance Explained") +
  theme_minimal()

print(scree_plot)

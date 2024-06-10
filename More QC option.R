##### More QC options, I applied them before performing EWAS analysis ####
library(minfi) 

# Load files
load('./RGset.rda') 
load('./Quan-norm.rda')

# Add SNP information to the data
objectWithSNPinfo <- addSnpInfo(object) 

# Drop probes that contain either an SNP at the CpG interrogation or at the single nucleotide extension
objectSNPQCed <- dropLociWithSnps(objectWithSNPinfo, snps = c("SBE", "CpG", "Probe"), maf = 0.05) 
rm(objectWithSNPinfo);  

# A detection p-value is returned for every genomic position in every sample
detP <- detectionP(RGset) 
Match1 <- match(colnames(objectSNPQCed), colnames(detP)) 
Match2 <- match(rownames(objectSNPQCed), rownames(detP)) 
detPSNPQCed <- detP[Match2[!is.na(Match2)], Match1[!is.na(Match1)]] 
rm(Match1, Match2, detP)

# Positions with non-significant p-values (typically >0.01) should not be trusted
failed <- detPSNPQCed > 0.01 
rm(detPSNPQCed)

# Get beta values
beta <- getBeta(objectSNPQCed) 

# Drop probes that failed quality control via the detection p-value in greater than 20% of samples
failedCG02 <- rowMeans(failed) > 0.2 

# Get the list of non-variable CpG sites i.e. those where beta values for all samples are ≤20% or ≥80%
ProbeInvar <- (rowSums(beta<=0.2)==ncol(beta))|(rowSums(beta>=0.8)==ncol(beta))  

# Mark probes with either all beta value <=0.2 or all beta value>=0.80 and generate a list of probes marked as invariant 
ListInvarProbe <- rownames(beta)[which(ProbeInvar)]  
rm(ProbeInvar)

# Remove sex chromosome probes
keepIndex=!seqnames(objectSNPQCed) %in% c("chrX", "chrY")  
rm(objectSNPQCed)

# Remove failed probes
keepIndex <- keepIndex&(!failedCG02) 
rm(failedCG02)

# Remove all probes with detected P-value > 0.01
beta[failed] <- NA  
rm(failed)

# Remove marked failed probes
betaQC <- beta[which(keepIndex),]  
rm(beta, keepIndex)

# Reformat beta value matrix for downstream analysis
Methy <- as.data.frame(t(betaQC)) 

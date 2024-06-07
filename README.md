# DNA-methylation-preprocessing
Preprocessing and QC steps for DNAm arrayed with EPIC v2 for IMAGEN, STRATIFY and ESTRA cohorts.

## Workflow of minfi package and different outputs
![image](https://github.com/XinyangYu918/DNA-methylation-preprocessing/assets/52769576/bfac942c-14d8-4a95-98a1-127ad3d1dd73)

## Outputs from preprocessing and QC
1. RGset.rda # This file contains raw intensity data (both methylated and unmethylated signals) from the microarray before any preprocessing or normalization steps.
2. beta_Quantile.rda # This file contains beta values that have undergone quantile normalization.
3. Quantile-norm.rda # This file contains data that have been processed with quantile normalization.
4. fast_svd.rda # This file contains results from a fast singular value decomposition (SVD) analysis, which is performed to identify and correct for batch effects or other confounding variables.
5. cellcount.rda # This file contains estimated cell-type proportions for each sample.
and all QC plots (see details from Preprocess and QC for DNAm data.R)

* For case-control cohorts, Fun-norm.rda (preprocess using preprocessFunnorm) is also provided. 

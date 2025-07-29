# DNA Methylation Preprocessing

This repository contains preprocessing and quality control (QC) workflows for DNA methylation (DNAm) data arrayed using the **Illumina EPIC v2 array**. The data were derived from the **IMAGEN**, **STRATIFY**, and **ESTRA** cohorts.

---

## Workflow Overview
Preprocessing and QC are performed using the [`minfi`](https://bioconductor.org/packages/release/bioc/html/minfi.html) package. 
The overall workflow includes raw data import, normalisation, batch correction, cell type deconvolution, and quality control. 
See full implementation in this script:  
[`Preprocess and QC for DNAm data.R`](https://github.com/XinyangYu918/DNA-methylation-preprocessing/blob/main/Preprocess%20and%20QC%20for%20DNAm%20data.R)

<img src="https://github.com/XinyangYu918/DNA-methylation-preprocessing/assets/52769576/bfac942c-14d8-4a95-98a1-127ad3d1dd73" alt="Workflow Image" width="400"/>

---

## Output Files
| File | Description |
|------|-------------|
| `RGset.rda` | Raw intensity data (methylated and unmethylated signals) prior to normalisation or preprocessing. |
| `beta_Quantile.rda` | Beta values obtained after quantile normalisation. |
| `Quantile-norm.rda` | Quantile-normalised intensity data. |
| `fast_svd.rda` | Results from fast singular value decomposition (SVD), used to detect and adjust for batch effects or confounding factors. |
| `cellcount.rda` | Estimated cell-type proportions per sample. |

**QC plots** are included. For documentation, see:  
- [`Methylation QC protocol_21062024.docx`](./Methylation%20QC%20protocol_21062024.docx)  
- [`Methylation QC_07062024.pptx`](./Methylation%20QC_07062024.pptx)

---

## Additional Probe-Level QC (Recommended)
Additional QC procedures are described in [`More QC option.R`](./More%20QC%20option.R), including:

1. **SNP-based probe filtering**  
   - Remove probes with SNPs at CpG, single base extension (SBE), or probe body (`MAF > 0.05`)

2. **Detection p-value filtering**  
   - Remove probes with detection p-values > 0.01  
   - Exclude probes failing in >20% of samples

3. **Non-variable CpG filtering**  
   - Remove probes with invariant methylation (beta ≤ 0.2 or ≥ 0.8 in all samples)

4. **Sex chromosome probe removal**  
   - Exclude probes on chrX and chrY

5. **Missing value handling**  
   - Set failed probes to `NA` in the beta matrix  
   - Retain only high-confidence probes in downstream analysis

---

## Probe Filtering Resources
To remove cross-reactive and polymorphic probes from your analysis, refer to:  
[https://github.com/markgene/maxprobes](https://github.com/markgene/maxprobes)


# HFpEF Survival Analysis

### Introduction
We performed an RNA-seq study to determine what genes predict survival in patients with HFpEF. Our results can be found in our [preprint](wikipedia.org). Here, we share all of the materials needed to reproduce our analysis.

### Data
All of the relevant files to reproduce the analysis can be downloaded at our [Synapse](https://www.synapse.org/Synapse:syn68157084/files/). The purpose of most of the files can be found in the comments in the `HFpEF_code.R` itself. However, if you want to skip all of the steps in the code and just work directly with the finalized objects, you can just use the `HFpEF_workspace.RData` file that is on Synapse.

The most relevant file on the Github is `HFpEF_code.R`, which contains the code to reproduce findings in the manuscript. Please note that there are some small functions that I drew from other packages and workspaces on my Github (mostly for converting between Ensembl gene and gene symbol), so the code doesn't work entirely off download. Again, if you want to avoid that headache, you can just directly use the Rdata file from the Synapse.

Otherwise, the other files on the Synapse include our clinical parameter data, as well as the list of what I used as the most important clinical parameters to be tested in our random forest models (`best_clin_param.txt`). At some point, I also did a quick dive into each gene identified as a predictor of HFpEF survival - some of my musings are included in the file called `genes.docx`.


### Dependencies
Most of the libraries used in our codebase can be found from CRAN or Bioconductor.

Please feel free to email or raise an issue if any of the code doesn't work as claimed!

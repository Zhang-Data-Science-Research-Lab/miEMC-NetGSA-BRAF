
# miEMC-NetGSA

The code files in this repository generate the results for the paper
"miRNA-Gene Activity Interaction Networks (miGAIN): Integrated joint
models of miRNA-gene targeting and disturbance in signaling pathways"
by Henry Linder and Yuping Zhang. This repository is hosted at the
following locations:

* [github.com/Zhang-Data-Science-Research-Lab/miEMC-NetGSA-BRAF](https://github.com/Zhang-Data-Science-Research-Lab/miEMC-NetGSA-BRAF)
* [github.com/mhlinder/miEMC-NetGSA-BRAF](https://github.com/mhlinder/miEMC-NetGSA-BRAF)

In order to run the code, the following resources are necessary:

* TCGA data omics data, downloaded and formatted using
  [TCGA-Assembler](http://www.compgenome.org/TCGA-Assembler/)
  software, version 2. The data must include:
  * Gene expression
  * miRNA expression
  * Gene methylation
  * Gene copy number The results in the paper were obtained using data
  published via the [NCI Genomic Data
  Commons](https://portal.gdc.cancer.gov/) Legacy Archive.
* The [mirDIP: microRNA Data Integration
  Portal](http://ophid.utoronto.ca/mirDIP/) database of miRNA-gene
  target interactions.

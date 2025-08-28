# snRNA-Seq Analysis Code for Pretangle Tau Model of Rat Hippocampus

This repository hosts all the code used for analyzing the single-nucleus RNA sequencing (snRNA-Seq) datasets in this study.

## Workflow

![Alt text](https://github.com/ziahasanz/hippocampus_paper/blob/main/snRNA-Seq%20data%20analysis%20worflow%20for%20hippocampus.png)
<p align="right"> <i>Image Created in <a href="https://BioRender.com">BioRender</a></i></p>

1. **Raw data processing**  
   FASTQ files for each sample were processed using **Cell Ranger v7.2** (10x Genomics).

2. **Cell Ranger outputs**  
   The `filtered_matrix.h5` files generated for each sample were used as input for downstream analysis.

3. **Quality control**  
   Each sample’s `filtered_matrix.h5` file was processed individually for quality control (QC) using **Trailmaker**, hosted by Parse Biosciences (https://www.parsebiosciences.com/data-analysis/).

4. **Aggregation**  
   The QC-passed files were aggregated, downloaded, and converted into a single AnnData file:  
   **`hippo_allcelltypesData.h5ad`**.

5. **Downstream analysis**  
   The aggregated file was used for clustering, differential gene expression analysis, enrichment analysis, and figure generation.

---

## Data Availability
- The **raw sequencing data** are deposited in the **Gene Expression Omnibus (GEO)** under accession number **GSE306235** (private reviewer access; will be made public upon acceptance).  
- The **processed datasets** (AnnData `.h5ad`, DEGs, and GO enrichment files) are available in the **Zenodo record** (private reviewer link: https://zenodo.org/

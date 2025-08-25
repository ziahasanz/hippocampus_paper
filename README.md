This repository hosts all codes used for analyzing snRNA-Seq datasets in this study.
The FASTq files of each sample are processed using Cell Ranger v7.2 (10x Genomics).
The Cell Ranger output files (filtered_matrix.h5) for each sample are utilized for the next step in the analysis.
Each sample's filtered_matrix.h5 files are individually processed in the Trailmaker, hosted by Parse Biosciences (https://www.parsebiosciences.com/data-analysis/).
The processed files for each sample are aggregated, downloaded and converted into a single anndata file (hippo_allcelltypesData.h5ad).
The aggregated file is then used for downstream analysis and figure generation.

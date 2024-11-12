# Ribosome Profiling Analysis Pipeline

This repository contains the analysis pipeline for ribosome profiling and RNA-seq data examining translational regulation during early C. elegans development in wild-type and OMA-1(zu405) mutant embryos.

## Overview

The analysis examines:
- Quality control metrics for ribosome profiling data
- Differential expression analysis between developmental stages
- Translational efficiency calculations
- Clustering analysis of gene expression patterns
- RNA binding protein motif analysis in wild-type and OMA-1(zu405) mutants

## Data Access

### Pre-processed Data
For ease of use, this repository includes pre-processed data files generated using the RiboFlow pipeline:
- `ITP.ribo`: Processed ribosome profiling data covering wild-type and OMA-1(zu405) mutant embryos
- `RNASEQ.ribo`: Stage-matched (but not sample-matched) RNA-seq data

### Experimental Design
The dataset includes:
- Wild-type samples at 1-cell, 2-cell, 4-cell, and 8-cell stages
- OMA-1(zu405) mutant samples at 1-cell, 2-cell, and 4-cell stages
- Multiple replicates per condition

### Raw Data
Raw sequencing data is available from GEO under accession number [GSE281412]. This includes:
- Ribosome profiling data for wild-type and OMA-1(zu405) mutant embryos
- RNA-seq data from equivalent developmental stages
- Multiple replicates for each developmental stage

## Dependencies

### R Packages
```r
library(ribor)          # For .Ribo file handling
library(cowplot)        # For plot layouts
library(tidyverse)      # For data manipulation
library(edgeR)          # For differential expression
library(biomaRt)        # For gene annotation
library(clusterProfiler) # For GO analysis
library(org.Ce.eg.db)   # C. elegans annotation database
library(pheatmap)       # For heatmaps
library(ggpubr)         # For publication-ready plots
library(S4Vectors)      # For data structures
library(sva)            # For batch correction
library(scales)         # For plot scaling
library(GO.db)          # For GO analysis
```

## File Organization
```
Celegans_translation_git/
├── input_data/
│   ├── ITP.ribo                                     # Processed ribosome profiling data
│   ├── RNASEQ.ribo                                 # Processed RNA-seq data
│   ├── UTR3_list.fa                                # 3'UTR sequences
│   ├── genetics.114.168823-7.xls                   # OMA-1 pulldown dataset
│   ├── embj2010334-sup-0001.xls                    # GLD-1 pulldown dataset
│   ├── enrichment_df_modified.csv                  # Modified enrichment data
│   └── multimodal_enrichment_df_modified.csv       # Modified multimodal enrichment data
├── script/
│   ├── Analysis_and_figures.qmd                    # Main analysis script
│   └── riboflow_scripts/
│       ├── project_umi_RNAseq_riboflow.yaml        # RNA-seq processing configuration
│       └── project_umi_riboITP_riboflow.yaml       # Ribo-seq processing configuration
└── output/
    ├── Figures/
    │   ├── Figure_1/
    │   ├── Figure_2/
    │   ├── Figure_3/
    │   ├── Figure_4/
    │   ├── Figure_5/
    │   ├── Figure_6/
    │   ├── Figure_7/
    │   ├── Figure_S1/
    │   ├── Figure_S2/
    │   ├── Figure_S3/
    │   ├── Figure_S4/
    │   └── Figure_S5/
    └── Tables/
```

## Analysis Steps

### 1. Quality Control
- CDS percentage calculation
- Start/stop site coverage analysis
- Read length distribution analysis
- Replicate correlation assessment
- Batch effect analysis and correction

### 2. Data Processing
```r
# Filter by CPM thresholds:
# RNA-seq: 10 CPM in ≥18 samples
# Ribo-seq: 3 CPM in ≥10 samples

# Apply batch correction using ComBat_seq
# Covariates: cell stage and genotype (WT/zu405)
```

### 3. Differential Expression Analysis
- Stage-wise comparisons using edgeR
- Separate analysis for:
  - RNA abundance
  - Ribosome occupancy
  - Translational efficiency

### 4. Clustering Analysis
- K-means clustering of gene expression patterns
- GO enrichment analysis of clusters
- Visualization of cluster profiles

### 5. RNA Binding Protein Analysis
- Motif identification in 3'UTRs
- Integration with pulldown data
- Linear modeling of binding effects

## Running the Analysis

### Quick Start
1. Clone this repository:
```bash
git clone https://github.com/yss322/Low-input-ribosome-profiling-of-Celegans-embryogenesis.git
cd Celegans_translation_git
```

2. Install required packages:
```r
# Install CRAN packages
install.packages(c("tidyverse", "cowplot", "ggpubr", "pheatmap", "scales"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ribor", "edgeR", "biomaRt", "clusterProfiler", 
                      "org.Ce.eg.db", "S4Vectors", "sva", "GO.db"))
```

3. Open `script/Analysis_and_figures.qmd` in RStudio

4. Click "Run All" in RStudio (or use Ctrl+Alt+R / Cmd+Alt+R)
   - The script automatically sets the working directory
   - All required input files are pre-organized in `input_data/`
   - Output will be generated in the `output/` directory

### Data Sources
- `genetics.114.168823-7.xls`: OMA-1 pulldown data (Spike et al. 2014)
- `embj2010334-sup-0001.xls`: GLD-1 pulldown data (Wright et al. 2011)
- Raw sequencing data: Available at GEO [INSERT NUMBER]

### RiboFlow Configuration
For users interested in processing raw data, RiboFlow configuration files are provided in `script/riboflow_scripts/`:
- `project_umi_riboITP_riboflow.yaml`: Parameters for ribosome profiling data
- `project_umi_RNAseq_riboflow.yaml`: Parameters for RNA-seq data

## Notes

- Histone genes are removed due to lack of polyA tails in C. elegans
- BatchQC was used to determine optimal batch correction parameters
- RNA binding protein motifs are based on published consensus sequences
- The zu405 allele is a temperature-sensitive allele of OMA-1

## Output Files

The analysis generates:
- Publication-ready figures in organized subdirectories
- Data tables containing:
  - Raw and processed count data
  - Differential expression results
  - GO enrichment analyses
  - Motif analysis results

## References

1. RiboFlow: Ozadam et al. (2020) Bioinformatics
2. EdgeR: Robinson et al. (2010)
3. ComBat-seq: Zhang et al. (2020)
4. OMA-1 pulldown data: Spike et al. (2014)
5. GLD-1 pulldown data: Wright et al. (2011)

## Contact

For questions or issues, please open an issue in this repository.

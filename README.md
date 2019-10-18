# Funm6AViewer
Visualization of single base differential m6A methylation sites and functional DmM genes.

# Installation

Funm6AViewer depends on GenomicFeatures, Guitar, trackViewer, DESeq2, STRINGdb, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db R packages and please make sure they are installed before installing Funm6AViewer.

1. Install the required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

BiocManager::install("Guitar", version = 3.8)

BiocManager::install("trackViewer")

BiocManager::install("DESeq2")

BiocManager::install("STRINGdb")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

BiocManager::install("org.Hs.eg.db")

2. Install Funm6AViewer

library("devtools")

install_github("NWPU-903PR/Funm6AViewer")

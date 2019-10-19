# Funm6AViewer
Visualization of single base differential m6A methylation sites and functional DmM genes.

# Installation

Funm6AViewer depends on GenomicFeatures, Guitar, trackViewer, DESeq2, STRINGdb, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db R packages and please make sure they are installed before installing Funm6AViewer.

1. Install the required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

BiocManager::install("Guitar")

BiocManager::install("trackViewer")

BiocManager::install("DESeq2")

BiocManager::install("STRINGdb")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

BiocManager::install("org.Hs.eg.db")

2. Install Funm6AViewer

library("devtools")

install_github("NWPU-903PR/Funm6AViewer")


---
title: "Funm6AViewer_Usage"
author: "Songyao"
date: "2019/10/19"
output:
  word_document: default
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

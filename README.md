# Funm6AViewer
Visualization of single base differential m6A methylation sites and functional DmM genes.

## 1. Installation

Funm6AViewer depends on GenomicFeatures, Guitar, trackViewer, DESeq2, STRINGdb, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db R packages and please make sure they are installed before installing Funm6AViewer.

### 1.1 Install the required packages

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
BiocManager::install("Guitar")
BiocManager::install("trackViewer")
BiocManager::install("DESeq2")
BiocManager::install("STRINGdb")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("org.Hs.eg.db")
```

### 1.2 Install Funm6AViewer

```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("NWPU-903PR/Funm6AViewer")
```

## 2. Data required

Funm6AViewer adopted FunDMDeep-m6A to idenify functional DmM genes which required 4 PPI networks. Then to run Funm6AViewer, users need to firstly download this data from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset

## 3. One step usage of Funm6AViewer

`funm6aviewer` takes single base DmM sites information, gene DE information as input and output:

1. Functional DmM genes (FDmMGenes);

2. DmM sites distribution on RNA;

3. Counts of DmM sites on different RNA regions;

4. DmM sites and reads coverage on interested gene;

5. Function enrichment of FDmMGenes;

6. Context specific function of interested genes;

7. DmMGene's MSB score along with DE score;

8. Network of FDmMGene's MSB neighbors.

Following is an example to achieve these using `funm6aviewer`:

Get input data:
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)
```

`dminfo` contains the position annotation
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)
```

   chr chromStart  chromEnd  name     score strand    log2fd
1 chr1  155160832 155160833  4582 0.9137714      - 2.7950754
2 chr1  171505224 171505225 23215 0.9431386      + 0.7589479
3 chr1  241767682 241767683 23596 0.8125095      - 1.0680801
4 chr1  243418399 243418400  9859 0.9652586      - 1.7805035
5 chr1    8073372   8073373 54206 0.8477583      - 2.5678589
6 chr1    8073689   8073690 54206 0.8089832      - 1.1375522

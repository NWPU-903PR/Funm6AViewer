# Funm6AViewer
Identification and visualization of  functional differential m6A methylation genes (FDmMGenes) and single base DmM sites.

## 1. Installation

Funm6AViewer depends on GenomicFeatures, Guitar, trackViewer, DESeq2, STRINGdb, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db R packages and please make sure they are installed before installing Funm6AViewer. An R version >= 3.6 is required.

Install the required packages
```{r, eval=FALSE}
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

Install Funm6AViewer
```{r, eval=FALSE}
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

```{r}
library(Funm6AViewer)
```

Get input data:
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)
```

`dminfo` contains the position annotation and log2 foldchange of DmM sites. It can be extracted from the result of DMDeep-m6A package using `summarydmdeepm6A`. Alternatively, users can use any other method to make it as the following formate: 
```{r}
head(dminfo)
```
The 'name' column can be entrez gene ID or gene symbol.

`deinfo` contains the differential expresion p-value and fdr for genes. It can be made using `makegrreadsfrombam` and `getdeinfo`, or users can use any other method to make it as the following formate:
```{r}
head(deinfo)
```
The 'name' column can be entrez gene ID or gene symbol.

`bamreadsgr` can be generated using `makegrreadsfrombam` from the MeRIP-Seq data in bam formate.
```{r, eval=FALSE}
bamreadsgr <- system.file("extdata", "bamgrlist_toy.RData", package="Funm6AViewer")
load(bamreadsgr)

siggene <- c("CCNT1", "MYC", "BCL2")
permutime <- 1000
```

The `datapath` is the filepath where the required PPI data saved and the `enrich_input_directory` is the filepath passed to string_db, the GO and KEGG function annotation data will be downloaded to this path. All these required data can be downloaded from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
```{r, eval=FALSE}
datapath <- "F:/Funm6A_package/data"
enrich_input_directory <- "F:/Funm6A_package/data"

savepath <- getwd()
re <- funm6aviewer(dminfo, deinfo, grlist, intrested_gene =  siggene, permutime = permutime, version = "11",
                   datapath = datapath, enrich_input_directory = enrich_input_directory, savepath = savepath)
```
The `version` parameter is passed to STRINGdb, it depends on which version is supported by STRINGdb. The results will be saved to `savepath`.

## 4. DmM sites plot for interested gene 

Users can use `dmsiteplot` to visaulize the DmM sites on their interested genes and their isoforms.

Get input:
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)

siggene <- c("MYC")
```

Make plot:
```{r}
re <- dmsiteplot(dminfo = dminfo, intrested_gene = siggene)
```

## 5. Reads coverage plot for interested gene

Users can use `coverageplot` to visaulize the reads coverage of DmM sites on their interested genes.

Get input:
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)

bamreadsgr <- system.file("extdata", "bamgrlist_toy.RData", package="Funm6AViewer")
load(bamreadsgr)

siggene <- c("MYC")
```

Make plot:
```{r}
re <- coverageplot(dminfo = dminfo, grlist = grlist, intrested_gene = siggene)
```

## 6. FunDMDeep-m6A

If users only hase a list of DmM genes and their DE information, then they can use `fdmdeepm6A` to identify functional DmM genes (FDmMGenes).

`DMgene` is a group of DmM genes, can be gene symbol or entrez gene ID.
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)

DMgene <- unique(dminfo$name)
descore <- getdescore(deinfo)
```

The `datapath` is the filepath where the required PPI data saved and it can be downloaded from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
```{r}
datapath <- "F:/Funm6A_package/data"
permutime <- 1000
```

Identify FDmMGenes:
```{r}
re <- fdmdeepm6A(DMgene = DMgene, descore = descore, datapath = datapath, permutime = permutime)
```

Plot interested genes' MSB score:
```{r}
siggene <- c("CCNT1", "MYC", "BCL2")
siggenescoreplot(fdmgene = re, siggene = siggene)
```

## 7. Context specific function annotation of interested FDmMGenes

Users can visualize the context specific function of interested FDmMGenes identified by FunDMDeep-m6A using `siggenepathplot`. 

`fdmgene` is a group of identified FDmMGnes. `siggene` is interested FDmMGne.
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
fdmgene <- unique(dminfo$name)

siggene <- c("MYC")
```

The `input_directory` is the filepath passed to string_db, the GO and KEGG function annotation data will be downloaded to this path. Users can also donwloaded the annotation data previously from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset and set the `input_directory` as where you save the data.
```{r}
input_directory <- "F:/Funm6A_package/data"

re <- siggenepathplot(fdmgene = fdmgene, intrested_gene = siggene, 
                      version = "11", input_directory = input_directory)
```

## 8. MSB net plot for interested FDmMGenes

Users can visualize the MSB neighbours of interested FDmMGenes identified by FunDMDeep-m6A using `msbnetplot`.

`dmgene` is a group of DmM gnes. `siggene` is interested FDmMGnes.
```{r}
dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)

dmgene <- unique(dminfo$name)
descore <- getdescore(deinfo)

siggene <- c("CCNT1", "MYC", "BCL2")
```

The `datapath` is the filepath where the required PPI data saved and it can be downloaded from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
```{r}
datapath <- "F:/Funm6A_package/data"
```

Plot for one FDmMGne:
```{r}
re <- msbnetplot(genesymbol = siggene[1], dmgene = dmgene, descore = descore, datapath = datapath)
```

plot for several FDmMGnes:
```{r}
re <- msbnetplot(genesymbol = siggene, dmgene = dmgene, descore = descore, datapath = datapath,
                 savename = "InterestedGene")
```

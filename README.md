# Funm6AViewer
Identification and visualization of  functional differential m6A methylation genes (FDmMGenes) and single base DmM sites.

## 1. Installation

Funm6AViewer depends on GenomicFeatures, Guitar, trackViewer, DESeq2, STRINGdb, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db R packages and please make sure they are installed before installing Funm6AViewer. An R version >= 3.6 is required.

Install the required packages
```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "Guitar", "trackViewer", "DESeq2", "STRINGdb",    
                       "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"), version = "3.10")
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

`dminfo` contains the position annotation and log2 foldchange of DmM sites. It can be extracted from the result of DMDeep-m6A package using `summarydmdeepm6A` (see [9.2](https://github.com/NWPU-903PR/Funm6AViewer#92-making-input-for-funm6aviewer) for more details). Alternatively, users can use any other method to make it as the following formate: 
```{r}
head(dminfo)
##    chr chromStart  chromEnd  name     score strand    log2fd
## 1 chr1  155160832 155160833  4582 0.9137714      - 2.7950754
## 2 chr1  171505224 171505225 23215 0.9431386      + 0.7589479
## 3 chr1  241767682 241767683 23596 0.8125095      - 1.0680801
## 4 chr1  243418399 243418400  9859 0.9652586      - 1.7805035
## 5 chr1    8073372   8073373 54206 0.8477583      - 2.5678589
## 6 chr1    8073689   8073690 54206 0.8089832      - 1.1375522
```
The 'name' column can be entrez gene ID or gene symbol.

`deinfo` contains the differential expresion p-value and fdr for genes. It can be made using `makegrreadsfrombam` and `getdeinfo` (see [9.2](https://github.com/NWPU-903PR/Funm6AViewer#92-making-input-for-funm6aviewer) for more details), or users can use any other method to make it as the following formate:
```{r}
head(deinfo)
##        name         pval         padj
## 1         1 7.578860e-01 8.990406e-01
## 2       100 6.958592e-01 8.695820e-01
## 3      1000 4.155368e-06 6.420489e-05
## 4     10000 2.043250e-02 9.424864e-02
## 5 100009676 4.524888e-01 7.148682e-01
## 6     10001 6.708161e-01 8.566831e-01
```
The 'name' column can be entrez gene ID or gene symbol.

`bamreadsgr` can be generated using `makegrreadsfrombam` from the MeRIP-Seq data in bam formate (see [9.2](https://github.com/NWPU-903PR/Funm6AViewer#92-making-input-for-funm6aviewer) for more details).
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
re <- funm6aviewer(dminfo, deinfo, grlist, intrested_gene =  siggene, permutime = permutime, version = "10",
                   datapath = datapath, enrich_input_directory = enrich_input_directory, savepath = savepath)
```
The `version` parameter is passed to STRINGdb, it depends on which version is supported by STRINGdb. The results will be saved to `savepath`.

If you are using other genomes, you need to install txdb and orgdb annotation for the corresponding genome. Taking mouse mm9 genome as an example, you should firstly install the genome annotation:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))    
    install.packages("BiocManager")
    
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm9.knownGene",    
                       "org.Mm.eg.db"))
```                    
And assign the annotation to `txdb`, `orgdb` and `orgsymbol` when running funm6aviewer:
```{r}
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
re <- funm6aviewer(dminfo = dminfo,    
                   deinfo = deinfo,    
                   grlist = grlist,    
                   intrested_gene =  siggene,    
                   txdb = TxDb.Mmusculus.UCSC.mm9.knownGene,    
                   orgdb = org.Mm.eg.db,    
                   orgsymbol = org.Mm.egSYMBOL,    
                   permutime = permutime,    
                   datapath = datapath,    
                   enrich_input_directory = enrich_input_directory,    
                   savepath = savepath)
```

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

## 9. Functional m6A analysis pipline from bam files

We introduced this pipeline using a pretend example MeRIP-Seq data which contains 2 replicates in bam format for each condition (untreated and treated), named as:

                             "untreated_input_rep1.bam","untreated_ip_rep1.bam",     
                             "untreated_input_rep2.bam", "untreated_ip_rep2.bam",     
                               "treated_input_rep1.bam", "treated_ip_replicate1.bam",     
                               "treated_input_rep2.bam", "treated_ip_rep2.bam".    
                               
These bam format files can be obtained by aligning the raw sequenced MeRIP-Seq data to genome, i.e., human hg19 using splice junction mapping tools for RNA-Seq reads, like TopHat, HISAT or STAR.

### 9.1 Calling differential m6A methylation sites using DMDeepm6A

Users can firstly call single base differential m6A methylation (DmM) sites using DMDeepm6A R package (https://github.com/NWPU-903PR/DMDeepm6A1.0) as following:
```{r}
library(DMDeepm6A)
ip_bams <- c("treated_ip_rep1.bam"," treated_ip_rep2.bam",     
             "untreated_ip_rep1.bam ", "untreated_ip_rep2.bam")     
input_bams <- c("treated_input_rep1.bam","treated_input_rep2.bam",     
                "untreated_input_rep1.bam", "untreated_input_rep2.bam")    
sample_condition <- c("treated", "treated", "untreated", "untreated")    
```
Assuming the aligned genome is hg19, users can call DmM sites as following:
```{r}
output_filepath <- getwd()    
re <- dmdeepm6A(ip_bams = ip_bams,    
                input_bams = input_bams,    
                sample_conditions = sample_condition,    
                output_filepath = output_filepath,    
                experiment_name = "DMDeepm6A_out")
```
See `?dmdeepm6A` for more details if you are using other genomes. The `experiment_name` is the name of the file folder where saved the result of `dmdeepm6A`.

### 9.2 Making input for Funm6AViewer

Funm6AViewer takes single base DmM sites information, gene DE information as input and a list of `GRanges` converted from MeRIP-Seq bam files using `makegrreadsfrombam` if users would like to see the reads coverage of interested DmM sites on gene. Single base DmM sites information named `dminfo` contains the position annotation and log2 foldchange of DmM sites. It can be extracted from the result of `DMDeepm6A` using `summarydmdeepm6A` as following:
```{r}
dminfo <- summarydmdeepm6A(dmpath = " DMDeepm6A_out", sigthresh = 0.05)
```
`dminfo` contains the following information:
```{r}
head(dminfo)
##    chr chromStart  chromEnd  name     score strand    log2fd
## 1 chr1  155160832 155160833  4582 0.9137714      - 2.7950754
## 2 chr1  171505224 171505225 23215 0.9431386      + 0.7589479
## 3 chr1  241767682 241767683 23596 0.8125095      - 1.0680801
## 4 chr1  243418399 243418400  9859 0.9652586      - 1.7805035
## 5 chr1    8073372   8073373 54206 0.8477583      - 2.5678589
## 6 chr1    8073689   8073690 54206 0.8089832      - 1.1375522
```
The ‘name’ column can be entrez gene ID or gene symbol.

Gene DE information named `deinfo` contains the differential expresion p-value and fdr for genes. It can be made using `makegrreadsfrombam` and `getdeinfo` from MeRIP-seq input samples as following:
```{r}
ip_bams <- c("treated_ip_rep1.bam"," treated_ip_rep2.bam",    
             "untreated_ip_rep1.bam ", "untreated_ip_rep2.bam")      
input_bams <- c("treated_input_rep1.bam","treated_input_rep2.bam",    
                "untreated_input_rep1.bam", "untreated_input_rep2.bam")    
sample_condition <- c("treated", "treated", "untreated", "untreated")    
savepath <- getwd()    
grlist <- makegrreadsfrombam(IP_bams = ip_bams,    
                             Input_bams = input_bams,    
                             condition = sample_condition,    
                             minimal_alignment_MAPQ = 30,    
                             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,    
                             savepath = savepath)    
deinfo <- getdeinfo(grlist = grlist,    
                    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,    
                    savepath = savepath)    
```

The converted `grlist` will be saved to `savepath` named as " bamgrlist.RData ". Users can directly load it to use next time. If you are using other genomes please install the corresponding txdb annotation similar to `TxDb.Hsapiens.UCSC.hg19.knownGene` of the genome and assign it to `txdb`.
`deinfo` contains the following information:
```{r}
head(deinfo)
##        name         pval         padj
## 1         1 7.578860e-01 8.990406e-01
## 2       100 6.958592e-01 8.695820e-01
## 3      1000 4.155368e-06 6.420489e-05
## 4     10000 2.043250e-02 9.424864e-02
## 5 100009676 4.524888e-01 7.148682e-01
## 6     10001 6.708161e-01 8.566831e-01
```
The ‘name’ column can be entrez gene ID or gene symbol.

### 9.3 One step usage of Funm6AViewer

Following is an example to achieve Functional m6A analysis using `funm6aviewer` with previously generated `dminfo`, `deinfo` and `grlist`:
```{r}
siggene <- c("CCNT1", "MYC", "BCL2")
permutime <- 100*length(unique(dminfo$name))
```
The `datapath` is the file path where the required PPI data saved and the `enrich_input_directory` is the file path passed to `string_db`, the GO and KEGG function annotation data will be downloaded and saved to this path. All these required data can be downloaded from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
```{r}
datapath <- "~./Funm6AViewer_data"
enrich_input_directory <- "~./Funm6AViewer_data"
savepath <- getwd()
re <- funm6aviewer(dminfo = dminfo,    
                   deinfo = deinfo,    
                   grlist = grlist,    
                   intrested_gene =  siggene,    
                   permutime = permutime,    
                   datapath = datapath,    
                   enrich_input_directory = enrich_input_directory,    
                   savepath = savepath)
```
The results will be saved to `savepath`. 

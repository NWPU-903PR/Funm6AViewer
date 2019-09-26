\name{funm6aviewer}
\alias{funm6aviewer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of single base DmM sites and FDmMGenes.
}
\description{
This function is used to visualize single base differential m6A methylation sites and functional DmM genes.
}
\usage{
funm6aviewer(dminfo,
             deinfo,
             bamgrlist = NA,
             siggene = NA,
             entrezid = NA,
             fungenethr = 0.3,
             permutime = 10^4,
             descoretype = "pval",
             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
             orgdb = org.Hs.eg.db,
             orgsymbol = org.Hs.egSYMBOL,
             datapath = NA,
             rescore_thr = 0.8,
             descore_thr = 0.8,
             enrich_input_directory = "",
             version = "10",
             species = 9606,
             bp_fdr_thr = 0.05,
             kegg_fdr_thr = 0.05,
             savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dminfo}{
A dataframe of DmM sites information. It can be generated using summarydmdeepm6A from DMDeepm6A package result.
}
  \item{deinfo}{
A dataframe of gene differential expression information
}
  \item{bamgrlist}{
A list containing the reads count of IP/Input bam files under different condition and the id of each files. It can be generated using  makegrreadsfrombam function. This is necessary if users want to watch the reads coverage for interested genes.
}
  \item{siggene}{
A vector of the gene symbol of interested genes to visualize
}
  \item{entrezid}{
A vector of the entrez id of interested genes to visualize
}
  \item{fungenethr}{
The threshold of the FDR used to determine whether a DmMGene is a FDmMGene
}
  \item{permutime}{
The permutation times to calculate an empirical p-value for DmMGenes to be FDmMGenes
}
  \item{descoretype}{
The value used to calculate the DE score. It can be "pval" or "padj" denotes use the p-value or FDR of genes to calculate the DE score.
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if use other genomes instead of hg19.
}
  \item{orgdb}{
The annotation of genome. You need to input it similar to "org.Hs.eg.db" if DefaultGenome is use other genomes instead of human genome.
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
  \item{datapath}{
The file path where the network information required of FunDMDeepm6A, usually do not need to input, only if you prefer to use your own PPI networks
}
  \item{rescore_thr}{
The threshold of relative MSB score percentage used to determine FDmMGenes with relatively higher MSB score
}
  \item{descore_thr}{
The threshold of DE score percentage used to determine FDmMGenes with relatively higher MSB score
}
  \item{enrich_input_directory}{
The parameter passed to STRINGdb to determine whether all the database files will be downloaded into this directory and the package can then be used off-line
}
  \item{version}{
The parameter passed to STRINGdb to determine the version used
}
  \item{species}{
The parameter passed to STRINGdb to determine the species to do enrichment analysis, in defalt it is human
}
  \item{bp_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched BP terms
}
  \item{kegg_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched KEGG pathways
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
funm6aviewer is the main function of the package which is used to viraulize the single base DmM sites, identify and visualize functional DmM genes in one button.
The permutime will influence the accuracy of identified FDmMGenes, more permutation times will generate more reliable result while takes longer time. In default, it is 10^4 and this is usually sufficient and comsumes less time, users can set it larger if the scale of DmM genes is larger.
}
\value{
By default, funm6aviewer will output results of

1. Functional DmM genes information;

2. DmM sites distribution on RNA;

3. Counts of DmM sites on different RNA regions;

4. DmM sites and reads coverage on interested gene;

5. Function enrichment of FDmMGenes;

6. DmMGene's MSB score along with DE score;

7. Network of FDmMGene's MSB neighbors.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

## example using default hg19 genome

dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)

bamreadsgr <- system.file("extdata", "bamgrlist_toy.RData", package="Funm6AViewer")
load(bamreadsgr)

siggene <- c("CCNT1", "MYC", "BCL2")
permutime <- 1000
enrich_input_directory <- system.file("extdata", package="Funm6AViewer")

re <- funm6aviewer(dminfo, deinfo, grlist, siggene =  siggene, permutime = permutime,
                   enrich_input_directory = enrich_input_directory)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

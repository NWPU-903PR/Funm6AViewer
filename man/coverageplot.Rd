\name{coverageplot}
\alias{coverageplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of reads coverage on gene.
}
\description{
Visualization of reads coverage on interested gene.
}
\usage{
coverageplot(dminfo,
             grlist,
             intrested_gene,
             track_height = NA,
             savepath = NA,
             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
             orgsymbol = org.Hs.egSYMBOL,
             orgdb = org.Hs.eg.db,
             zoom_region = NA,
             label_log2FD = FALSE,
             label_log2FD_rot = 45,
             add_guideline = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dminfo}{
A dataframe of DmM sites information. It can be generated using summarydmdeepm6A from DMDeepm6A package result.
}
  \item{grlist}{
A list containing the reads count of IP/Input bam files under different condition and the id of each files. It can be generated using  makegrreadsfrombam function. This is necessary if users want to watch the reads coverage for interested genes.
}
  \item{intrested_gene}{
A vector of interested gene to visualize, can be entrez gene id or gene symbol
}
  \item{track_height}{
A numeric number to set the max ylim of tracks showing reads count coverage on gene
}
  \item{savepath}{
The file path where to save the result
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if use other genomes instead of hg19.
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
  \item{orgdb}{
The annotation of genome. You need to input it similar to "org.Hs.eg.db" if DefaultGenome is use other genomes instead of human genome.
}
  \item{zoom_region}{
The region to be zoomed near the differential methylation sites. It should be a integer, e.g., 100 means 100 bp up- and down-stream of the sites will be zoomed.
}
  \item{label_log2FD}{
A logical value determines whether the log2 fold change should be labeled with differential methylation sites, FALSE in default.
}
  \item{label_log2FD_rot}{
The angle to rotate the text of log2 fold change.
}
  \item{add_guideline}{
A logical value determines whether to add guideline for each differential methylation site to better observe the fragments coverage near the sites, FALSE in default.
}
}
\details{
This function is used to visualize reads coverage and DmM sites on interested gene.
}
\value{
By default, coverageplot will output results of DmM sites and reads coverage on interested gene.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

## reads coverage plot for interested genes

dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)

bamreadsgr <- system.file("extdata", "bamgrlist_toy.RData", package="Funm6AViewer")
load(bamreadsgr)

siggene <- c("CCNT1", "MYC", "BCL2")

re <- coverageplot(dminfo = dminfo, grlist = grlist, intrested_gene = siggene)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

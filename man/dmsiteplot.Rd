\name{dmsiteplot}
\alias{dmsiteplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of DmM sites on gene's isoforms.
}
\description{
Visualization of DmM sites on gene's isoforms.
}
\usage{
dmsiteplot(dminfo,
           intrested_gene,
           savepath = NA,
           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
           orgsymbol = org.Hs.egSYMBOL,
           orgdb = org.Hs.eg.db)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dminfo}{
A dataframe of DmM sites information. It can be generated using summarydmdeepm6A from DMDeepm6A package result.
}
  \item{intrested_gene}{
A vector of interested gene to visualize, can be entrez gene id or gene symbol
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
}
\details{
This function is used to visualize DmM sites on interested gene's isoforms.
}
\value{
By default, dmsiteplot will output results of DmM sites on interested gene's isoforms.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

## dm sites plot for interested genes

dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)

siggene <- c("CCNT1", "MYC", "BCL2")

re <- dmsiteplot(dminfo = dminfo, intrested_gene = siggene)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

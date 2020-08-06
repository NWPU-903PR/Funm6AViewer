\name{siggenepathplot}
\alias{siggenepathplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of BP and KEGG pathway enrichment of functional DmM genes with 5'UTR, CDS and 3'UTR DmM sites.
}
\description{
Visualization of BP and KEGG pathway enrichment of functional DmM genes with 5'UTR, CDS and 3'UTR DmM sites.
}
\usage{
siggenepathplot(fdmgene,
                intrested_gene,
                bp_fdr_thr = 0.05,
                kegg_fdr_thr = 0.05,
                orgsymbol = org.Hs.egSYMBOL,
                input_directory = "",
                version = "10",
                species = 9606,
                savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdmgene}{
A vector of functional DmM genes, can be gene entrez id or gene symbol
}
  \item{intrested_gene}{
A vector of interested genes to visualize, can be gene entrez id or gene symbol
}
  \item{bp_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched BP terms
}
  \item{kegg_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched KEGG pathways
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
  \item{input_directory}{
The parameter passed to STRINGdb to determine whether all the database files will be downloaded into this directory and the package can then be used off-line
}
  \item{version}{
The parameter passed to STRINGdb to determine the version used
}
  \item{species}{
The parameter passed to STRINGdb to determine the species to do enrichment analysis, in defalt it is human
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
The enricnment analysis is done using STRINGdb.
}
\value{
By default, siggenepathplot will output results of

1. Plot of the context specific BP and KEGG pathway of interested genes;

2. xls formate of BP and KEGG pathway enrichment terms.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

## context specific function annotation for interested genes

dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
fdmgene <- unique(dminfo$name)

siggene <- c("CCNT1", "MYC", "BCL2")

# The input_directory is the parameter passed to string_db, the GO and KEGG function
# annotation data will be downloaded to the setted filepath.
# You can also donwloaded the annotation data previously
# from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
# and set the input_directory as where you save the data.
input_directory <- "E:/Funm6A_package/data"

re <- siggenepathplot(fdmgene = fdmgene, intrested_gene = siggene,
                      input_directory = input_directory)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

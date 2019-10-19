\name{enrichmentplot}
\alias{enrichmentplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of BP and KEGG pathway enrichment of functional DmM genes with 5'UTR, CDS and 3'UTR DmM sites.
}
\description{
Visualization of BP and KEGG pathway enrichment of functional DmM genes with 5'UTR, CDS and 3'UTR DmM sites.
}
\usage{
enrichmentplot(fdmgene,
               sigthr = 0.3,
               bp_fdr_thr = 0.05,
               kegg_fdr_thr = 0.05,
               input_directory = "",
               version = "10",
               species = 9606,
               savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdmgene}{
A dataframe of functional DmM genes' information. It is generated using fdmdeepm6A function.
}
  \item{sigthr}{
The threshold of the FDR used to determine whether a DmMGene is a FDmMGene
}
  \item{bp_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched BP terms
}
  \item{kegg_fdr_thr}{
The parameter passed to STRINGdb to determine the threshold of fdr to determine significant enriched KEGG pathways
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
By default, enrichmentplot will output results of

1. Heat map of BP and KEGG pathway enrichment of FDmMGenes;

2. xls formate of the heatmap information.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

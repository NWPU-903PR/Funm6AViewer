\name{siggenescoreplot}
\alias{siggenescoreplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of MSB scores of interested genes.
}
\description{
Visualization of MSB scores of interested genes.
}
\usage{
siggenescoreplot(fdmgene,
                 siggene,
                 plotname = "FocusedGene",
                 savepath = NA,
                 sigthresh = 0.3,
                 rescore_thr = 0.8,
                 descore_thr = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdmgene}{
A dataframe of functional DmM genes' information. It is generated using fdmdeepm6A function.
}
  \item{siggene}{
A vector of the gene symbol of interested genes to visualize
}
  \item{plotname}{
The name of the plot
}
  \item{savepath}{
The file path where to save the result
}
  \item{sigthresh}{
The threshold of the FDR used to determine whether a DmMGene is a FDmMGene
}
  \item{rescore_thr}{
The threshold of relative MSB score percentage used to determine FDmMGenes with relatively higher MSB score
}
  \item{descore_thr}{
The threshold of DE score percentage used to determine FDmMGenes with relatively higher MSB score
}
}
\details{
Visualization of MSB scores along with DE scores of interested genes and FDmMGenes with relatively higher MSB scores and DE scores.
}
\value{
By default, siggenescoreplot will output result of DmMGene’s MSB score along with DE score.
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

\name{getdescore}
\alias{getdescore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get DE score from DE information.
}
\description{
Get DE score from DE information.
}
\usage{
getdescore(deseqre,
           scoretype = "pval",
           orgsymbol = org.Hs.egSYMBOL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deseqre}{
A dataframe of DE information. It can be generated using getdeinfo from input bam files.
}
  \item{scoretype}{
The value used to calculate the DE score. It can be "pval" or "padj" denotes use the p-value or FDR of genes to calculate the DE score.
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
}
\details{
Get DE score from DE information.
}
\value{
By default, getdescore will output a named vector of DE scores. The name is the gene symbol.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

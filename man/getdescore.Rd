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
           scoretype = "pval")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deseqre}{
A dataframe of DE information. It can be generated using getdeinfo from input bam files.
}
  \item{scoretype}{
The value used to calculate the DE score. It can be "pval", "padj", or "log2FoldChange" denotes use the p-value, FDR or log2FoldChange of genes to calculate the DE score.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

\name{getdeinfo}
\alias{getdeinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get gene differential expression information from input samples.
}
\description{
Get gene differential expression information from input samples.
}
\usage{
getdeinfo(grlist,
          txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
          savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grlist}{
A list containing the reads count of IP/Input bam files under different condition and the id of each files. It can be generated using  makegrreadsfrombam function.
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if use other genomes instead of hg19.
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
The DE analysis is done using DESeq2.
}
\value{
By default, getdeinfo will output result of gene differential expresion information in both xls and data.frame which can used to calculate the DE score.
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

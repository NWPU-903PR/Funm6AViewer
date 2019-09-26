\name{makegrreadsfrombam}
\alias{makegrreadsfrombam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get reads count from bam files of MeRIP-Seq data.
}
\description{
This function will read the bam files to GRanges and save them with their condition id as RData.
}
\usage{
makegrreadsfrombam(IP_bams,
                   Input_bams,
                   condition,
                   minimal_alignment_MAPQ = 30,
                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP_bams}{
A vector of characters denote the file path of IP samples in bam formate
}
  \item{Input_bams}{
A vector of characters denote the file path of input samples in bam formate, should be paired with the IP samples
}
  \item{condition}{
A vector of characters denote the condition of IP and Input samples. Should be the same length as the ip_bams and the values should be "untreated" or "treated".
}
  \item{minimal_alignment_MAPQ}{
parameter passed to GenomicAlignments to contol the reads quality.
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if use other genomes instead of hg19.
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
This function will read the bam files to GRanges and save them with their condition id as RData.
}
\value{
By default, makegrreadsfrombam will output a list containing the reads count of IP/Input bam files under different condition and the id of each files.
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

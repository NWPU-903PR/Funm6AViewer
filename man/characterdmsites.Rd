\name{characterdmsites}
\alias{characterdmsites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of characteristics for all DmM sites
}
\description{
This function is used to visualize characteristics of all DmM sites.
}
\usage{
characterdmsites(dminfo,
                 txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 orgsymbol = org.Hs.egSYMBOL,
                 savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dminfo}{
A dataframe of DmM sites information. It can be generated using summarydmdeepm6A from DMDeepm6A package result.
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if use other genomes instead of hg19.
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
This function is used to visualize characteristics of all DmM sites.
}
\value{
The outputs of the function are:

1. DmM sites distribution on RNA;

2. Counts of DmM sites on different RNA regions;

3. Counts of hyper- and hypo- DmM sites on different RNA regions.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }% use one of  RShowDoc("KEYWORDS")

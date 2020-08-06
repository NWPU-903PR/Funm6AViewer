\name{fdmdeepm6A}
\alias{fdmdeepm6A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional DmM gene identification.
}
\description{
This function is used to identify functional DmM genes from a list of DmM genes using FunDMDeepm6A approach.
}
\usage{
fdmdeepm6A(dminfo,
           descore,
           top_alph = 0.8,
           fungenethr = 0.3,
           datapath = NA,
           UTR5only = FALSE,
           orgsymbol = org.Hs.egSYMBOL,
           savepath = NA,
           savename = "Funm6AGene",
           permutime = NA,
           no_cores = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dminfo}{
A dataframe of DmM sites information. It can be generated using summarydmdeepm6A from DMDeepm6A package result.
}
  \item{descore}{
A vector of named numbers denote the differntial expression score of all genes, it can be genarated using getdescore function
}
  \item{top_alph}{
The top percentage threshold used to aggregate the ranks from each PPI network. In defalt, the threshold is set as 0.8 which means only MSB scores larger than 80 pencentage of DE scores will contribute to the functional rank. A larger top_alph (must no more than 1) makes the functional ranks aggregation test more rigorous.
}
  \item{fungenethr}{
The threshold of the FDR used to determine whether a DmMGene is a FDmMGene
}
  \item{datapath}{
The file path where the network information required of FunDMDeepm6A, usually do not need to input, only if you prefer to use your own PPI networks
}
  \item{UTR5only}{
Whether the input DmM genes only harbor DmM sites on 5'UTR, it can be TRUE, FALSE or a vector of 1 and 0 which has the same lenght of DMgene to denote each DmM gene only harbor DmM sites on 5'UTR or not
}
  \item{orgsymbol}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is use other genomes instead of human genome.
}
  \item{savepath}{
The file path where to save the result
}
  \item{savename}{
The name of the xls formate result which contain the information of FDmMGenes
}
  \item{permutime}{
The permutation times to calculate an empirical p-value for DmMGenes
}
  \item{no_cores}{
The number of cores used to run the function in parallel
}
}
\details{
The permutime will influence the accuracy of identified FDmMGenes, more permutation times will generate more reliable result while takes longer time. In default, it is 100 times of the number of DmM genes and this is adequate but time consuming.
}
\value{
By default, fdmdeepm6A will output result as data.frame and save an xls formate file containg the MSB scores and the significance of a DmMGene to be a functional DmMGene.
}
\references{
Funm6AViewer: Visualization of single base differential m6A methylation sites and functional DmM genes.
}
\author{
Songyao Zhang
}
\examples{

## FunDMDeep-m6A

dminfo <- system.file("extdata", "DMinfo_toy.xls", package="Funm6AViewer")
deinfo <- system.file("extdata", "DEinfo_toy.xls", package="Funm6AViewer")

dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- read.delim(deinfo, header = TRUE, stringsAsFactors = FALSE)
descore <- getdescore(deinfo)

# The datapath is the filepath where the required PPI data saved, they can be downloaded
# from https://pan.baidu.com/s/1qOGG57OgxmrTwSbbBEeQ2w&shfl=sharepset
datapath <- "E:/Funm6A_package/data"
permutime <- 1000

re <- fdmdeepm6A(dminfo = dminfo, descore = descore, datapath = datapath, permutime = permutime)

# plot interested gene MSB score
siggene <- c("CCNT1", "MYC", "BCL2")
siggenescoreplot(fdmgene = re, siggene = siggene)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }

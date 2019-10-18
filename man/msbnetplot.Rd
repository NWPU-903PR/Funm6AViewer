\name{msbnetplot}
\alias{msbnetplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of MSB neighbors of interested genes.
}
\description{
Visualization of MSB neighbors of interested genes.
}
\usage{
msbnetplot(genesymbol,
           dmgene,
           descore,
           datapath = NA,
           savepath = NA,
           savename = NA,
           plotsize = NA,
           labeloff = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genesymbol}{
A vector of the gene symbol of interested genes to visualize
}
  \item{dmgene}{
A vector of the gene symbol of DmM genes
}
  \item{descore}{
A vector of named numbers denote the differntial expression score of all genes, it can be genarated using getdescore function
}
  \item{datapath}{
The file path where the network information required of FunDMDeepm6A, usually do not need to input, only if you prefer to use your own PPI networks
}
  \item{savepath}{
The file path where to save the result
}
  \item{savename}{
The name of the saved figure
}
  \item{plotsize}{
The size of the saved figure, by defalt, it is 6 when the number of MSB neighbors is no more than 50.
}
  \item{labeloff}{
Logical parameter that determine whether to display the gene name of MSB neighbors.
}
}
\details{
msbnetplot is used to Visualize of MSB neighbors of one or more interested genes.
}
\value{
By default, msbnetplot will output result of network of one or more FDmMGenes' MSB neighbors.
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

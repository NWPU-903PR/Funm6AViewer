\name{summarydmdeepm6A}
\alias{summarydmdeepm6A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the result of DMDeepm6A.
}
\description{
Summarize the result of DMDeepm6A.
}
\usage{
summarydmdeepm6A(dmpath,
                 sigthresh = 0.05,
                 savepath = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmpath}{
The filepath of DMDeepm6A package results.
}
  \item{sigthresh}{
The threshold used to determine whether a candidate DmM site is significant DmM site.
}
  \item{savepath}{
The file path where to save the result
}
}
\details{
This function is used to summary the result of DMDeepm6A package.
}
\value{
The output of this function is used as the input to Funm6AViewer package.
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

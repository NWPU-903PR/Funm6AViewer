
## get de score function

getdescore <- function(deseqre,
                       scoretype = "pval") {

  if (scoretype == "pval") {
    descore <- deseqre$pval
  } else {
    descore <- deseqre$padj
  }

  dename <- deseqre$name
  names(descore) <- as.character(dename)

  descore <- descore[!is.na(descore)]
  descore <- -log10(descore)
  cf <- quantile(descore, 0.99)
  descore[descore > cf] <- cf

  return(descore)
}

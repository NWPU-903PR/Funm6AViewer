
## get de score function

getdescore <- function(deseqre,
                       scoretype = "pval") {

  names(deseqre) <- c("name", "pval", "padj", "log2FoldChange")
  if (scoretype == "pval") {
    descore <- as.numeric(deseqre$pval)
  }
  if (scoretype == "padj") {
    descore <- as.numeric(deseqre$padj)
  }
  if (scoretype == "log2FoldChange") {
    descore <- as.numeric(deseqre$log2FoldChange)
  }

  dename <- deseqre$name
  names(descore) <- as.character(dename)

  descore[is.na(descore)] <- 1
  if (scoretype != "log2FoldChange") {
    descore <- -log10(descore)
  }
  # de_sig <- descore[descore >= -log10(0.05)]
  cf <- quantile(descore, 0.99)
  descore[descore > cf] <- cf

  return(descore)
}

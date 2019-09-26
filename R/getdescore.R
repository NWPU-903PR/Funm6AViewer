
## get de score function

getdescore <- function(deseqre,
                       scoretype = "pval",
                       orgsymbol = org.Hs.egSYMBOL) {

  if (scoretype == "pval") {
    descore <- deseqre$pval
  } else {
    descore <- deseqre$padj
  }

  dename <- .entrez2symbol(orgsymbol, deseqre$name)
  names(descore) <- as.character(dename)

  descore <- descore[!is.na(descore)]
  descore <- -log10(descore)
  cf <- quantile(descore, 0.99)
  descore[descore > cf] <- cf

  return(descore)
}


## get gene symbol

.entrez2symbol <- function(x, id) {

  mapped_genes <- mappedkeys(x)
  result <- as.list(x[mapped_genes])
  entrez_id <- as.numeric(names(result))  # entrez ID
  gene_symbol <- as.character(result)     # gene symbol
  ID_convert <- data.frame(entrez_id,gene_symbol)
  gene_symbol <- ID_convert$gene_symbol[match(id, ID_convert$entrez_id)]

  return(as.character(gene_symbol))
}

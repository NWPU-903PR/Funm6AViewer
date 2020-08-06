

getdeinfo <- function(grlist,
                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                      savepath = NA) {

  if (is.na(savepath)) { savepath <- getwd() }
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

  ## make genes
  gr <- exonsBy(txdb, by = "gene")
  # gr <- genes(txdb)
  # gr1 <- transcriptsBy(txdb, by = "gene")
  # gr0 <- gr1[!is.element(names(gr1), names(gr))]
  #
  # gr0 <- union(gr0, gr0)
  # ind <- names(gr0)
  # gr0 <- unlist(gr0)
  # ind <- match(ind, names(gr0))
  # gr0 <- gr0[ind]
  #
  # gr <- c(gr, gr0)
  # gr <- gr[names(gr1)]


  ## get reads
  len <- length(grlist)
  grid <- grlist[[len]]
  grlist <- grlist[-len]
  grlist <- grlist[c(grid$untreated_input, grid$treated_input)]

  ## get sample condition
  sample_condition <- c(rep("untreated", length(grid$untreated_input)), rep("treated", length(grid$treated_input)))

  ## count reads for genes
  se <- lapply(grlist, function(x, gr){countOverlaps(gr, x, ignore.strand = TRUE)}, gr)
  re <- matrix(unlist(se), ncol = length(se))
  rownames(re) <- names(se[[1]])
  re <- re[rowSums(re) > 10,]

  ## DE analysis
  condition <- factor(sample_condition, levels = c("untreated", "treated"))
  dds <- DESeqDataSetFromMatrix(re, DataFrame(condition), ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)

  resultsNames(dds)
  resLFC <- lfcShrink(dds, coef=2, type="apeglm")

  xls <- data.frame(name = row.names(re),
                    pval = res$pvalue,
                    padj = res$padj,
                    log2FoldChange = resLFC$log2FoldChange)

  write.table(xls, file =  paste(savepath, "DEGeneSummary.xls" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  return(xls)
}


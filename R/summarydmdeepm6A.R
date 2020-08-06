
## get m6A info from DMDeepm6A result

summarydmdeepm6A <- function(dmpath,
                             sigthresh = 0.05,
                             savepath = NA) {

  hyper <- paste(dmpath, "HyperMethySite.xls", sep = "/")
  hypo <- paste(dmpath, "HypoMethySite.xls", sep = "/")
  diff <- paste(dmpath, "CandidateDiffMethySite.xls", sep = "/")

  ## hyper
  hypergr <- .readxsltogr(hyper, "hyper")

  ## hypo
  hypogr <- .readxsltogr(hypo, "hypo")

  ## diff
  diffgr <- .readxsltogr(diff, "diff", sigthresh = sigthresh)

  dmgr <- rbind(hypergr, hypogr, diffgr)

  if (is.na(savepath)) {savepath <- getwd()}
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}
  write.table(dmgr, file =  paste(savepath, "DMDeepm6AResSummary.xls" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  return(dmgr)
}

.readxsltogr <- function(xls, xlstype, sigthresh = NA) {

  xlstable <- read.table(xls, header = T, stringsAsFactors = F)
  if (!is.na(sigthresh)) {
    xlstable <- xlstable[xlstable$diff.padj <= sigthresh,]
  }

  gr <- xlstable[,1:6]

  if (xlstype == "hyper") {
    gr$log2fd <- log2(xlstable$fold_enrchment)
  }

  if (xlstype == "hypo") {
    gr$log2fd <- -log2(xlstable$fold_enrchment)
  }

  if (xlstype == "diff") {
    gr$log2fd <- xlstable$diff.log2.OR
  }

  return(gr)
}


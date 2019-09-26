
enrichmentplot <- function(fdmgene,
                           sigthr = 0.3,
                           bp_fdr_thr = 0.05,
                           kegg_fdr_thr = 0.05,
                           input_directory = "",
                           version = "10",
                           species = 9606,
                           savepath = NA) {

  if(is.na(savepath)) {savepath <- getwd()}

  string_db <- STRINGdb$new( version = version, species = species,
                             score_threshold=0, input_directory = input_directory)

  fgutr3 <- fdmgene[fdmgene$padj <= sigthr & fdmgene$SitePosition == "UTR3",]
  fgutr5 <- fdmgene[fdmgene$padj <= sigthr & fdmgene$SitePosition == "UTR5",]
  fgcds <- fdmgene[fdmgene$padj <= sigthr & fdmgene$SitePosition == "CDS",]

  fgutr3 <- fgutr3[order(fgutr3$padj),]
  fgutr5 <- fgutr5[order(fgutr5$padj),]
  fgcds <- fgcds[order(fgcds$padj),]

  hits_fgutr3 <- string_db$map(fgutr3, "DMgene", removeUnmappedRows = TRUE)
  hits_fgutr3 <- hits_fgutr3$STRING_id

  hits_fgutr5 <- string_db$map(fgutr5, "DMgene", removeUnmappedRows = TRUE)
  hits_fgutr5 <- hits_fgutr5$STRING_id

  hits_fgcds <- string_db$map(fgcds, "DMgene", removeUnmappedRows = TRUE)
  hits_fgcds <- hits_fgcds$STRING_id

  e1 <- string_db$enrichment_heatmap( list( hits_fgutr3, hits_fgutr5, hits_fgcds),
                                      list("UTR3","UTR5", "CDS"),
                                      fdr_threshold = bp_fdr_thr, enrichmentType = "Process", title="BP",
                                      output_file = paste(savepath, "BP_enrichment.pdf", sep = "/"))

  e2 <- string_db$enrichment_heatmap( list( hits_fgutr3, hits_fgutr5, hits_fgcds),
                                      list("UTR3","UTR5", "CDS"),
                                      fdr_threshold = kegg_fdr_thr, enrichmentType = "KEGG", title="KEGG",
                                      output_file = paste(savepath, "KEGG_enrichment.pdf", sep = "/"))
  e1 <- cbind(rownames(e1), e1)
  colnames(e1)[1] <- "BP"
  write.table(e1, file =  paste(savepath, "BP_enrichment.xls", sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  e2 <- cbind(rownames(e2), e2)
  colnames(e2)[1] <- "KEGG"
  write.table(e2, file =  paste(savepath, "KEGG_enrichment.xls", sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  return(list(BP = e1, KEGG = e2))
}










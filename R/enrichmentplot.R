
enrichmentplot <- function(fdmgene,
                           sigthr = 0.3,
                           bp_fdr_thr = 0.05,
                           kegg_fdr_thr = 0.05,
                           input_directory = "",
                           version = "10",
                           species = 9606,
                           savepath = NA) {

  if(is.na(savepath)) {savepath <- getwd()}
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

  string_db <- STRINGdb$new( version = version, species = species,
                             score_threshold=0, input_directory = input_directory)

  if (is.element("SitePosition", names(fdmgene))) {

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
                                        list("UTR3","UTR5", "CDS"), iea = FALSE,
                                        fdr_threshold = bp_fdr_thr, enrichmentType = "Process", title="BP",
                                        output_file = paste(savepath, "BP_enrichment.pdf", sep = "/"))

    e2 <- string_db$enrichment_heatmap( list( hits_fgutr3, hits_fgutr5, hits_fgcds),
                                        list("UTR3","UTR5", "CDS"), iea = FALSE,
                                        fdr_threshold = kegg_fdr_thr, enrichmentType = "KEGG", title="KEGG",
                                        output_file = paste(savepath, "KEGG_enrichment.pdf", sep = "/"))
    e1 <- cbind(rownames(e1), e1)
    colnames(e1)[1] <- "BP"
    e2 <- cbind(rownames(e2), e2)
    colnames(e2)[1] <- "KEGG"

    .writeenrichment(hits_fgutr3, string_db,
                     bp_fdr_thr = bp_fdr_thr,
                     kegg_fdr_thr = kegg_fdr_thr,
                     BPsavename = paste(savepath, "BP_enrichment_UTR3fdmgene.xls", sep = "/"),
                     KEGGsavename = paste(savepath, "KEGG_enrichment_UTR3fdmgene.xls", sep = "/"))

    .writeenrichment(hits_fgutr5, string_db,
                     bp_fdr_thr = bp_fdr_thr,
                     kegg_fdr_thr = kegg_fdr_thr,
                     BPsavename = paste(savepath, "BP_enrichment_UTR5fdmgene.xls", sep = "/"),
                     KEGGsavename = paste(savepath, "KEGG_enrichment_UTR5fdmgene.xls", sep = "/"))

    .writeenrichment(hits_fgcds, string_db,
                     bp_fdr_thr = bp_fdr_thr,
                     kegg_fdr_thr = kegg_fdr_thr,
                     BPsavename = paste(savepath, "BP_enrichment_CDSfdmgene.xls", sep = "/"),
                     KEGGsavename = paste(savepath, "KEGG_enrichment_CDSfdmgene.xls", sep = "/"))

    re <- list(BP = e1, KEGG = e2)
  } else {

    fg <- fdmgene[fdmgene$padj <= sigthr,]

    fg <- fg[order(fg$padj),]

    hits<- string_db$map(fg, "DMgene", removeUnmappedRows = TRUE)
    hits <- hits$STRING_id

    .writeenrichment(hits, string_db,
                     bp_fdr_thr = bp_fdr_thr,
                     kegg_fdr_thr = kegg_fdr_thr,
                     BPsavename = paste(savepath, "BP_enrichment.xls", sep = "/"),
                     KEGGsavename = paste(savepath, "KEGG_enrichment.xls", sep = "/"))

    print(paste("BP and KEGG enrichment results are saved in BP_enrichment.xls and KEGG_enrichment.xls under folder", savepath))

    re <- fg

  }



  return(re)
}


.writeenrichment <- function(hits,
                             string_db,
                             bp_fdr_thr = 0.05,
                             kegg_fdr_thr = 0.05,
                             BPsavename = "BP_enrichment.xls",
                             KEGGsavename = "KEGG_enrichment.xls") {

  ## GO

  enrichmentGO <- string_db$get_enrichment(hits, category = "Process", methodMT = "fdr", iea = FALSE)
  enrichmentGO <- enrichmentGO[enrichmentGO$pvalue_fdr <= bp_fdr_thr,]

  gene_id <- string_db$get_term_proteins(term_ids = enrichmentGO$term_id, string_ids = hits, enableIEA = FALSE)
  ind <- tapply(gene_id$preferred_name, gene_id$term_id, function(x){paste(x, collapse = ", ")})
  enriched_genes <- ind[match(enrichmentGO$term_id, names(ind))]

  enrichmentGO <- cbind(enrichmentGO, enriched_genes)

  ## KEGG

  enrichmentKEGG <- string_db$get_enrichment(hits, category = "KEGG", methodMT = "fdr", iea = FALSE)
  enrichmentKEGG <- enrichmentKEGG[enrichmentKEGG$pvalue_fdr <= kegg_fdr_thr,]

  gene_id <- string_db$get_term_proteins(term_ids = enrichmentKEGG$term_id, string_ids = hits, enableIEA = FALSE)
  ind <- tapply(gene_id$preferred_name, gene_id$term_id, function(x){paste(x, collapse = ", ")})
  enriched_genes <- ind[match(enrichmentKEGG$term_id, names(ind))]

  enrichmentKEGG <- cbind(enrichmentKEGG, enriched_genes)

  ## save result

  write.table(enrichmentGO, file =  BPsavename,
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.table(enrichmentKEGG, file =  KEGGsavename,
              sep = "\t", row.names = FALSE, quote = FALSE)

}


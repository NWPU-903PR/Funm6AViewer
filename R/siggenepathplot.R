
siggenepathplot <- function(fdmgene,
                            intrested_gene,
                            bp_fdr_thr = 0.05,
                            kegg_fdr_thr = 0.05,
                            orgsymbol = org.Hs.egSYMBOL,
                            input_directory = "",
                            version = "10",
                            species = 9606,
                            savepath = NA) {

  fdmgene <- .getgenesymbol(orgsymbol, fdmgene)
  fdmgene <- fdmgene$genesymbol

  siggene <- .getgenesymbol(orgsymbol, intrested_gene)
  siggene <- siggene$genesymbol

  if(is.na(savepath)) {savepath <- getwd()}
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

  string_db <- STRINGdb$new( version = version, species = species,
                             score_threshold=0, input_directory = input_directory)

  fdmgene <- data.frame(fdmgene)
  hits_fg <- string_db$map(fdmgene, "fdmgene", removeUnmappedRows = TRUE)
  hits <- unique(hits_fg$STRING_id)

  ## enrichment
  enrichmentGO <- string_db$get_enrichment(hits, category = "Process", methodMT = "fdr", iea = FALSE)
  enrichmentGO <- enrichmentGO[enrichmentGO$pvalue_fdr <= bp_fdr_thr,]

  enrichmentKEGG <- string_db$get_enrichment(hits, category = "KEGG", methodMT = "fdr", iea = FALSE)
  enrichmentKEGG <- enrichmentKEGG[enrichmentKEGG$pvalue_fdr <= kegg_fdr_thr,]

  sighits <- hits_fg[is.element(hits_fg$fdmgene, siggene),]
  sighits <- unique(sighits)

  for (i in 1:length(siggene)) {
    genename <- sighits$fdmgene[i]

    enrichmentGO_sig <- string_db$get_enrichment(sighits$STRING_id[i], category = "Process", methodMT = "fdr", iea = FALSE)
    enrichmentGO_sig <- enrichmentGO[is.element(enrichmentGO$term_id, enrichmentGO_sig$term_id),]
    n <- nrow(enrichmentGO_sig)
    n <- min(n, 30)
    enrichmentGO_sig <- enrichmentGO_sig[1:n,]

    savename <- paste(savepath, "BP", sep = "/")
    if (n > 0) {
      .pathwayplot(genename, enrichmentGO_sig, savename)
    }


    enrichmentKEGG_sig <- string_db$get_enrichment(sighits$STRING_id[i], category = "KEGG", methodMT = "fdr", iea = FALSE)
    enrichmentKEGG_sig <- enrichmentKEGG[is.element(enrichmentKEGG$term_id, enrichmentKEGG_sig$term_id),]
    n <- nrow(enrichmentKEGG_sig)
    n <- min(n, 30)
    enrichmentKEGG_sig <- enrichmentKEGG_sig[1:n,]

    savename <- paste(savepath, "KEGG", sep = "/")
    if (n > 0) {
      .pathwayplot(genename, enrichmentKEGG_sig, savename)
    }
  }

  .writeenrichment(hits, string_db,
                   bp_fdr_thr = bp_fdr_thr,
                   kegg_fdr_thr = kegg_fdr_thr,
                   BPsavename = paste(savepath, "BP_enrichment_Allfdmgenes.xls", sep = "/"),
                   KEGGsavename = paste(savepath, "KEGG_enrichment_Allfdmgenes.xls", sep = "/"))

}




.pathwayplot <- function(genename, entichmentpathways, savename) {

  n <- nrow(entichmentpathways)
  ## plot net

  ed <- cbind(rep(genename, n), entichmentpathways$term_description, -log10(entichmentpathways$pvalue_fdr))
  ed <- data.frame(ed)
  names(ed) <- c("from", "to", "weight")

  nodename <- c(genename, entichmentpathways$term_description)
  nodesize <- -log2(entichmentpathways$pvalue_fdr)
  nodesize <- nodesize/max(nodesize)*10
  nodesize <- c(10, nodesize)

  dmnode <- data.frame(id = nodename, nodesize = nodesize)

  net <- graph_from_data_frame(d=ed, vertices=dmnode, directed=F)

  # Generate colors
  V(net)$color <- c("orange", rep("lightgreen", n))
  V(net)$frame.color <- NA

  # set node size
  V(net)$size <- V(net)$nodesize
  # E(net)$width <- as.numeric(E(net)$weight)/5


  l <- cbind(1:vcount(net), c(1, vcount(net):2))
  plot(net, vertex.label.font = 4, vertex.label.color = "gray30",
       vertex.label.cex = .7, edge.color = "gray75", layout = l)

  pdf(file = paste(savename, "_", genename, "_Function.pdf", sep = ""), width = 6, height = 6)
  plot(net, vertex.label.font = 4, vertex.label.color = "gray30",
       vertex.label.cex = .7, edge.color = "gray75", layout = l)
  dev.off()

}























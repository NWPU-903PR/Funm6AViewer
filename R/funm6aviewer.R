
funm6aviewer <- function(dminfo,
                         deinfo,
                         bamgrlist = NA,
                         siggene = NA,
                         entrezid = NA,
                         fungenethr = 0.3,
                         permutime = 10^4,
                         descoretype = "pval",
                         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         orgdb = org.Hs.eg.db,
                         orgsymbol = org.Hs.egSYMBOL,
                         datapath = NA,
                         rescore_thr = 0.8,
                         descore_thr = 0.8,
                         enrich_input_directory = "",
                         version = "10",
                         species = 9606,
                         bp_fdr_thr = 0.05,
                         kegg_fdr_thr = 0.05,
                         savepath = NA) {



  if (is.na(savepath)) { savepath <- getwd() }
  if (is.na(datapath)) { datapath <- system.file("extdata", package="Funm6AViewer") }
  if (!is.na(entrezid[1])) { siggene <- dminfo$genesymbol[match(entrezid, dminfo$name)] }

  descore <- getdescore(deinfo, scoretype = descoretype, orgsymbol = orgsymbol)

  ## plot character

  print("Visualizing DmM sites characteristics...")

  chara <- characterdmsites(dminfo, txdb = txdb, savepath = savepath)

  if (!is.na(siggene[1])) {
    savepath1 <- paste(savepath, "DmMSiteOnGene", sep = "/")
    dir.create(savepath1)
    ## plot dmsites
    for (i in 1:length(siggene)) {
      vp <- dmsiteplot(dminfo, genesymbol = siggene[i], savepath = savepath1, txdb = txdb, orgdb = orgdb)
    }

    ## plot coverage
    if (!is.na(bamgrlist[1])) {
      for (i in 1:length(siggene)) {
        vp <- coverageplot(dminfo, bamgrlist, genesymbol = siggene[i], savepath = savepath1, txdb = txdb, orgdb = orgdb)
      }
    }
  }

  ## FDMDeepm6A

  rm(bamgrlist)
  gc()

  ## get dmgene
  gene_utr3 <- unique(dminfo$genesymbol[dminfo$UTR3 != 0])
  gene_cds <- unique(dminfo$genesymbol[dminfo$CDS != 0])
  gene_utr5 <- unique(dminfo$genesymbol[dminfo$UTR5 != 0])
  dmgene <- list(utr3 = gene_utr3, utr5 = gene_utr5, cds = gene_cds)

  savepath2 <- paste(savepath, "FunDMDeepm6ARe", sep = "/")
  dir.create(savepath2)

  print("Running FunDMDeepm6A for genes with UTR3 DmM sites...")
  ## utr3
  DMgene <- dmgene$utr3
  fdm_utr3 <- fdmdeepm6A(DMgene, descore, datapath = datapath, savepath = savepath2, permutime = permutime,
                         savename = "Funm6AGene_utr3.xls")

  print("Running FunDMDeepm6A for genes with CDS DmM sites...")
  ## cds
  DMgene <- dmgene$cds
  fdm_cds <- fdmdeepm6A(DMgene, descore, datapath = datapath, savepath = savepath2, permutime = permutime,
                        savename = "Funm6AGene_cds.xls")

  print("Running FunDMDeepm6A for genes with UTR5 DmM sites...")
  ## utr5
  DMgene <- dmgene$utr5
  fdm_utr5 <- fdmdeepm6A(DMgene, descore, datapath = datapath, savepath = savepath2, permutime = permutime,
                         savename = "Funm6AGene_utr5.xls", UTR5only = TRUE)

  ## combine result
  fdm_utr3$UTR5only <- "UTR3"
  fdm_utr5$UTR5only <- "UTR5"
  fdm_cds$UTR5only <- "CDS"

  fdmgene <- rbind(fdm_utr5, fdm_cds, fdm_utr3)
  names(fdmgene)[7] <- "SitePosition"

  write.table(fdmgene, file =  paste(savepath2, "Funm6AGene.xls", sep = "/"), sep = "\t", row.names = FALSE, quote = FALSE)

  print("Visualizing functional DmM genes...")
  if (!is.na(siggene[1])) {
    savepath3 <- paste(savepath, "MSBScorePlot", sep = "/")
    dir.create(savepath3)

    ## plot MSB score
    siggenescoreplot(fdm_utr3, siggene, plotname = "UTR3_Gene", sigthresh = fungenethr, savepath = savepath3,
                     rescore_thr = rescore_thr, descore_thr = descore_thr)
    siggenescoreplot(fdm_utr5, siggene, plotname = "UTR5_Gene", sigthresh = fungenethr, savepath = savepath3,
                     rescore_thr = rescore_thr, descore_thr = descore_thr)
    siggenescoreplot(fdm_cds, siggene, plotname = "CDS_Gene", sigthresh = fungenethr, savepath = savepath3,
                     rescore_thr = rescore_thr, descore_thr = descore_thr)
  }

  if (!is.na(siggene[1])) {
    savepath4 <- paste(savepath, "MSBNetPlot", sep = "/")
    dir.create(savepath4)

    ## plot MSB net
    dmgene <- unique(as.character(fdmgene$DMgene))
    netl <- lapply(as.list(siggene), msbnetplot, dmgene, descore, datapath, savepath4)
    net <- msbnetplot(siggene, dmgene, descore, datapath, savepath = savepath4, savename = "SigGene", labeloff = T)
  }

  ## plot enrichment
  savepath5 <- paste(savepath, "FunctionEnrichmentPlot", sep = "/")
  dir.create(savepath5)

  funenrich <- enrichmentplot(fdmgene, sigthr = fungenethr, bp_fdr_thr = bp_fdr_thr, kegg_fdr_thr = kegg_fdr_thr,
                              input_directory = enrich_input_directory, savepath = savepath5)

  re <- list(dmgene = dminfo,
             descore = descore,
             fundmgene_utr3 = fdm_utr3,
             fundmgene_utr5 = fdm_utr5,
             fundmgene_cds = fdm_cds,
             fundmgene = fdmgene,
             funenrichment = funenrich)

  return(re)
}

















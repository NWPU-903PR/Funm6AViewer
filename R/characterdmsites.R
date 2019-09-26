

characterdmsites <- function(dmgr,
                             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             savepath = NA) {

  if (is.na(savepath)) {savepath <- getwd()}

  gc_txdb <- makeGuitarCoordsFromTxDb(txdb, noBins = 10)

  gr <-  GRanges(seqnames = dmgr$chr,
                 IRanges(start = dmgr$chromEnd, width = 1),
                 strand = dmgr$strand)
  mcols(gr) <- dmgr

  ## guitar plot
  gr_hyper <- gr[gr$foldenrich == "hyper"]
  gr_hypo <- gr[gr$foldenrich == "hypo"]

  saveToPDFprefix <- paste(savepath, "DmMSiteDistribution", sep = "/")
  grplot <- list(diffpeak = gr, hyper = gr_hyper, hypo = gr_hypo)
  GuitarPlot(grplot,
             GuitarCoordsFromTxDb = gc_txdb,
             saveToPDFprefix = saveToPDFprefix)


  ## pichart plot
  pidat <- data.frame(mcols(gr))

  df <- data.frame(group = factor(c("UTR3","CDS","UTR5", "LncRNA"), levels = c("UTR5", "CDS", "UTR3", "LncRNA")),
                   value = c(sum(pidat$UTR3), sum(pidat$CDS), sum(pidat$UTR5), sum(pidat$LncRNA)))

  bp <- ggplot(df,aes(x="",y=value,fill=group)) + geom_bar(stat="identity", width=1)

  pie1 <- bp + coord_polar("y", start = 0) + labs(x="") +
    scale_fill_brewer(palette="Dark2") +
    geom_text(aes(label = paste(value, "\n", scales::percent(value/sum(value)), sep = "")),
              size=4, position=position_stack(vjust = 0.5)) +
    theme_minimal() +
    theme(axis.title=element_blank(),
          axis.ticks=element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank())
  print(pie1)

  pdf(file = paste(savepath,  "SitePositionCounts.pdf", sep = "/"), width = 6, height = 4)
  print(pie1)
  dev.off()

  ## hyper hypo
  df <- data.frame(group = factor(rep(c("UTR3", "CDS", "UTR5", "LncRNA"), each = 2),
                                  levels = c("UTR3", "UTR5", "CDS", "LncRNA")),
                   value = c(sum(pidat$UTR3 == 1 & pidat$foldenrich == "hyper"),
                             sum(pidat$UTR3 == 1 & pidat$foldenrich == "hypo"),
                             sum(pidat$CDS == 1 & pidat$foldenrich == "hyper"),
                             sum(pidat$CDS == 1 & pidat$foldenrich == "hypo"),
                             sum(pidat$UTR5 == 1 & pidat$foldenrich == "hyper"),
                             sum(pidat$UTR5 == 1 & pidat$foldenrich == "hypo"),
                             sum(pidat$LncRNA == 1 & pidat$foldenrich == "hyper"),
                             sum(pidat$LncRNA == 1 & pidat$foldenrich == "hypo")),
                   fd = rep(c("hyper", "hypo"), 4))

  bp <- ggplot(df,aes(x=group, y=value, fill=fd)) + geom_bar(stat="identity", width=0.9)
  pie2 <- bp + coord_polar("x", start = 0) + scale_fill_brewer(palette="Dark2", direction = -1) +
    geom_text(aes(label = value), size = 2.5, position=position_stack(vjust = 0.8)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_line(colour = "lightblue", size = 1),
          legend.title = element_blank())
  print(pie2)

  pdf(file = paste(savepath,  "HyperHypoSitePositionCounts.pdf", sep = "/"), width = 6, height = 4)
  print(pie2)
  dev.off()

  return(list(SitePositionCounts = pie1, HyperHypoSitePositionCounts = pie2))
}




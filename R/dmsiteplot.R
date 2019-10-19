
## plot dm gene sites

dmsiteplot <- function(dminfo,
                       intrested_gene,
                       savepath = NA,
                       txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                       orgsymbol = org.Hs.egSYMBOL,
                       orgdb = org.Hs.eg.db) {

  ## prepare parameter

  if (ncol(dminfo) < 9) {
    genesymbol <- .getgenesymbol(orgsymbol, dminfo$name)
    genesymbol <- genesymbol$genesymbol
    dminfo$genesymbol <- genesymbol
    dminfo$foldenrich <- "hyper"
    dminfo$foldenrich[dminfo$log2fd < 0] <- "hypo"
    dminfo <- .getpeakposition(dminfo, txdb)
  }

  dmanno <- dminfo
  dmgr <-  GRanges(seqnames = dminfo$chr,
                   IRanges(start = dminfo$chromEnd, width = 1),
                   strand = dminfo$strand)
  mcols(dmgr) <- dmanno[, c(4, 5, 7, 8, 9)]

  for (i in 1:length(intrested_gene)) {
    siggene <- .getgenesymbol(orgsymbol, intrested_gene[i])
    sampleid <- siggene$entrezeid
    genesymbol <- siggene$genesymbol

    if (is.na(savepath)) {
      savepath <- getwd()
    }

    ## get site info
    examsite <- dmgr[which(dmgr$genesymbol == genesymbol)]

    ## get gene track
    gr <- genes(txdb, filter = list(gene_id = sampleid))
    grs <- transcriptsBy(txdb, by = "gene")
    grs <- grs[names(grs) == sampleid]
    grs <- unlist(grs)
    trs <- geneModelFromTxdb(txdb, orgdb,gr=gr)
    trs <- trs[is.element(names(trs), grs$tx_name)]

    if (length(trs) > 8) {
      width_trs <- lapply(trs, function(x){sum(width(x$dat))})
      width_trs <- unlist(width_trs)
      thr <- sort(width_trs, decreasing = T)[8]
      trs <- trs[width_trs >= thr]
    }

    ## get plot tracklist
    examsite$score <- abs(examsite$log2fd)
    examsite$color <- "#CC9999"
    examsite$color[examsite$foldenrich == "hypo"] <- "#666699"
    examsite$border <- "#999999"
    examsite$cex = 0.8
    examsite$label.parameter.rot <- 45
    examsite$feature.height <- 0

    lollipopData <- new("track", dat=examsite[examsite$foldenrich == "hyper"],
                        dat2=examsite[examsite$foldenrich == "hypo"], type="lollipopData")

    ## set plot parameter
    optSty <- optimizeStyle(trackList(trs, lollipopData))
    trackList <- optSty$tracks
    viewerStyle <- optSty$style

    ## set plot position
    setTrackViewerStyleParam(viewerStyle, "margin", c(.05, .12, .05, .1))

    ## set height and color
    len <- length(trs)
    for (i in 1:len) {
      setTrackStyleParam(trackList[[i]], "color", c("#D2691E"))
      setTrackStyleParam(trackList[[i]], "height",0.05)
    }
    setTrackStyleParam(trackList[[length(trackList)]], "height",(0.9 - length(trs)*0.05)*0.8)

    ## set name
    names(trackList) <- c(names(trs),"DmMSite")

    ## set y
    for (i in 1:(length(trs))) {
      setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=.6))
      setTrackStyleParam(trackList[[i]], "ylabpos", "upstream")
    }
    setTrackStyleParam(trackList[[length(trackList)]], "ylabgp", list(cex=0.8))

    if (!is.element("hypo", examsite$foldenrich)) {
      setTrackStyleParam(trackList[[length(trackList)]], "marginBottom",
                         (0.9 - length(trs)*0.05)*0.4 + length(trs)*0.05 +0.05)
    }

    ## plot
    vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle)

    addArrowMark(list(x=unit(-0.02, "npc"),
                      y=unit(length(trs)*0.05 + 0.08, "npc")),
                 label="Hypo",
                 angle = 10,
                 length = unit(0.12, "inches"),
                 col="#666699",
                 quadrant = 3,
                 vp=vp)

    addArrowMark(list(x=unit(-0.02, "npc"),
                      y=unit(0.88, "npc")),
                 label="Hyper",
                 angle = 10,
                 length = unit(0.12, "inches"),
                 col="#CC9999",
                 quadrant = 2,
                 vp=vp)

    pdf(file = paste(savepath, "/", genesymbol, "_DmMSites.pdf", sep = ""), width = 6, height = 4)

    viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle)

    addArrowMark(list(x=unit(-0.02, "npc"),
                      y=unit(length(trs)*0.05 + 0.08, "npc")),
                 label="Hypo",
                 angle = 10,
                 length = unit(0.12, "inches"),
                 col="#666699",
                 quadrant = 3,
                 vp=vp)

    addArrowMark(list(x=unit(-0.02, "npc"),
                      y=unit(0.88, "npc")),
                 label="Hyper",
                 angle = 10,
                 length = unit(0.12, "inches"),
                 col="#CC9999",
                 quadrant = 2,
                 vp=vp)

    dev.off()
  }

  return(vp)
}


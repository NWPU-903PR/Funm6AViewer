
## plot dm gene sites

coverageplot <- function(dmgr,
                         grlist,
                         sampleid = NA, 
                         genesymbol = NA, 
                         savepath = NA,
                         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         orgdb = org.Hs.eg.db) {
  
  ## prepare parameter
  
  grid <- grlist[[length(grlist)]]
  grlist <- grlist[-length(grlist)]
  
  dmanno <- dmgr
  dmgr <-  GRanges(seqnames = dmgr$chr,
                   IRanges(start = dmgr$chromEnd, width = 1),
                   strand = dmgr$strand)
  mcols(dmgr) <- dmanno[, c(4, 5, 7, 8, 9)]
  
  if (!is.na(genesymbol)) {
    sampleid <- dmgr$name[match(toupper(genesymbol), toupper(dmgr$genesymbol))]
  }
  
  if (!is.na(sampleid)) {
    genesymbol <- dmgr$genesymbol[match(sampleid, dmgr$name)]
  }
  
  if (is.na(savepath)) {
    savepath <- getwd()
  }
  
  ## get site info
  examsite <- dmgr[dmgr$name == sampleid]
  
  ## get gene track
  gr <- genes(txdb, filter = list(gene_id = sampleid))
  grs <- transcriptsBy(txdb, by = "gene")
  grs <- grs[names(grs) == sampleid]
  grs <- unlist(grs)
  trs <- geneModelFromTxdb(txdb, orgdb,gr=gr)
  trs <- trs[is.element(names(trs), grs$tx_name)]
  
  width_trs <- lapply(trs, function(x){sum(width(x$dat))})
  width_trs <- unlist(width_trs)
  trs <- trs[width_trs == max(width_trs)]
  
  ## get plot tracklist
  examsite$score <- examsite$score/5
  examsite$color <- "#CC9999"
  examsite$color[examsite$foldenrich == "hypo"] <- "#666699"
  examsite$border <- "#999999"
  examsite$cex = 0.8
  examsite$feature.height <- 0
  
  lollipopData <- new("track", dat=examsite[examsite$foldenrich == "hyper"],
                      dat2=examsite[examsite$foldenrich == "hypo"], type="lollipopData")
  
  
  ## sizefactor
  totalreads <- lapply(grlist, length)
  totalreads <- unlist(totalreads)
  
  sf <- totalreads/median(totalreads)
  
  ## get coverage
  grcl <- list()
  for (i in 1:length(grlist)) {
    grc <- grlist[[i]][countOverlaps(grlist[[i]], gr, ignore.strand = T) > 0]
    dat <- coverageGR(grc)
    dat$score <- dat$score/sf[i]
    grcl[[i]] <- dat
  }
  
  len <- length(grid$ip)
  ymax <- lapply(grcl, function(x){max(x$score)})
  ymax <- unlist(ymax)[1:len]
  ymax <- round(max(ymax)*1.05)
  
  ## make tack
  
  trackList <- trackList(trs, lollipopData)
  for (i in 1:len) {
    track0 <- new("track", dat=grcl[[i]], dat2 = grcl[[len + i]], format = "BED", type = "data")
    trackList[[2 + i]] <- track0
  }
  
  trname <- rep("treated", len)
  trname[grid$treated_ip] <- paste("Treated", 1:length(grid$treated_ip))
  trname[grid$untreated_ip] <- paste("Untreated", 1:length(grid$untreated_ip))
  
  ## set name
  names(trackList) <- c(genesymbol, "DmMSite", trname)
  
  
  ## set plot parameter
  optSty <- optimizeStyle(trackList)
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  
  ## set plot position
  setTrackViewerStyleParam(viewerStyle, "margin", c(.05, .12, .05, .08))
  
  ## set height and color
  setTrackStyleParam(trackList[[1]], "color", c("#D2691E"))
  setTrackStyleParam(trackList[[1]], "height", 0.08)
  setTrackStyleParam(trackList[[2]], "height", 0.18)
  
  for (i in 1:len) {
    setTrackStyleParam(trackList[[i + 2]], "color", c("red", "blue"))
    setTrackStyleParam(trackList[[i + 2]], "height", 0.6/len)
  }
  
  ## set y
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=1))
  setTrackStyleParam(trackList[[1]], "ylabpos", "upstream")
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1))
  
  for (i in 1:len) {
    setTrackStyleParam(trackList[[i + 2]], "ylabgp", list(cex=1))
    setTrackYaxisParam(trackList[[i + 2]], "main", FALSE)
    setTrackStyleParam(trackList[[i + 2]], "ylim", c(0, ymax))
  }
  
  if (!is.element("hypo", examsite$foldenrich)) {
    setTrackStyleParam(trackList[[2]], "marginBottom", 0.65)
  }
  
  trackList <- trackList[c(1, 3:(len + 2), 2)]
  
  ## plot
  vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")
  
  pdf(file = paste(savepath, "/", genesymbol, "_ReadsCoverage.pdf", sep = ""), width = 12, height = 8)
  vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")
  dev.off()
  
  return(vp)
}






## plot dm gene sites

coverageplot <- function(dminfo,
                         grlist,
                         intrested_gene,
                         track_height = NA,
                         savepath = NA,
                         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         orgsymbol = org.Hs.egSYMBOL,
                         orgdb = org.Hs.eg.db) {

  ## prepare parameter

  if (ncol(dminfo) < 9) {
    names(dminfo) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "log2fd")
    genesymbol <- .getgenesymbol(orgsymbol, dminfo$name)
    genesymbol <- genesymbol$genesymbol
    dminfo$genesymbol <- genesymbol
    dminfo$foldenrich <- "hyper"
    dminfo$foldenrich[dminfo$log2fd < 0] <- "hypo"
    dminfo <- .getpeakposition(dminfo, txdb)
  }

  grid <- grlist[[length(grlist)]]
  grlist <- grlist[-length(grlist)]

  dmanno <- dminfo
  dmgr <-  GRanges(seqnames = dminfo$chr,
                   IRanges(start = (dminfo$chromStart + 1), end = dminfo$chromEnd),
                   strand = dminfo$strand)

  mcols(dmgr) <- dmanno[, c(4, 5, 7, 8, 9)]

  for (h in 1:length(intrested_gene)) {

    siggene <- .getgenesymbol(orgsymbol, intrested_gene[h])
    sampleid <- siggene$entrezeid
    genesymbol <- siggene$genesymbol

    if (is.na(savepath)) {
      savepath <- getwd()
    }
    if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

    ## get gene track
    gr <- genes(txdb, filter = list(gene_id = sampleid))
    grs <- transcriptsBy(txdb, by = "gene")
    grs <- grs[names(grs) == sampleid]
    grs <- unlist(grs)

    if ((length(gr) == 0) & (length(grs) != 0)) {

      ind <- as.character(seqnames(grs))
      grs <- grs[seqnames(grs) == ind[1]]

      ind <- as.character(strand(grs))
      grs <- grs[strand(grs) == ind[1]]

      gr <- union(grs, grs)
      gr$gene_id <- sampleid

    }

    ## get site info
    examsite <- dmgr[which(countOverlaps(dmgr, gr, ignore.strand = TRUE) != 0)]

    if (length(gr) != 0) {

      ## make transcripts track
      trs <- geneModelFromTxdb(txdb, orgdb,gr=gr)
      trs <- trs[is.element(names(trs), grs$tx_name)]

      width_trs <- lapply(trs, function(x){sum(width(x$dat))})
      width_trs <- unlist(width_trs)
      trs <- trs[width_trs == max(width_trs)]

      if (length(examsite) == 0) {
        examsite <- gr
        nodmpeak <- TRUE
      } else {
        nodmpeak <- FALSE
      }

      if (width(examsite)[1] == 1) {

        ## get plot tracklist
        examsite$score <- abs(examsite$log2fd)
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
          grc <- grlist[[i]]
          ind <- countOverlaps(grc, gr, ignore.strand = T) != 0
          grc <- grc[ind]
          dat <- coverageGR(grc)
          dat$score <- dat$score/sf[i]
          grcl[[i]] <- dat
        }

        ## make track ymax
        len <- length(grid$ip)
        ymax <- lapply(grcl, function(x){max(x$score)})
        ymax <- unlist(ymax)[1:len]
        ymax <- max(c(20, round(max(ymax)*1.05)))

        if (!is.na(track_height)) {ymax <- as.numeric(track_height)}

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

        trlun <- grid$untreated_ip + 2
        trltr <- grid$treated_ip + 2
        trackList <- trackList[c(1, trlun, trltr, 2)]

        ## plot
        vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")

        pdf(file = paste(savepath, "/", genesymbol, "_ReadsCoverage.pdf", sep = ""), width = 12, height = 8)
        vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")
        dev.off()

      } else {

        if (nodmpeak) {
          examsite$score <- 0.01
          strand(examsite) <- "*"
        } else {
          examsite$score[examsite$foldenrich == "hyper"] <- 1
          examsite$score[examsite$foldenrich == "hypo"] <- -1
          strand(examsite) <- "*"
        }

        lollipopData <- new("track", dat=examsite, format = "BED", type = "data")

        ## sizefactor
        totalreads <- lapply(grlist, length)
        totalreads <- unlist(totalreads)

        sf <- totalreads/median(totalreads)

        ## get coverage
        grcl <- list()
        for (i in 1:length(grlist)) {
          grc <- grlist[[i]]
          ind <- countOverlaps(grc, gr, ignore.strand = T) != 0
          grc <- grc[ind]
          dat <- coverageGR(grc)
          dat$score <- dat$score/sf[i]
          grcl[[i]] <- dat
        }

        ## make track ymax
        len <- length(grid$ip)
        ymax <- lapply(grcl, function(x){max(x$score)})
        ymax <- unlist(ymax)[1:len]
        ymax <- max(c(20, round(max(ymax)*1.05)))

        if (!is.na(track_height)) {ymax <- as.numeric(track_height)}

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
        names(trackList) <- c(genesymbol, "DmMPeak", trname)


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
        setTrackStyleParam(trackList[[2]], "ylim", c(-8, 8))
        setTrackYaxisParam(trackList[[2]], "draw", FALSE)

        for (i in 1:len) {
          setTrackStyleParam(trackList[[i + 2]], "ylabgp", list(cex=1))
          setTrackYaxisParam(trackList[[i + 2]], "main", FALSE)
          setTrackStyleParam(trackList[[i + 2]], "ylim", c(0, ymax))
        }

        ## set x
        if (nodmpeak) {
          setTrackXscaleParam(trackList[[2]], "label", "No DM peak detected")
          setTrackXscaleParam(trackList[[2]], "draw", TRUE)
        }

        trlun <- grid$untreated_ip + 2
        trltr <- grid$treated_ip + 2
        trackList <- trackList[c(1, trlun, trltr, 2)]

        ## plot
        vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")

        if (!nodmpeak) {
          ind <- resize(examsite, width = 1, fix = "center")
          indx <- (start(ind) - start(gr))/width(gr)
          indy <- rep(0.78, length(ind))
          indy[ind$foldenrich == "hypo"] <- 0.75
          indl <- as.character(ind$foldenrich)
          ind <- ind$foldenrich == "hyper"

          if (sum(ind) > 0) {
            addArrowMark(list(x=unit(indx[ind], "npc"),
                              y=unit(indy[ind], "npc")),
                         label = indl[ind],
                         angle = 10,
                         length = unit(0.12, "inches"),
                         col = "#CC9999",
                         quadrant = 2,
                         vp=vp)
          }

          if (sum(!ind) > 0) {
            addArrowMark(list(x=unit(indx[!ind], "npc"),
                              y=unit(indy[!ind], "npc")),
                         label = indl[!ind],
                         angle = 10,
                         length = unit(0.12, "inches"),
                         col = "#666699",
                         quadrant = 3,
                         vp=vp)
          }
        }


        ## save plot
        pdf(file = paste(savepath, "/", genesymbol, "_ReadsCoverage.pdf", sep = ""), width = 12, height = 8)
        vp <- viewTracks(trackList, gr=gr, autoOptimizeStyle=TRUE, newpage=TRUE, viewerStyle=viewerStyle, operator="+")

        if (!nodmpeak) {
          ind <- resize(examsite, width = 1, fix = "center")
          indx <- (start(ind) - start(gr))/width(gr)
          indy <- rep(0.78, length(ind))
          indy[ind$foldenrich == "hypo"] <- 0.75
          indl <- as.character(ind$foldenrich)
          ind <- ind$foldenrich == "hyper"

          if (sum(ind) > 0) {
            addArrowMark(list(x=unit(indx[ind], "npc"),
                              y=unit(indy[ind], "npc")),
                         label = indl[ind],
                         angle = 10,
                         length = unit(0.12, "inches"),
                         col = "#CC9999",
                         quadrant = 2,
                         vp=vp)
          }

          if (sum(!ind) > 0) {
            addArrowMark(list(x=unit(indx[!ind], "npc"),
                              y=unit(indy[!ind], "npc")),
                         label = indl[!ind],
                         angle = 10,
                         length = unit(0.12, "inches"),
                         col = "#666699",
                         quadrant = 3,
                         vp=vp)
          }
        }
        dev.off()

      }

    } else {
      print(paste(intrested_gene[h], "is not official gene symbol or entrez gene id! Please input either of them!"))
    }

  }

}





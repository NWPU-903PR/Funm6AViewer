# make grlist of reads from bam files

makegrreadsfrombam <- function(IP_bams,
                               Input_bams,
                               condition,
                               minimal_alignment_MAPQ = 30,
                               fragment_length = 100,
                               txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                               savepath = NA) {

  if (is.na(savepath)) { savepath <- getwd() }
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

  bams <- c(IP_bams, Input_bams)
  trid <- which(condition == "treated")
  unid <- which(condition == "untreated")
  sampleid <- list(ip = 1:length(IP_bams), input = 1:length(IP_bams) + length(IP_bams),
                   treated = c(trid, trid + length(IP_bams)), untreated = c(unid, unid + length(IP_bams)),
                   treated_ip = trid, treated_input = trid + length(IP_bams),
                   untreated_ip = unid, untreated_input = unid + length(IP_bams))

  grlist <- list()
  for (i in 1:length(bams)) {

    bamfile = bams[i]
    print(paste("Processing bam file", i))

    # prepare bam parameters
    what <- c("rname","strand", "pos","mapq","qwidth")
    param <- ScanBamParam(what=what)

    # read bam file
    gr0 <- readGAlignments(bamfile, param = param)
    ba <- mcols(gr0)
    gr0 <- granges(gr0)

    # MAPQ filter
    ba$mapq[which(is.na(ba$mapq))] <- 255

    ind <- which(ba$mapq > minimal_alignment_MAPQ)
    gr0 <- gr0[ind]
    ba <- ba[ind,]

    id_filter <- (!is.na(ba$rname)) & (!is.na(ba$pos)) & (!is.na(ba$strand))
    gr0 <- gr0[id_filter]

    rm(ba)
    gc()

    # shift
    gr <- resize(gr0, width = fragment_length, fix = "start", ignore.strand = F)

    rm(gr0)
    gc()

    ## get transcripts
    genegr <- transcripts(txdb)

    len <- width(genegr)
    tx <- resize(genegr, width = len + 2000, fix = "center")

    # unimaped reads
    gr <- gr[countOverlaps(gr, tx, ignore.strand = T) > 0]

    grlist[[i]] <- gr
  }

  grlist[[i+1]] <- sampleid

  save(grlist, file = paste(savepath, "bamsgrlist.RData", sep = "/"))

  return(grlist)

}




















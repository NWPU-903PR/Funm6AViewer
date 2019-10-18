# make grlist of reads from bam files

makegrreadsfrombam <- function(IP_bams,
                               Input_bams,
                               condition,
                               minimal_alignment_MAPQ = 30,
                               txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                               savepath = NA) {

  if (is.na(savepath)) { savepath <- getwd() }

  bams <- c(IP_bams, Input_bams)
  trid <- which(condition == "treated")
  unid <- which(condition == "untreated")
  sampleid <- list(ip = 1:length(IP_bams), input = 1:length(IP_bams) + length(IP_bams),
                   treated = c(trid, trid + length(IP_bams)), untreated = c(unid, unid + length(IP_bams)),
                   treated_ip = trid, treated_input = trid + length(IP_bams),
                   untreated_ip = unid, untreated_input = unid + length(IP_bams))

  grlist <- list()
  for (i in 1:length(bams)) {

    file = bams[i]
    print(paste("Processing bam file", i))

    # prepare bam parameters
    what <- c("rname","strand", "pos","mapq","qwidth")
    param <- ScanBamParam(what=what)

    # read bam file
    ba <- scanBam(file, param=param)
    ba <- ba[[1]]

    # MAPQ filter
    ba$mapq[which(is.na(ba$mapq))] <- 255
    ba$rname <- ba$rname[ba$mapq > minimal_alignment_MAPQ]
    ba$strand <- ba$strand[ba$mapq > minimal_alignment_MAPQ]
    ba$pos <- ba$pos[ba$mapq > minimal_alignment_MAPQ]
    ba$qwidth <- ba$qwidth[ba$mapq > minimal_alignment_MAPQ]
    ba$mapq <- ba$mapq[ba$mapq > minimal_alignment_MAPQ]

    gr <- GRanges(seqnames = ba$rname,
                  ranges = IRanges(start=ba$pos, width = ba$qwidth),
                  strand = ba$strand)

    genegr <- genes(txdb)

    # unimaped reads
    gr <- gr[countOverlaps(gr, genegr, ignore.strand = T) > 0]

    grlist[[i]] <- gr
  }

  grlist[[i+1]] <- sampleid

  save(grlist, file = paste(savepath, "bamsgrlist.RData", sep = "/"))

  return(grlist)

}




















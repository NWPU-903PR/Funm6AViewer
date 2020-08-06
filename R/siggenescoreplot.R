
siggenescoreplot <- function(fdmgene,
                             siggene,
                             plotname = "FocusedGene",
                             savepath = NA,
                             sigthresh = 0.3,
                             rescore_thr = 0.8,
                             descore_thr = 0.8) {

  if (is.na(savepath)) {savepath <- getwd()}
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}

  msbscore <- fdmgene$MSBscore
  descore <- fdmgene$DEscore
  rescore <- fdmgene$Relative_MSBscore

  ind <- is.element(fdmgene$DMgene, siggene)

  cf <- quantile(rescore, rescore_thr)
  cf0 <- quantile(descore, descore_thr)
  fdind <- (rescore >= cf) & (fdmgene$padj <= sigthresh) & (descore > cf0)

  group <- rep("OtherGene", length(ind))
  group[fdind] <- "RelativeSigGene"
  group[ind] <- "FocusedGene"
  shape <- fdmgene$padj <= sigthresh

  Gene <- fdmgene$DMgene
  DEScore <- as.numeric(descore)
  MSBScore <- as.numeric(msbscore)
  Significant <- ind
  GeneType <- factor(group, levels = c("FocusedGene", "RelativeSigGene", "OtherGene"))
  FDmMGene <- shape

  dat <- data.frame(Gene = Gene,
                    DEScore = DEScore,
                    MSBScore = MSBScore,
                    Significant = Significant,
                    fdind = fdind,
                    GeneType = GeneType,
                    FDmMGene = FDmMGene)
  dat <- rbind(dat[!ind,], dat[ind,])
  dat <- dat[dat$MSBScore!=0,]

  dat <- dat[order(dat$MSBScore, decreasing = T),]
  dat <- dat[!duplicated(dat$Gene),]

  p <- ggplot(dat, aes(x = DEScore, y = MSBScore)) +
    geom_point(aes(color = GeneType, shape = FDmMGene), size = 2) +
    scale_color_manual(values = c("red", "blue", "grey")) +
    theme_bw(base_size = 14) + theme(legend.position = "bottom") +
    ggtitle(plotname) +
    geom_text_repel(
      data = subset(dat, Significant| fdind  == TRUE),
      aes(label = Gene),
      size = 4,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))

  pdf(file = paste(savepath, "/", plotname, "_MSBScore.pdf", sep = ""), width = 12, height = 10)
  print(p)
  dev.off()

  return(p)
}

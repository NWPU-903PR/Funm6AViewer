
## fdmdeepm6A

fdmdeepm6A <- function(DMgene,
                       descore,
                       datapath = NA,
                       UTR5only = FALSE,
                       orgsymbol = org.Hs.egSYMBOL,
                       savepath = NA,
                       savename = "Funm6AGene.xls",
                       permutime = 100*length(DMgene),
                       no_cores = NA) {

  if (is.na(datapath)) { datapath <- system.file("extdata", package="Funm6AViewer") }

  DMgene <- .getgenesymbol(orgsymbol, DMgene)
  DMgene <- DMgene$genesymbol
  DMgene <- DMgene[!is.na(DMgene)]

  if (sum(UTR5only) == 0) {UTR5only <- rep(0, length(DMgene))}
  if (sum(UTR5only) == 1) {UTR5only <- rep(1, length(DMgene))}

  degenes <- .getgenesymbol(orgsymbol, names(descore))
  names(descore) <- degenes$genesymbol

  load(paste(datapath, "net.RData", sep = "/"))
  load(paste(datapath, "netmotif.RData", sep = "/"))

  rankpath <- c(paste(datapath, "hint_rankall.RData", sep = "/"),
                paste(datapath, "biogrid_rankall.RData", sep = "/"),
                paste(datapath, "iRef_rankall.RData", sep = "/"),
                paste(datapath, "multinet_rankall.RData", sep = "/"))

  ## get MSB
  dmdescore <- list()
  for (i in 1:4) {
    print(names(net)[i])
    load(rankpath[i])
    netid <- net[[i]]
    motifid <- sigmotif[[i]]
    rank <- as.matrix(rank)

    dmdescore0 <- .DEScore(DMgene, descore, netid, motifid, rank, UTR5only)
    dmdescore[[i]] <- dmdescore0

    # dat <- data.frame(x = c(dmdescore0$DEScore, dmdescore0$RawDEScore),
    #                   group = rep(c("DEScore", "RawDEScore"), each = nrow(dmdescore0)))
    # print(ggplot(dat,aes(x=x, fill=group, color = group)) + geom_density(alpha=0.25))
  }

  names(dmdescore) <- names(net)

  ## get all dmde score
  msbscore <- descore[match(toupper(DMgene), toupper(names(descore)))]
  names(msbscore) <- DMgene
  msbscore <- .getmsbscore(msbscore, dmdescore)

  ## get rank score

  msbscoresf <- colSums(msbscore)
  standard_library_size <- exp(mean(log(msbscoresf)))
  msbscoresf <- msbscoresf/standard_library_size
  msbscore <- t(t(msbscore)/msbscoresf)

  print(paste(permutime, "times random test for ranks, this may take a few hours..."))
  if (!is.na(no_cores)) {print(paste(no_cores, "cores are used."))}

  rankscore <- .getrankscore(msbscore)

  alphapath <- system.file("extdata", "alphascore.R", package="Funm6AViewer")
  rankpval <- .rankrandomtest(msbscore, permutime, no_cores, alphapath)

  rankfdr <- p.adjust(rankpval, method = "BH")

  dat <- data.frame(x = c(rankpval, rankfdr),
                    group = rep(c("pval", "padj"), each = length(rankpval)))
  # print(ggplot(dat,aes(x=x, fill=group, color = group)) + geom_density(alpha=0.25))

  rmsbscore <- rowMeans(msbscore)

  id <- match(toupper(DMgene), toupper(names(descore)))
  genedescore <- descore[id]
  genedescore[is.na(genedescore)] <- 0

  relativemsbscore <- rmsbscore - genedescore*(1 - UTR5only)

  xls <- data.frame(DMgene = DMgene,
                    MSBscore = rmsbscore,
                    DEscore = genedescore,
                    Relative_MSBscore = relativemsbscore,
                    pval = rankpval,
                    padj = rankfdr)

  colnames(msbscore) <- paste(names(net), "Score", sep = ".")
  xls <- cbind(xls, UTR5only, msbscore)

  if (is.na(savepath)) {savepath <- getwd()}
  write.table(xls, file =  paste(savepath, savename, sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  return(xls)
}

## DEScore

.DEScore <- function(DMgene, descore, netid, motifid, rank, UTR5only) {

  DMind <- is.element(toupper(netid$gene_name), toupper(DMgene))
  Funind <- is.element(toupper(netid$gene_name), toupper(netid$pathGene))

  DMind <- which(DMind == 1)
  Funind <- which(Funind == 1)

  id <- match(toupper(netid$gene_name), toupper(names(descore)))
  genedescore <- descore[id]
  genedescore[is.na(genedescore)] <- 0

  id <- match(toupper(netid$gene_name)[DMind], toupper(DMgene))
  UTR5only <- UTR5only[id]

  rank <- rank[,DMind]
  for (i in 1:length(DMind)) {
    if (UTR5only[i] == 1) {rank[DMind[i],i] <- 0} else {rank[DMind[i],i] <- 1}
    rank[,i] <- rank[,i]*genedescore
  }

  motifind <- lapply(motifid, function(x,y){sum(is.element(x,y))}, DMind)
  motifid <- motifid[unlist(motifind) != 0]

  len <- unlist(lapply(motifid, length))
  motifind <- lapply(motifid, function(x, y, z){sum(is.element(x,y) & is.element(x,z))}, DMind, Funind)
  motifind <- unlist(motifind) != len
  motifid <- motifid[motifind]

  len <- unlist(lapply(motifid, length))
  motifind <- unlist(motifid)
  id <- rep(1:length(len), times = len)

  dmmotif <- lapply(as.list(DMind), function(x, id, motifind) {unique(id[motifind == x])}, id, motifind)
  dmmotif <- lapply(dmmotif, function(x, motifid) {unique(unlist(motifid[x]))}, motifid)

  len <- unlist(lapply(dmmotif, length))
  dmmotif[len == 0] <- DMind[len == 0]

  dmdescore <- numeric()
  rdmdescore <- dmdescore
  for (i in 1:length(dmmotif)) {
    dmdescore[i] <- sum(rank[dmmotif[[i]],i])
    rdmdescore[i] <- rank[DMind[i],i]
  }

  dmdescore <- data.frame(DEScore = dmdescore, RawDEScore = rdmdescore)
  row.names(dmdescore) <- netid$gene_name[DMind]

  return(dmdescore)
}


## MSB score
.getmsbscore <- function(msbscore, dmdescore) {

  len <- length(msbscore)
  msbscorem <- matrix(rep(msbscore, 4), nrow = len, ncol = 4)
  for (i in 1:4) {
    msbscore0 <- dmdescore[[i]]
    ind <- match(toupper(row.names(msbscore0)), toupper(names(msbscore)))
    rankm0 <- msbscorem[,i]
    rankm0[ind] <- msbscore0$DEScore
    rankm0[is.na(rankm0)] <- 0
    msbscorem[,i] <- rankm0
  }

  row.names(msbscorem) <- names(msbscore)

  return(msbscorem)
}

.alphascore <- function(x, alph) {
  r <- x[1:4]
  s <- x[5:8]
  sor <- sort(r, index.return = TRUE)
  r <- r[sor$ix]
  s <- s[sor$ix]

  bscore <- betaScores(r)
  if (sum(s >= alph) == 0) {
    bscore <- 1
  } else {
    bscore <- min(bscore[s >= alph])
  }

  return(bscore)
}

.getrankscore <- function(msbscore) {

  len <- nrow(msbscore)
  rankm <- msbscore
  for (i in 1:4) {
    rankm0 <- rankm[,i]
    x <- order(rankm0, decreasing = TRUE)
    rankm0[x] <- (1:len)/len
    rankm[,i] <- rankm0
  }

  rankm <- cbind(rankm, msbscore)
  bscore <- apply(rankm, 1, .alphascore, -log10(0.05))
  names(bscore) <- row.names(msbscore)

  return(bscore)
}


## random test for rank

.rankrandomtest <- function(msbscore, m, no_cores, alphapath) {

  rankscore0 <- .getrankscore(msbscore)
  len <- length(rankscore0)

  if (!is.na(no_cores)) {
    cl <- makeCluster(no_cores)
    rmscore <- parLapply(cl, 1:m, .getrandomscore, msbscore, alphapath)
    stopCluster(cl)
  } else {
    rmscore <- lapply(1:m, .getrandomscore, msbscore, alphapath)
  }

  rmscorex <- matrix(unlist(rmscore), nrow = len)
  rankscore <- cbind(rankscore0, rmscorex)

  pval <- apply(rankscore, 1, function(x, m){sum(x[2:(m-1)] <= x[1])/m}, m)

  return(pval)
}

.getrandomscore <- function(x, msbscore, alphapath) {

  source(alphapath)
  len <- nrow(msbscore)

  randomscorem <- as.vector(msbscore)
  randomscorem <- sample(randomscorem, length(randomscorem))
  randomscorem <- matrix(randomscorem, nrow = len, ncol = 4)

  rankm <- randomscorem
  for (i in 1:4) {
    rankm0 <- rankm[,i]
    x <- order(rankm0, decreasing = TRUE)
    rankm0[x] <- (1:len)/len
    rankm[,i] <- rankm0
  }

  rankm <- cbind(rankm, randomscorem)
  bscore <- apply(rankm, 1, .alphascorex, -log10(0.05))
  names(bscore) <- row.names(msbscore)

  return(bscore)
}








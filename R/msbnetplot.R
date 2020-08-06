

msbnetplot <- function(genesymbol,
                       dmgene,
                       descore,
                       orgsymbol = org.Hs.egSYMBOL,
                       datapath = NA,
                       savepath = NA,
                       savename = NA,
                       plotsize = NA,
                       labeloff = FALSE) {

  descore <- abs(descore)

  genesymbol <- .getgenesymbol(orgsymbol, genesymbol)
  genesymbol <- genesymbol$genesymbol

  dmgene <- .getgenesymbol(orgsymbol, dmgene)
  dmgene <- dmgene$genesymbol
  dmgene <- dmgene[!is.na(dmgene)]

  degenes <- .getgenesymbol(orgsymbol, names(descore))
  names(descore) <- degenes$genesymbol

  if (is.na(datapath)) { datapath <- system.file("extdata", package="Funm6AViewer") }
  if (is.na(savepath)) {savepath <- getwd()}
  if (!dir.exists(savepath)) {dir.create(savepath, recursive = T)}
  if (is.na(savename)) {savename <- genesymbol[1]}

  net <- NULL
  sigmotif <- NULL

  load(paste(datapath, "net.RData", sep = "/"))
  load(paste(datapath, "netmotif.RData", sep = "/"))

  dmnodel <- character()
  dmedgel <- data.frame()
  pathind <- character()
  for (i in 1:4) {
    net0 <- net[[i]]
    sigmotif0 <- sigmotif[[i]]
    pathind <- c(pathind, net0$pathGene)

    symboleid <- match(toupper(genesymbol), toupper(net0$gene_name))

    len <- unlist(lapply(sigmotif0, length))
    motifind <- unlist(sigmotif0)
    id <- rep(1:length(len), times = len)

    dmmotif <- lapply(as.list(symboleid), function(x, id, motifind) {unique(id[motifind == x])}, id, motifind)
    dmmotif <- lapply(dmmotif, function(x, motifid) {unique(unlist(motifid[x]))}, sigmotif0)

    dmnode <- unique(unlist(dmmotif))
    dmedge <- net0$idlist0
    dmedge <- dmedge[is.element(dmedge[,1], dmnode) & is.element(dmedge[,2], dmnode),]

    dmnode <- net0$gene_name[dmnode]
    dmedge <- cbind(net0$gene_name[dmedge[,1]], net0$gene_name[dmedge[,2]])
    dmnodel <- c(dmnodel, dmnode)
    dmedgel <- rbind(dmedgel, dmedge)
  }
  pathind <- unique(pathind)

  dmnode <- unique(dmnodel)

  if (length(dmnode) != 0) {

    len <- length(dmnode)
    ag <- matrix(0, nrow = len, ncol=len)
    ed <- cbind(match(dmedgel[,1], dmnode), match(dmedgel[,2], dmnode))
    for(h in 1:nrow(ed)) {
      ag[ed[h,1],ed[h,2]] <- 1 + ag[ed[h,1],ed[h,2]]
      ag[ed[h,2],ed[h,1]] <- 1 + ag[ed[h,2],ed[h,1]]
    }

    rownames(ag) <- dmnode
    s <- rowSums(ag) != 0
    ag <- ag[s,s]
    ag[lower.tri(ag)] <- 0

    dmnode <- rownames(ag)
    edind <- which(ag!=0, arr.ind = T)
    w <- ag[edind]
    ed <- cbind(dmnode[edind[,1]], dmnode[edind[,2]], w)
    ed <- data.frame(ed)
    names(ed) <- c("from", "to", "weight")

    nodesize <- as.numeric(descore[match(toupper(dmnode), toupper(names(descore)))])
    nodesize[is.na(nodesize)] <- min(nodesize, na.rm = T)
    nodesize[nodesize == 0] <- max(nodesize) + 1
    nodesize[nodesize == max(nodesize)] <- min(nodesize)

    nodetype <- rep("Other", length(dmnode))
    nodetype[is.element(toupper(dmnode), toupper(pathind))] <- "PathGene"
    nodetype[is.element(toupper(dmnode), toupper(dmgene))] <- "DmMGene"
    nodetype[is.element(toupper(dmnode), toupper(genesymbol))] <- "FunDmMGene"

    dmnode <- data.frame(id = dmnode, nodesize = nodesize, nodetype = nodetype)

    ## generat igraph
    net <- graph_from_data_frame(d=ed, vertices=dmnode, directed=F)

    # Generate colors based on node type:
    colrs <- c("tomato", "gray70", "lightgreen", "orange")
    V(net)$color <- colrs[factor(V(net)$nodetype, levels = c("FunDmMGene", "Other", "PathGene", "DmMGene"))]
    V(net)$frame.color <- NA


    # set node size based on nodesize:
    nodesize <- V(net)$nodesize/min(V(net)$nodesize) + 1
    V(net)$size <- log(nodesize)/max(log(nodesize))*10

    vertexlabelfont <- 4
    if (labeloff | nrow(dmnode) > 150) {
      V(net)$label <- rep(NA, length(V(net)$name))
      V(net)$label[V(net)$color == "tomato"] <-  V(net)$name[V(net)$color == "tomato"]
      vertexlabelfont <- 16
    }

    # Set edge width based on weight:
    E(net)$width <- as.numeric(E(net)$weight)/2

    # We can even set the network layout:
    l <- layout_with_fr(net)
    plot(net, vertex.label.font = vertexlabelfont, vertex.label.color = "gray30",
         vertex.label.cex = .7, edge.color = "gray75", layout = l)
    legend(x=-1, y=-1.1, c("FocusedGene", "OtherGene", "PathGene", "DmMGene"), pch=21,
           col = NA, pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=2)

    if (is.na(plotsize)) {
      if (nrow(dmnode) <= 50) {plotsize <- 6}
      if (nrow(dmnode) > 50) {plotsize <- 8}
      if (nrow(dmnode) > 100) {plotsize <- 10}
    }

    pdf(file = paste(savepath, "/", savename, "_MSBNet.pdf", sep = ""), width = plotsize, height = plotsize)
    plot(net, vertex.label.font = vertexlabelfont, vertex.label.color = "gray30",
         vertex.label.cex = .7, edge.color = "gray75", layout = l)
    legend(x=-1, y=-1.1, c("FocusedGene", "OtherGene", "PathGene", "DmMGene"), pch=21,
           col = NA, pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=2)
    dev.off()

    return(net)
  }

}


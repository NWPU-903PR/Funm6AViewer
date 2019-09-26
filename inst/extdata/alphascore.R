
.alphascorex <- function(x, alph) {

  library(RobustRankAggreg)
  r <- x[1:4]
  s <- x[5:8]
  sor <- sort(r, index.return = T)
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

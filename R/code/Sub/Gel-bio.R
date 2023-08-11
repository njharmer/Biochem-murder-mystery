# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021

# Define function to provide an SDS-PAGE gel clue; mark is the marker sizes; sus is a matrix of the suspects' protein sizes
# The gel is giving separation in a range from 250 kDa to 3.2 kDa over a total of 1000 pixels.
# Gel equation used is log10(size in kDa) = 2.5 -(0.5*RBPB).
# Band intensity reflects quantity (currently at two levels, could be increased).
# A random amount of smiling is produced that will cause a change of up to 10% at the extremes.
# Band loading is also dodgy in some places.
# Returns a gel image.
# Requires that the elliptic package has been called.

Gel <- function(mark,sus,names=""){
  wells <- length(sus[,1])+1
  gel <- matrix(rep(0,(110000*wells)), nrow=100*wells)
  smile <- sample(2:10,1)
  for (i in 1:wells) {
    gel[(30+((i-1)*100)):(80+((i-1)*100)),1020:1100] <- 1
  }
  load <- round(rnorm(1,11,5)); if (load < 1){load <- 1}
  width <- round(rnorm(1,80,10))
  for (i in 1:length(mark)) {
    DBPB <- 1000-round(500*(2.5-(log10(mark[i]))))
    for (j in load:(load+width-1)){
      DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-j)/(100*wells/2)))/100)))
      for (k in (DBPBi-8):(DBPBi+8)){
        gel[j,k] <- limit(gel[j,k] + 0.015*(limit(8-(abs(k-DBPBi)),upper=6))*
                            (limit((width/2)-(abs(load+(width/2)-j)),10)),1)
      }
    }
  }
  for (i in (1:(wells-1))) {
    load <- round(rnorm(1,(11+(100*i)),5))
    width <- round(rnorm(1,80,10)); if (load+width-1 > 100*wells){load <- (100*wells)-width}
    for (j in 2:(length(sus[1,]))) {
      DBPB <- 1000-round(500*(2.5-log10(sus[i,j])))
      for (k in load:(load+width-1)){
        DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-k)/(100*wells/2)))/100)))
        for (l in (DBPBi-8):(DBPBi+8)){
          gel[k,l] <- limit(gel[k,l] + 0.015*(limit(8-(abs(l-DBPBi)),upper=6))*
                              (limit((width/2)-(abs(load+(width/2)-k)),10)),1)
        }
      }
    }
  }
  
  # Write out an image and return it as the output
  lab <- rep("",wells) ; lab[1] <- "markers"
  for (i in 2:(wells)) { lab[i] <- names[sus[(i-1),1]]}
  image(gel, axes=F, col = hcl.colors(n=10, palette = "Blues 3", rev = TRUE))
  mtext(text=lab, side=1, line=0.3, at=seq((0.5/wells),(1-(0.5/wells)),(1/wells)), las=2, cex=1)
  output <- recordPlot()
}

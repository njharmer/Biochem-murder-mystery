# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021, 2022

# Define function to provide an SDS-PAGE gel clue; mark is the marker sizes; sus is a matrix of the suspects' protein sizes
# The gel is giving separation in a range from 250 kDa to 3.2 kDa over a total of 100 pixels.
# Gel equation used is log10(size in kDa) = 2.5 -(0.5*RBPB)
# Returns a gel image.
Gel <- function(mark,sus,names=""){
  wells <- length(sus[,1])+1
  gel <- matrix(rep(0,(1100*wells)), nrow=10*wells)
  for (i in 1:wells) {
    gel[(3+((i-1)*10)):(8+((i-1)*10)),102:110] <- 1
  }
  for (i in 1:length(mark)) {
    DBPB <- 100-round(50*(2.5-(log10(mark[i]))))
    gel[2:9,(DBPB-1):(DBPB+1)] <- 0.5; gel[3:8,DBPB] <- 1
  }
  for (i in (1:(wells-1))) {
    for (j in 2:(length(sus[1,]))) {
      DBPB <- 100-round(50*(2.5-log10(sus[i,j])))
      gel[(2+(10*i)):(9+(10*i)),(DBPB-1):(DBPB+1)] <- 0.5; gel[(3+(10*i)):(8+(10*i)),DBPB] <- 1
    }
  }
  
  # Write out an image and return it as the output
  lab <- rep("",wells) ; lab[1] <- "markers"; lab[wells] <- "Impostor"
  for (i in 2:(wells-1)) { lab[i] <- names[sus[(i-1),1]]}
  image(gel, axes=F, col = hcl.colors(n=10, palette = "Blues 3", rev = TRUE))
  mtext(text=lab, side=1, line=0.3, at=seq((0.5/wells),(1-(0.5/wells)),(1/wells)), las=2, cex=1)
  output <- recordPlot()
}

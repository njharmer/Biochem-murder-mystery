# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021

# Define function to provide a western blot clue; mark is the marker sizes; sus is a matrix of the suspects' protein sizes
# probe is either "J" or "K".
# The second and fifth markers will light up with probe J; the first and fourth markers will light up with probe K.
# Only the largest suspect protein will light up (so no additional information is given about suspects).
# The gel is giving separation in a range from 250 kDa to 3.2 kDa over a total of 100 pixels.
# Gel equation used is log10(size in kDa) = 2.5 -(0.5*RBPB)
West <- function(mark,wsus,probe){
  wells <- length(wsus[,1])+1
  wes <- matrix(rep(0,(1100*wells)), nrow=wells*10)
  for (i in 1:wells) {
    wes[(3+((i-1)*10)):(8+((i-1)*10)),102:110] <- 1
  }
  if(probe == "J"){
    DBPB <- 100-round(50*(2.5-(log10(mark[2])))); wes[2:9,(DBPB-1):(DBPB+1)] <- 0.5; wes[3:8,DBPB] <- 1
    DBPB <- 100-round(50*(2.5-(log10(mark[5])))); wes[2:9,(DBPB-1):(DBPB+1)] <- 0.5; wes[3:8,DBPB] <- 1
    DBPB <- 100-round(50*(2.5-(log10(Prots[1]))))
    for (i in 1:(wells-1)) {
      if(substr(wsus[i,2],1,1) == "J") {
        wes[(2+(10*i)):(9+(10*i)),(DBPB-1):(DBPB+1)] <- 0.5; wes[(3+(10*i)):(8+(10*i)),DBPB] <- 1
      }
    }
  }

  if(probe == "K"){
    DBPB <- 100-round(50*(2.5-(log10(mark[1])))); wes[2:9,(DBPB-1):(DBPB+1)] <- 0.5; wes[3:8,DBPB] <- 1
    DBPB <- 100-round(50*(2.5-(log10(mark[4])))); wes[2:9,(DBPB-1):(DBPB+1)] <- 0.5; wes[3:8,DBPB] <- 1
    DBPB <- 100-round(50*(2.5-(log10(Prots[1]))))
    for (i in 1:(wells-1)) {
      if(substr(wsus[i,2],nchar(wsus[i,2]),nchar(wsus[i,2])) == "K") {
        wes[(2+(10*i)):(9+(10*i)),(DBPB-1):(DBPB+1)] <- 0.5; wes[(3+(10*i)):(8+(10*i)),DBPB] <- 1
      }
    }
  }
  # Write out an image and return it as the output
  lab <- rep("",wells) ; lab[1] <- "markers"; lab[wells] <- "Impostor"
  for (i in 2:(wells-1)) { lab[i] <- wsus[(i-1),1]}
  image(wes, axes=F, col = hcl.colors(n=10, palette = "Blue-Yellow", rev = TRUE), xlab = paste0("Probe ",probe))
  mtext(text=lab, side=1, line=0.3, at=seq((0.5/wells),(1-(0.5/wells)),(1/wells)), las=2, cex=1.2)
  output <- recordPlot()
  
}

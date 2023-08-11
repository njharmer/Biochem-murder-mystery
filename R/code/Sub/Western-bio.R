# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2022

# Define function to provide a western blot clue; mark is the marker sizes; sus is a matrix of the suspects' protein sizes
# probe is either "J" or "K".
# The second and fifth markers will light up with probe J; the first and fourth markers will light up with probe K.
# Only the largest suspect protein will light up (so no additional information is given about suspects).
# The gel is giving separation in a range from 250 kDa to 3.2 kDa over a total of 100 pixels.
# Gel equation used is log10(size in kDa) = 2.5 -(0.5*RBPB)
# Band intensity reflects quantity (currently at two levels, could be increased).
# A random amount of smiling is produced that will cause a change of up to 10% at the extremes.
# Band loading is also dodgy in some places.
# Returns a gel image.
# Requires that the elliptic package has been called.

West <- function(mark,wsus,probe="J"){
  wells <- length(wsus[,1])+1
  wes <- matrix(rep(0,(110000*wells)), nrow=wells*100)
  smile <- sample(2:10,1)
  for (i in 1:wells) {
    wes[(30+((i-1)*100)):(80+((i-1)*100)),1020:1100] <- 1
  }
  if(probe == "J"){
    markj <- c(mark[2],mark[5])
    load <- round(rnorm(1,11,5)); if (load < 1){load <- 1}
    width <- round(rnorm(1,80,10))
    for (i in 1:2){
      DBPB <- 1000-round(500*(2.5-(log10(markj[i])))) 
      for (j in load:(load+width-1)){
        DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-j)/(100*wells/2)))/100)))
        for (k in (DBPBi-8):(DBPBi+8)){
          wes[j,k] <- limit(wes[j,k] + 0.015*(limit(8-(abs(k-DBPBi)),upper=6))*
                              (limit((width/2)-(abs(load+(width/2)-j)),10)),1)
        }
      }
    }  
    
    
    DBPB <- 1000-round(500*(2.5-(log10(Prots[1]))))
    for (i in 1:(wells-1)) {
      if(substr(wsus[i,2],1,1) == "J") {
        load <- round(rnorm(1,(11+(100*i)),5))
        width <- round(rnorm(1,80,10)); if (load+width-1 > 100*wells){load <- (100*wells)-width}
        for (j in load:(load+width-1)){
          DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-j)/(100*wells/2)))/100)))
          for (k in (DBPBi-8):(DBPBi+8)){
            wes[j,k] <- limit(wes[j,k] + 0.015*(limit(8-(abs(k-DBPBi)),upper=6))*
                                (limit((width/2)-(abs(load+(width/2)-j)),10)),1)
          }
        }  
      }
    }
  }
  
  if(probe == "K"){
    markk <- c(mark[1],mark[4])
    load <- round(rnorm(1,11,5)); if (load < 1){load <- 1}
    width <- round(rnorm(1,80,10))
    for (i in 1:2){
      DBPB <- 1000-round(500*(2.5-(log10(markk[i])))) 
      for (j in load:(load+width-1)){
        DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-j)/(100*wells/2)))/100)))
        for (k in (DBPBi-8):(DBPBi+8)){
          wes[j,k] <- limit(wes[j,k] + 0.015*(limit(8-(abs(k-DBPBi)),upper=6))*
                              (limit((width/2)-(abs(load+(width/2)-j)),10)),1)
        }
      }
    }  
    
    DBPB <- 1000-round(500*(2.5-(log10(Prots[1]))))
    for (i in 1:(wells-1)) {
      if(substr(wsus[i,2],nchar(wsus[i,2]),nchar(wsus[i,2])) == "K") {
        load <- round(rnorm(1,(11+(100*i)),5))
        width <- round(rnorm(1,80,10)); if (load+width-1 > 100*wells){load <- (100*wells)-width}
        for (j in load:(load+width-1)){
          DBPBi <- DBPB + round(DBPB*((smile*(1-cos(((100*wells/2)-j)/(100*wells/2)))/100)))
          for (k in (DBPBi-8):(DBPBi+8)){
            wes[j,k] <- limit(wes[j,k] + 0.015*(limit(8-(abs(k-DBPBi)),upper=6))*
                                (limit((width/2)-(abs(load+(width/2)-j)),10)),1)
          }
        }  
      }
    }
  }
  # Write out an image and return it as the output
  lab <- rep("",wells) ; lab[1] <- "markers"
  for (i in 2:(wells)) { lab[i] <- wsus[(i-1),1]}
  image(wes, axes=F, col = hcl.colors(n=10, palette = "Blue-Yellow", rev = TRUE), xlab = paste0("Probe ",probe))
  mtext(text=lab, side=1, line=0.3, at=seq((0.5/wells),(1-(0.5/wells)),(1/wells)), las=2, cex=1)
  output <- recordPlot()
  
}

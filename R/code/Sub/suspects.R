# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021

# Define function to generate suspects
# x is the number of suspects to generate, tot is the number of Sprites, vic is the list of victims.
# returns a list of suspects
suspects <- function(x,tot,vic) {
  unique <- FALSE
  while (unique == FALSE) {
    unique <- TRUE
    sus <- sample(1:tot,x)
    for (i in 1:x) {
      for (j in 1:length(vic)) {
        if (sus[i] == vic[j]) {unique <- FALSE}
      }
    }
  }
  suspects <- sus[sort.list(sus)]
}

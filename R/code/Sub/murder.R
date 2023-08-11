# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021

# Define function for victim killing; tot is the number of sprites,
# vic is the list of victims; imp is the impostor no.
murder <- function(tot, vic, imp){
  unique <- FALSE
  while (unique == FALSE) {
    victim <- sample(1:tot,1)
    if (victim != imp) {unique <- TRUE}
    for (i in 1:length(vic)) {
      if (victim == vic[i]) {unique <- FALSE}
    }
  }
  output <- victim
}

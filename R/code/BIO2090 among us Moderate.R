# Script to produce data for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021-2023

# Install packages if necessary
# To run this script, it is first necessary to install the following four packages:
# tidyverse, rcartocolor, elliptic, kableExtra
# These have dependencies that will be automatically installed.
# On the first use of these scripts, uncomment the following line to install the necessary packages.
# install.packages(c('tidyverse','rcartocolor','elliptic','kableExtra'))

# Load packages
library(tidyverse)
library(rcartocolor)
library(elliptic)
library(kableExtra)
setwd("C:/Users/njh209/OneDrive - University of Exeter/Documents/R/Among Us")

# Source functions
source("code/Sub/murder.R")
source("code/Sub/suspects.R")
source("code/Sub/Gel.R")
source("code/Sub/Western.R")
source("code/Sub/HPLC.R")

# Set up core data
Sprites <- data.frame(Names = c("Blackcurrant","Bubblegum","Burgundy","Cappuccino","Cocoa","Dalmation","Lemon","Lime","Nitrile",
                                "Ocean","Orange","Peach"),
                      chem1 = c(0),
                      chem2 = c(0),
                      chem3 = c(0),
                      prot1 = c(1),
                      prot2 = c(0),
                      prot3 = c(0),
                      glyco = c(0))
Glyco <- c("O", "J", "JK", "JL", "K", "JLK")
lectins <- matrix(c(0,0,1,0,1,1,1,0,0,1,1,1),nrow=2)
MS <- c(1,2,3,3,2,4)
PepMass <- 836; hexmass <- 162.1
Prots <- c(75,55,40,30,25)
Markers <- c(210, 105,78,55,45,34,17,16,7,4)
HPLC <- data.frame(chems = c("Benzene", "3-BP", "2,4-DBP", "2,4,6-TBP", "2B-4NP"),
                   logP = c(2.1,2.6,3.2,4.4,2.4),
                   D0 = c(-0.0677,-0.0803,-0.0161,-0.0138,-0.1565),
                   D1 = c(0.0216,0.0312,0.0113,0.0079,0.0589),
                   c0 = c(4.57,4.734,6.722,7.931,4.6),
                   c1 = c(-9.331,-11.58,-15.89,-17.53,-12.46),
                   c2 = c(3.602,5.33,8.399,9.416,6.315),
                   h0 = c(0.0962,0.2041,0.1420,0.2658,0.3055),
                   h1 = c(-0.0112, -0.0358, -0.0141,-0.0171,-0.0707))
victims <- c(0,0,0)
cluesize <- 6 # The final  clue is always the impostor sample so this needs to be greater than 1!
Error <- 4 # Chance that the impostor will have tampered with a sample; range 0-100, whole numbers only.
safe = carto_pal(12,"Safe") # Set up a colourblind friendly colour palette

# Assign features to Sprites; features are name(1); three chemicals (2-4); three proteins (5-7); one glyco group(8).
for (i in 1:length(Sprites$Names)) {
  Sprites[i,2] <- sample(1:5,1)
  unique <- FALSE
  while (unique == FALSE) {
    Sprites[i,3] <- sample(1:5,1)
    Sprites[i,4] <- sample(1:5,1)
    if ((Sprites[i,4] != Sprites[i,2]) && (Sprites[i,4] != Sprites[i,3]) && (Sprites[i,2] != Sprites[i,3])) {unique <- TRUE}
  }
  Sprites[i,6] <- sample(2:5,1); Sprites[i,7] <- sample(2:5,1); Sprites[i,8] <- sample(1:6,1)
}

# Generate an impostor
Impostor <- sample(1:length(Sprites$Names),1)

# Carry out the first murder

victims[1] <- murder(length(Sprites$Names),victims,Impostor)

# Obtain a gel clue. Includes the impostor sample. Defined probability that samples are mixed up.

susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
gelsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
for (i in 1:cluesize) {
  gelsus[i,1] <- susp[i] ; gelsus[i,2] <- Prots[1]; gelsus[i,3] <- Prots[Sprites$prot2[susp[i]]]
  gelsus[i,4] <- Prots[Sprites$prot3[susp[i]]]
  if(sample(0:99,1)<Error){gelsus[i,3] <- Prots[sample(1:5,1)]}
  if(sample(0:99,1)<Error){gelsus[i,4] <- Prots[sample(1:5,1)]}
}
png(filename="output/Gel1.png")
gel <- Gel(Markers, gelsus, names=Sprites$Names)
dev.off()

# Obtain a western clue. Includes the impostor sample. Defined probability that samples are mixed up.
  
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
wessus <- matrix(rep(0,cluesize*2), nrow=cluesize)
for (i in 1:cluesize) {
  wessus[i,1] <- Sprites$Names[susp[i]] ; wessus[i,2] <- Glyco[Sprites$glyco[susp[i]]]
  if(sample(0:99,1)<Error){wessus[i,2] <- Glyco[sample(1:6,1)]}
}
if(sample(1:2,1) == 1){probe <- "J"} else {probe <- "K"}
png(filename="output/West1.png")
west <- West(Markers,wessus,probe)
dev.off()  

# Obtain a mass spectrometry clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
mssus <- matrix(rep(0,cluesize*2), nrow=cluesize); colnames(mssus) <- c("Sprite","m/z")
for (i in 1:cluesize) {
  mssus[i,1] <- Sprites$Names[susp[i]] ; mssus[i,2] <- PepMass + (hexmass * MS[Sprites$glyco[susp[i]]])
  if(sample(0:99,1)<Error){mssus[i,2] <- PepMass + (hexmass * MS[sample(1:6,1)])}
}
mssus[cluesize,1] <- "Impostor"
mssus1 <- as.data.frame(mssus)

# Obtain an HPLC clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
hpsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
names <- rep(" ", cluesize)
for (i in 1:cluesize) {
  hpsus[i,1] <- susp[i] ; hpsus[i,2] <- Sprites$chem1[susp[i]]; hpsus[i,3] <- Sprites$chem2[susp[i]]
  hpsus[i,4] <- Sprites$chem3[susp[i]]
  if(sample(0:99,1)<Error){hpsus[i,2] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,3] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,4] <- sample(1:5,1)}
  names[i] <- Sprites[susp[i],1]; 
}
names[cluesize] <- "X: impostor"
phi1 <- 0.4 + (sample(0:4,1)*0.05)
hp1 <- hplc(hpsus,HPLC,names,phi1)


# Carry out the second murder

victims[2] <- murder(length(Sprites$Names),victims,Impostor)

# Obtain a gel clue. Includes the impostor sample. Defined probability that samples are mixed up.

susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
gelsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
for (i in 1:cluesize) {
  gelsus[i,1] <- susp[i] ; gelsus[i,2] <- Prots[1]; gelsus[i,3] <- Prots[Sprites$prot2[susp[i]]]
  gelsus[i,4] <- Prots[Sprites$prot3[susp[i]]]
  if(sample(0:99,1)<Error){gelsus[i,3] <- Prots[sample(1:5,1)]}
  if(sample(0:99,1)<Error){gelsus[i,4] <- Prots[sample(1:5,1)]}
}
png(filename="output/Gel2.png")
gel <- Gel(Markers, gelsus, names=Sprites$Names)
dev.off()  

# Obtain a western clue. Includes the impostor sample. Defined probability that samples are mixed up.

susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
wessus <- matrix(rep(0,cluesize*2), nrow=cluesize)
for (i in 1:cluesize) {
  wessus[i,1] <- Sprites$Names[susp[i]] ; wessus[i,2] <- Glyco[Sprites$glyco[susp[i]]]
  if(sample(0:99,1)<Error){wessus[i,2] <- Glyco[sample(1:6,1)]}
}
if(probe == "K"){probe <- "J"} else {probe <- "K"}
png(filename="output/West2.png")
west <- West(Markers,wessus,probe)
dev.off()  

# Obtain a mass spectrometry clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
mssus <- matrix(rep(0,cluesize*2), nrow=cluesize); colnames(mssus) <- c("Sprite","m/z")
for (i in 1:cluesize) {
  mssus[i,1] <- Sprites$Names[susp[i]] ; mssus[i,2] <- PepMass + (hexmass * MS[Sprites$glyco[susp[i]]])
  if(sample(0:99,1)<Error){mssus[i,2] <- PepMass + (hexmass * MS[sample(1:6,1)])}
}
mssus[cluesize,1] <- "Impostor"
mssus2 <- as.data.frame(mssus)


# Obtain an HPLC clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
hpsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
names <- rep(" ", cluesize)
for (i in 1:cluesize) {
  hpsus[i,1] <- susp[i] ; hpsus[i,2] <- Sprites$chem1[susp[i]]; hpsus[i,3] <- Sprites$chem2[susp[i]]
  hpsus[i,4] <- Sprites$chem3[susp[i]]
  if(sample(0:99,1)<Error){hpsus[i,2] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,3] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,4] <- sample(1:5,1)}
  names[i] <- Sprites$Names[susp[i]]; 
}
names[cluesize] <- "X: impostor"
phi2 <- 0.4 + (sample(0:4,1)*0.05)
hp2 <- hplc(hpsus,HPLC,names,phi2)


# Carry out the third murder

victims[3] <- murder(length(Sprites$Names),victims,Impostor)

# Obtain a gel clue. Includes the impostor sample. Defined probability that samples are mixed up.

susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
gelsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
for (i in 1:cluesize) {
  gelsus[i,1] <- susp[i] ; gelsus[i,2] <- Prots[1]; gelsus[i,3] <- Prots[Sprites$prot2[susp[i]]]
  gelsus[i,4] <- Prots[Sprites$prot3[susp[i]]]
  if(sample(0:99,1)<Error){gelsus[i,3] <- Prots[sample(1:5,1)]}
  if(sample(0:99,1)<Error){gelsus[i,4] <- Prots[sample(1:5,1)]}
}
png(filename="output/Gel3.png")
gel3 <- Gel(Markers, gelsus, names=Sprites$Names)
dev.off()  

# Obtain a western clue. Includes the impostor sample. Defined probability that samples are mixed up.

susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
wessus <- matrix(rep(0,cluesize*2), nrow=cluesize)
for (i in 1:cluesize) {
  wessus[i,1] <- Sprites$Names[susp[i]] ; wessus[i,2] <- Glyco[Sprites$glyco[susp[i]]]
  if(sample(0:99,1)<Error){wessus[i,2] <- Glyco[sample(1:6,1)]}
}
if(probe == "K"){probe <- "J"} else {probe <- "K"}
png(filename="output/West3.png")
west <- West(Markers,wessus,probe)
dev.off()  

# Obtain a mass spectrometry clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor 
mssus <- matrix(rep(0,cluesize*2), nrow=cluesize); colnames(mssus) <- c("Sprite","m/z")
for (i in 1:cluesize) {
  mssus[i,1] <- Sprites$Names[susp[i]] ; mssus[i,2] <- PepMass + (hexmass * MS[Sprites$glyco[susp[i]]])
  if(sample(0:99,1)<Error){mssus[i,2] <- PepMass + (hexmass * MS[sample(1:6,1)])}
}
mssus[cluesize,1] <- "Impostor"
mssus3 <- as.data.frame(mssus)


# Obtain an HPLC clue. Includes the impostor sample. Defined probability that samples are mixed up.
susp <- suspects((cluesize-1),length(Sprites$Names),victims); susp[cluesize] <- Impostor  
hpsus <- matrix(rep(0,cluesize*4), nrow=cluesize)
names <- rep(" ", cluesize)
for (i in 1:cluesize) {
  hpsus[i,1] <- susp[i] ; hpsus[i,2] <- Sprites$chem1[susp[i]]; hpsus[i,3] <- Sprites$chem2[susp[i]]
  hpsus[i,4] <- Sprites$chem3[susp[i]]
  if(sample(0:99,1)<Error){hpsus[i,2] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,3] <- sample(1:5,1)}
  if(sample(0:99,1)<Error){hpsus[i,4] <- sample(1:5,1)}
  names[i] <- Sprites$Names[susp[i]]; 
}
names[cluesize] <- "X: impostor"
phi3 <- 0.4 + (sample(0:4,1)*0.05)
hp3 <- hplc(hpsus,HPLC,names,phi3)

# Generate output file

rmarkdown::render("code/Output standard.rmd", output_dir = "output", 
                  output_file = paste0("output/Output moderate ",format(Sys.time(), "%y-%m-%d-%H-%M")))


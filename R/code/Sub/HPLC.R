# Function call for the BIO2090 "Among Us" event
# Nicholas Harmer, 2021, 2022

# Define function to provide an HPLC clue; HPLC is the data frame of data for the compounds.
# phi is the percentage of acetonitrile; dt is the resolution in the x-axis; names is the list of suspect names.
# The column model is from Fasoula, 2017. The first 15 min of the column is illustrated
# column length could be altered if desired or included as an additional variable.
# Requires that the "safe" color set has been set (e.g. from rcartocolor package)

hplc <- function(sus,data,names="",phi=0.5,dt=0.01){
  # Set up core parameters
  sample <- length(data[,1])
  t0 <- sample(9:12,1)/10; tmax <- 15; 
  tR <- rep(0,sample); D <- rep(0,sample); h <- rep(0,sample)
  # Set the retention time, width, and peak height for each sample
  for (i in 1:sample){
    tR[i] <- t0*(1+exp(data$c0[i]+(phi*data$c1[i])+(phi*phi*data$c2[i])))
    D[i] <- data$D0[i]+(data$D1[i]*tR[i]); if (D[i] < 0.01) {D[i] <- 0.01} 
    h[i] <- data$h0[i] + (data$h1[i]*tR[i])
  }
  # Set up a "standard" result matrix for each compound
  AbsR <- matrix(rep(0,sample*(1+(tmax/dt))),nrow=sample)
  for (i in 1:sample){
    for (j in 1:(1+(tmax/dt))){
      AbsR[i,j] <- 20*h[i]*exp(-((((j-1)*dt)-tR[i])^2)/D[i])/(phi^4)
      if (AbsR[i,j] < 0) {AbsR[i,j] <- 0}
    }
  }
  # Set up a data frame (results) for the results for each suspect. x = time (min), y = absorbance (set to zero)
  result <- data.frame(x = rep(seq(0,tmax,by=dt),length(sus[,1])),
                       y = c(0))
  samples <- rep(" ",(1+(tmax/dt))*length(sus[,1]))
  
  # Calculate the y value for each timepoint for each suspect. Loop over suspects (i); timepoints (j); compounds (k)
  for (i in 1:length(sus[,1])){
    for (j in (round((1/dt),0)):(1+(tmax/dt))){
      for (k in 2:length(sus[1,])){
        result[j+((i-1)*(1+(tmax/dt))),2] <- result[j+((i-1)*(1+(tmax/dt))),2]+AbsR[sus[i,k],(j-((i-1)*(round((0.1/dt),0))))]
      }
      result[j+((i-1)*(1+(tmax/dt))),2] <- result[j+((i-1)*(1+(tmax/dt))),2]+(40*(length(sus[,1])-i))
      samples[j+((i-1)*(1+(tmax/dt)))] = names[i]
    }
    # Correct sample to include the first 1/dt values.
    for (j in 1:((round((1/dt),0))-1)) {
      samples[j+((i-1)*(1+(tmax/dt)))] = names[i]
    }
  }
  result$sample <- samples
  
  # Draw an image using ggplot
  hp <- ggplot(data=result, aes(x = x, y = y))
  hp <- hp + geom_line(mapping = aes(color = sample), size = 1) +
    theme_bw(base_size=15) +
    labs(y="Absorbance / mAU", x="time / min") + coord_cartesian(xlim = c(1.6, 15)) +
    theme(axis.title = element_text(size = 20)) +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.y = element_blank()) +
    theme(legend.text = element_text(size = 20)) +
    theme(legend.title = element_text(size = 20)) +
    scale_colour_manual(values=safe)
  output <- hp
}
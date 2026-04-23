# Life_underground_rodents

Data and code related to the study “The many faces of life underground: imperfect convergent evolution of the cranial musculo
skeletal system in subterranean rodents.” 

The scripts require libraries scales, nlme, convevol and viridis, as well as their dependencies to be installed. To run the analyses and
create the plots included in the study (and more), open R (>= 4.5.0) set the cloned repository as working directory then run : 

source(“011_analysis_muscles_fossorials.R”) 

The first lines of the “011_analysis_muscles_fossorials.R” script can be modified to avoid creating the plots, to create plots with species names, to sink results in a text
file or not, and to change the number of simulation in function convSig. Reducing this latter number can drastically reduce
execution time of the script, which can tkae several hours to run if NSIM <- 1000

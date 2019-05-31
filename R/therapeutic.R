#' 
#' Data come from an experimental study (Rosa, Rosa, Sarner, and Barrett, 1998) 
#' that were also used in Howell (2012, p.196).
#' An experiment was conducted to investigate if Therapeutic Touch practitioners 
#' who were blindfolded can effectively identify which of their hands is below the experimenter¡¯s.
#' Twenty-eight practitioners were involved and tested 10 times in the experiment. 
#' Researchers expected an average of 5 correct answers from each practitioner 
#' as it is the number by chance if they do not outperform others.
#' 


## Load data from local source. Please set working directory where data stored.
## setwd(D:/...)
ttestdata<-read.csv("therapeutic.csv")
#ttestdata is a data frame with one column (correct answers) and 28 rows (number of practitioners).

#overall goal of experiment: check for difference between open and closed fjords 
#in taxa/biomass/environmental factors/size distribution.
#TO DO: enter other taxa than copepods, and size of copepods in excel, get CTD data


setwd("~/R/Vitenskapelig_metode/Project")
library(readr)
library(ggmap)
library(ggplot2)

#Import Data
DATA0 <- read.csv("Zooplankton_data.csv", head = T)
View(DATA0)
DATA0[is.na(DATA0)] <- 0 #TEMPORARY

#Making a new dataframe based on the mean of the sub samples
DATA_mean <- aggregate(Copepods ~ Site, data = DATA0, mean)#mean of sub-samples
DATA_temp <- (DATA0[c(1,4,7,10,13,16),4:11])
DATA <- cbind(DATA_mean, DATA_temp)
rm(DATA_temp, DATA_mean)#clean up environment

View(DATA)

#DATA ANALYSIS
plot(DATA0[,-c(1,10,11)])#general plot without coordinates and ID

#closer inspection of correlation between copepods and biomass
ggplot(data= DATA, aes(x = Copepods, y = Biomass)) +
  geom_point() +
  geom_smooth(method = lm)

cor.test(DATA$Copepods, DATA$Biomass)

ggplot(data= DATA0, aes(x = Copepods, y = Biomass)) +
  geom_point() +
  geom_smooth(method = lm)

cor.test(DATA0$Copepods, DATA0$Biomass)
cor.test(DATA[1:3,2],DATA[4:6,2])


plot(DATA0$Copepods ~ DATA0$Site)#Interesting plot

#watercolumn sample function x/m^3_________________________________________________
sample_density = function(x) {
  x <- (x * 200)/18.75# 200ml sample from half of the sample (2x), devided by waterculm volume
  
  return(x)
}

a <- sample_density(DATA$Biomass/1000)# biomass in seawater (g/m^2)

b <- sample_density(DATA$Copepods)# copepods in seawater (n/m^2)

plot(a,b)

#MAP___________________________________________________________________
library(ggmap)
pointLabels <- annotate("text",x=DATA0$lon,y=c(DATA0$lat),size=3,label=as.vector(DATA0$Site))
map <- qmplot(lon, lat, data = DATA0, maptype = "toner-lite", color = I("red")) + pointLabels
map#sampling sites OBS! - M3 and M2 have the wrong coordinates

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
DATA_mean <- aggregate(Copepods ~ Sample, data = DATA0, mean)#mean of sub-samples
DATA_temp <- (DATA0[c(1,4,7,10,13,16),4:12])
DATA <- cbind(DATA_mean, DATA_temp)
rm(DATA_temp, DATA_mean)#clean up environment

View(DATA)

#DATA ANALYSIS____________________________________________________________________
plot(DATA0[,-c(1,10,11)])#general plot without coordinates and ID

#closer inspection of correlation between copepods and biomass
#mean sample plot
ggplot(data= DATA, aes(x = Biomass, y = Copepods)) +
  geom_point() +
  geom_smooth(method = lm)

cor.test(DATA$Copepods, DATA$Biomass)

#all samples plot
ggplot(data= DATA0, aes(x = Biomass, y = Copepods)) +
  geom_point() +
  geom_smooth(method = lm) +
  annotate("text",x=DATA0$Biomass,y=c(DATA0$Copepods),size=3,label=as.vector(DATA0$Sample))


#simple linear regression 
fit = lm(Copepods ~ Biomass, data=DATA0)
fit0 = lm(Copepods ~ 0+Biomass, data=DATA0)
plot(DATA0$Copepods ~DATA0$Biomass)
abline(fit, col="red")
abline(fit0, col="green")


cor.test(DATA0$Copepods, DATA0$Biomass)#seems to be correlation between biomass and copepods counted


plot(DATA0$Copepods ~ DATA0$Sample)#sample summary

#watercolumn sample function x/m^3________________________________________________
sample_density = function(x) {
  x <- (x * 200)/18.75# 200ml sample from half of the sample (2x), devided by waterculm volume
  
  return(x)
}

a <- sample_density(DATA$Biomass/1000)# biomass in seawater (g/m^2)
b <- sample_density(DATA$Copepods)# copepods in seawater (n/m^2)
plot(a,b)

#MAP______________________________________________________________________________
pointLabels <- annotate("text",x=DATA0$lon,y=c(DATA0$lat),size=3,label=as.vector(DATA0$Sample))
map <- qmplot(lon, lat, data = DATA0, maptype = "toner-lite", color = I("red")) + pointLabels
map#sampling sites OBS! - M3 and M2 have the wrong coordinates

#Difference in Biomass between fjords_____________________________________________
plot(DATA0$Biomass ~ DATA0$Fjord)
plot(DATA0$Copepods ~DATA0$Fjord)
ggplot(DATA0, aes(x=Fjord, y=Biomass)) + 
  geom_boxplot()

wilcox.test(DATA0$Biomass~DATA0$Fjord, paired=F, exact=T)

mean_data <- aggregate(cbind(Biomass, Copepods) ~ Fjord, data = DATA0, mean)
t.data.frame(mean_data)
summary(DATA0$Biomass)
#test for difference in two different sites, use students t-test__________________


hist(DATA[1:3,3], freq = FALSE, xlab = "Biomass", main = "Distribution of Biomass", col = "lightgray", xlim = c(1700,2200), ylim = c(0,0.015),)
hist(DATA[4:6,3], add = T, freq = FALSE, col = "darkgray", xlim = c(1600,2300), ylim = c(0,0.015))
curve(dnorm(x, mean=mean(DATA[1:3,3]), sd = sd(DATA0[1:9,4])), add = TRUE, col="Orange", lwd = 2)
curve(dnorm(x, mean=mean(DATA[4:6,3]), sd = sd(DATA0[10:18,4])), add = TRUE, col="black", lwd = 2)

boxplot(DATA[1:3,3], DATA[4:6,3])#guess this gives more information than the data above

#t-test mess 
mean(DATA[1:3,3])#Mistfjord
sd(DATA[1:3,3])#Mistfjord
var(DATA[1:3,3])#Mistfjord

mean(DATA[4:6,3])#Nordfjord
sd(DATA0[4:6,3])#Nordfjord
var(DATA0[4:6,3])#Nordfjord



plot(DATA$Biomass ~ DATA$Fjord)

a <- t.test(Biomass ~ Fjord, data = DATA)#more noise than signal t < 1
str(a)
a

t.test(Biomass ~ Sample == c("M1", "M2", "M4"), data = DATA)#t-test for inside the fjords



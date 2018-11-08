#overall goal of experiment: check for difference between open and closed fjords 
#in taxa/biomass/environmental factors/size distribution.
#TO DO: enter other taxa than copepods, and size of copepods in excel, get CTD data

#Setup ####
setwd("~/R/zooplankton")

#LIBRARY 
library(readr)
library(reshape)
library(ggmap)
library(ggplot2)
library(oce)
library(effects)


#Import Data ####

#ZOOPLANKTON
DATA0 <- read.csv("Zooplankton_data.csv", head = TRUE)
View(DATA0)

#Run this code for scaling of data
MEAN <-function(y){do.call(rbind,as.list(by(y,DATA0$Site,mean)))} #new dataframe based on the mean of the sub samples
DATA_mean <- apply(DATA0[,4:17],2,MEAN)
DATA_temp <- (DATA0[c(1,4,7,10,13,15),c(1,2,3,18,19,20)])
DATA <- cbind(DATA_temp, DATA_mean)
rm(DATA_temp, DATA_mean, MEAN) #Clean up environment
View(DATA)


#watercolumn sample function x/m^3
sample_density_cop = function(x) {
  x <- (x * 200)/18.75# 200ml sample from half of the sample (2x), devided by waterculm volume
  
  return(x)
}
DATA[,7:20] <- sample_density_cop(DATA[,7:20])

#CTD
CTD <- read.csv("CTD_mist - Copy.txt", sep = ";") #Import CTD as txt file
CTD$X.1 <- NULL#remove NA value
CTD$X <- NULL#remove NA value
CTD<-CTD[!(CTD$Depth.u.==0),] #remove CTD entry from 0m depth
colnames(CTD) <- c("Ser", "Meas", "Salinity", "Temp", "Oxygen", "mg.l", "Fluoride", "Density",  "Depth", "Date", "Time") #rename columns
CTD$Seq <- with(CTD, ave(seq_along(Ser), Ser, FUN=seq_along))#seq for Ser
CTD$Ser <- paste("M", CTD$Ser)
CTD$Fjord[429:770] <- "saltfjord";CTD$Fjord[1:428] <- "mistfjord"
View(CTD)

CTD_mean <- aggregate(. ~Ser, data=CTD, mean, na.rm=TRUE)
View(CTD_mean)

#Taxanomic plot ####
library(reshape)
b <- melt(DATA[,-c(3,5:7,1,20)])

ggplot(data = b, aes(x=Site, y= value/100 , fill=variable)) + 
  geom_bar(stat="identity",position="fill", colour="black") +
  labs(y = "Relative abundance")

ggplot(data = b, aes(x="", fill=variable)) + 
  geom_bar(position="fill", colour="black") +
  labs(y = "Relative abundance")+ 
  facet_grid(. ~ Fjord)+
  coord_polar(theta = "y")

#CTD plots ####


#
ggplot(data = CTD, aes(x = Depth, y = Temp, by = Ser, color = Ser)) + 
  geom_line() +
  labs(x = "Depth", color = "Sites") + 
  facet_grid(. ~ Fjord)


ggplot(data = CTD, aes(x = Temp, y = Fluoride, by = Ser, color = Ser))+
  geom_point()+
  facet_grid(. ~ Ser) +
  ylim(0,3)

#MAP ####
pointLabels <- annotate("text",x=DATA0$Lon,y=c(DATA0$Lat),size=3,label=as.vector(DATA0$Site))
map <- qmplot(Lon, Lat, data = DATA0, zoom = 10, maptype = "toner-lite", color = I("red")) + pointLabels
map

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


MEAN <-function(y){do.call(rbind,as.list(by(y,DATA0$Site,mean)))} #new dataframe based on the mean of the sub samples
DATA_mean <- apply(DATA0[,4:17],2,MEAN)
DATA_temp <- (DATA0[c(1,4,7,10,13,15),c(1,2,3,18,19,20)])
DATA <- cbind(DATA_temp, DATA_mean)
rm(DATA_temp, DATA_mean, MEAN) #Clean up environment
View(DATA)
a <- DATA

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
CTD$Fjord[928:1846] <- "saltfjord";CTD$Fjord[1:927] <- "mistfjord"


CTD$Fjord[418:791] <- "saltfjord";CTD$Fjord[1:417] <- "mistfjord"
View(CTD)

CTD_mean <- aggregate(. ~Ser, data=CTD, mean, na.rm=TRUE)
View(CTD_mean)

#FIND CODE FOR DELETING DEPTH AFTER BOTTOM IS REACHED, TEMP SOLUTION: CHANGE DEPTH MANUALY IN TXT FILE


#Taxanomic plot ####
library(reshape)
b <- melt(DATA[,-c(3:6,1,8,20)])

ggplot(data = b, aes(x=Site, y= value/100 , fill=variable)) + 
  geom_bar(stat="identity",position="fill", colour="black") +
  labs(y = "Relative abundance")


#CTD plots ####


ggplot(data = CTD, aes(x = Seq, y = Temp, by = Ser, color = Ser)) + 
  #geom_point() + 
  #geom_line() + 
  geom_smooth() + 
  labs(x = "Time in sec")

ggplot(data = CTD, aes(x = Depth, y = Temp, by = Ser, color = Ser)) + 
  #geom_point() + 
  geom_line() +
  #geom_smooth() + 
  labs(x = "Depth", color = "Sites") + 
  facet_grid(. ~ Fjord)


#DATA ANALYSIS ####
plot(DATA0[,-c(1,5:16,19:20)])#general plot with copepods_all and zooplankton_all and without coordinates, ID


#simple linear regression ####
fit = lm(data=DATA0, Biomass ~ Copepods)
plot(DATA0$Biomass ~ DATA0$Copepods, xlim = c(0,450), ylim=c(1500,2400))
abline(fit, col="green")


cor.test(DATA0$Copepods, DATA0$Biomass)#seems to be correlation between biomass and copepods counted


plot(DATA0$Copepods ~ DATA0$Site)#sample summary

fit = lm(data = DATA0, Zooplankton_all ~ Copepods)
plot(data = DATA0, Zooplankton_all ~Copepods, xlim = c(0,450), ylim=c(0,800))
abline(fit, col="gray")
#COMPARE ZOOPLANKTON_ALL TO BIOMASS WITH THE SAME TEST AS ABOVE


#MAP ####
pointLabels <- annotate("text",x=DATA0$Lon,y=c(DATA0$Lat),size=3,label=as.vector(DATA0$Site))
map <- qmplot(Lon, Lat, data = DATA0, zoom = 10, maptype = "toner-lite", color = I("red")) + pointLabels
map

#Difference in Biomass between fjord ####
plot(DATA$Biomass ~ DATA$Fjord)
plot(DATA$Copepods ~DATA$Fjord)

wilcox.test(DATA$Biomass~DATA$Fjord, paired=F, exact=T)

mean_data <- aggregate(cbind(Biomass, Copepods) ~ Fjord, data = DATA, sum)
t.data.frame(mean_data)

sample_density_bio(mean_data$Biomass) #total Biomass mg/m^3 in the fjords


#ANOVA Test Biomass between fjords ####
#W shapiro test
fit <- aov(Biomass ~ Fjord, data = DATA)
summary.lm(fit)
summary(fit)

boxplot(Biomass~Fjord, data=DATA,ylab="Biomass",col="lightblue")

leveneTest(Biomass~Fjord, data=DATA) # Levene's test, F=6.3527, p=0.002259, i.e. variances significantly non-homogenous
max(by(DATA$Biomass,DATA$Fjord,sd))^2/min(by(DATA$Biomass,DATA$Fjord,sd))^2 # 3.25, a rule of thumb is that if the ratio of the highest to the lowest variance is < 5, variances are sufficiently close to equality
# since Levene's test indicates a significant deviation from homogeneity, we can use
car::Anova(fit,type="III",white.adjust=TRUE) 

#ANOVA test difference for sampling sites ####
fit <- aov(Zooplankton_all ~ Site, data = DATA0)
summary.lm(fit)
summary(fit)

plot(allEffects(fit,confidence.level = 0.95))
boxplot(Zooplankton_all~Site, data=DATA0,ylab="Biomass",col="lightblue")

#errr, 2-way anova zooplankton_all by fjord and biomass
fit <- aov(Zooplankton_all ~ Fjord + Biomass, data = DATA)
summary(fit)

plot(fit)
plot(allEffects(fit,confidence.level = 0.95))


#WRONG TEST - test for difference in two different sites, use students t-test ####
hist(DATA$Biomass[DATA$Fjord == "Mistfjord"], freq = FALSE, xlab = "Biomass", main = "Distribution of Biomass", col = "lightgray", xlim = c(-50,850), ylim = c(0,0.008),)
hist(DATA$Biomass[DATA$Fjord == "Nordfjord"], add = T, freq = FALSE, col = "darkgray", xlim = c(1600,2300), ylim = c(0,0.008))
curve(dnorm(x, mean=mean(DATA[1:3,3]), sd = sd(DATA[1:3,3])), add = TRUE, col="Orange", lwd = 2)
curve(dnorm(x, mean=mean(DATA[4:6,3]), sd = sd(DATA[4:6,3])), add = TRUE, col="black", lwd = 2)

#t-test mess 
mean(DATA[1:3,3])#Mistfjord
mean(DATA$Biomass[DATA$Fjord == "Mistfjord"])
sd(DATA[1:3,3])#Mistfjord
var(DATA[1:3,3])#Mistfjord

mean(DATA[4:6,3])#Nordfjord
sd(DATA0[4:6,3])#Nordfjord
var(DATA0[4:6,3])#Nordfjord

plot(DATA$Biomass ~ DATA$Fjord)

a <- t.test(Biomass ~ Fjord, data = DATA)#more noise than signal t < 1
str(a)

t.test(Biomass ~ Sample == c("M1", "M2", "M4"), data = DATA)#t-test for inside the fjords


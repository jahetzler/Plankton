#overall goal of experiment: check for difference between open and closed fjords 
#in taxa/biomass/environmental factors/size distribution.
#TO DO: enter other taxa than copepods, and size of copepods in excel, get CTD data

#Working dir ####
setwd
("~/R/Vitenskapelig_metode/Project")

library(readr)
library(ggmap)
library(ggplot2)
library(oce)
library(car)
#Import Data ####
DATA0 <- read.csv("Zooplankton_data_rev.csv", head = TRUE)
View(DATA0)
#DATA0[is.na(DATA0)] <- 0 #TEMPORARY

#Making a new dataframe based on the mean of the sub samples
DATA_mean <- aggregate(Copepods ~ Sample, data = DATA0, mean)#mean of sub-samples
DATA_temp <- (DATA0[c(1,4,7,10,13,15),c(4,5,8)])
DATA <- cbind(DATA_mean, DATA_temp)
rm(DATA_temp, DATA_mean)#clean up environment

View(DATA)
#Import CTD ####
a <- read.oce(system.file("extdata", "ctd.cnv", package="oce"))
plot(a)



#Import CTD as txt file
CTD <- read.csv("CTD_mist - Copy.txt", sep = ";")
CTD$X.1 <- NULL#remove NA value
CTD$X <- NULL#remove NA value
CTD<-CTD[!(CTD$Depth.u.==0),] #remove CTD entry from 0m depth

CTD_mean <- aggregate(. ~Ser, data=CTD, mean, na.rm=TRUE)
#
#DATA ANALYSIS ####
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


#simple linear regression ####
fit = lm(Copepods ~ Biomass, data=DATA0)
fit0 = lm(Copepods ~ 0+Biomass, data=DATA0)
plot(DATA0$Copepods ~DATA0$Biomass)
abline(fit, col="red")
abline(fit0, col="green")


cor.test(DATA0$Copepods, DATA0$Biomass)#seems to be correlation between biomass and copepods counted


plot(DATA0$Copepods ~ DATA0$Sample)#sample summary

#watercolumn sample function x/m^3 ####
sample_density_cop = function(x) {
  x <- (x * 200)/18.75# 200ml sample from half of the sample (2x), devided by waterculm volume
  
  return(x)
}
sample_density_bio = function(x){
  x <- (x*2)/18.75
  return(x)
}

Bio_cubic <- sample_density_bio(DATA$Biomass)# biomass in seawater (mg/m^2)
Cop_cubic <- sample_density_cop(DATA$Copepods)# copepods in seawater (n/m^2)
plot(Bio_cubic,Cop_cubic)

#MAP ####
pointLabels <- annotate("text",x=DATA0$Lon,y=c(DATA0$Lat),size=3,label=as.vector(DATA0$Sample))
map <- qmplot(Lon, Lat, data = DATA0, maptype = "toner-lite", color = I("red")) + pointLabels
map#sampling sites OBS! - M3 and M2 have the wrong coordinates

#Difference in Biomass between fjord ####
plot(DATA$Biomass ~ DATA$Fjord)
plot(DATA$Copepods ~DATA$Fjord)
ggplot(DATA, aes(x=Fjord, y=Biomass)) + 
  geom_boxplot()

wilcox.test(DATA$Biomass~DATA$Fjord, paired=F, exact=T)

mean_data <- aggregate(cbind(Biomass, Copepods) ~ Fjord, data = DATA, sum)
t.data.frame(mean_data)

sample_density_bio(mean_data$Biomass) #total Biomass mg/m^3 in the fjords

#test for difference in two different sites, use students t-test ####


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
a

t.test(Biomass ~ Sample == c("M1", "M2", "M4"), data = DATA)#t-test for inside the fjords



#ANOVA Test Biomass between fjords ####
layout(c(1,1))
plot(DATA$Biomass ~ DATA$Fjord)

#do W shapiro test
fit <- aov(Biomass ~ Fjord, data = DATA)
summary.lm(fit)
summary(fit)

boxplot(Biomass~Fjord, data=DATA,ylab="Biomass",col="lightblue")

leveneTest(Biomass~Fjord, data=DATA) # Levene's test, F=6.3527, p=0.002259, i.e. variances significantly non-homogenous
max(by(DATA$Biomass,DATA$Fjord,sd))^2/min(by(DATA$Biomass,DATA$Fjord,sd))^2 # 3.25, a rule of thumb is that if the ratio of the highest to the lowest variance is < 5, variances are sufficiently close to equality
# since Levene's test indicates a significant deviation from homogeneity, we can use
car::Anova(fit,type="III",white.adjust=TRUE) 

#ANOVA test difference for sampling sites ####
fit <- aov(Copepods ~ Sample, data = DATA0)
summary.lm(fit)
summary(fit)

library(effects)
plot(allEffects(fit,confidence.level = 0.95))

boxplot(Copepods~Sample, data=DATA0,ylab="Biomass",col="lightblue")

leveneTest(Biomass~Fjord, data=DATA) # Levene's test, F=6.3527, p=0.002259, i.e. variances significantly non-homogenous
max(by(DATA$Biomass,DATA$Fjord,sd))^2/min(by(DATA$Biomass,DATA$Fjord,sd))^2 # 3.25, a rule of thumb is that if the ratio of the highest to the lowest variance is < 5, variances are sufficiently close to equality
# since Levene's test indicates a significant deviation from homogeneity, we can use
car::Anova(fit,type="III",white.adjust=TRUE) 


#errr, 2-way anova copepods from fjord and sample
fit <- aov(Copepods~Fjord + Sample, data = DATA0)
summary(fit)

plot(fit)
plot(allEffects(fit,confidence.level = 0.95))


#check if there is correlation between biomass and number of copepods


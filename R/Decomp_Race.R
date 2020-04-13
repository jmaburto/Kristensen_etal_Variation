############# Written by JMA
############# 11/04/2020
rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(DemoDecomp)

# Loading data
load("Data/US_HIVmortality.RData")

Data <- data.table(HIV_mortality)
Data <- Data[order(Year,State,Gender,Race,Age_group.code)]
Data[HIV_deaths >0]


# source functions
source('R/Functions_1.R')

# Choose values for epsilon
epsilon <- c(1)
#x <- 0

#example of age-specific decomp
mx1<-HMDL[Year ==2010 & Sex == 'f' & PopName == 'ESP']$mx
mx2<-HMDL[Year ==2015 & Sex == 'f' & PopName == 'ESP']$mx

Efficiency.decomp    <- horiuchi(func = Efficiency.fun,pars1 = mx1,pars2 =mx2 ,N = 50 ,sex = 'f')
plot(Efficiency.decomp)
Inequality.decomp    <- horiuchi(func = Inequality.fun,pars1 = mx1,pars2 =mx2 ,N = 50 ,sex = 'f', epsilon = 1)
plot(Inequality.decomp)


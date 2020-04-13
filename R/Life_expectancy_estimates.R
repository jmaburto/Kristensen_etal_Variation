############# Written by JMA
############# 11/04/2020
rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)

# Loading data
load('Data/Ungrouped_estimates.RData')

# source functions
source('R/Functionsl.R')

# check life expectancy estimates
Single.mx[,Sex1 :=  tolower(as.character(Gender.code))]

Single.mx[Year == 2000 & State.code == 49]

eo.estimates <- Single.mx[,list(eo = LifeExpectancy(Mx,sex = ifelse(Sex1[1]== 'm','m','f'))),
                          by = list(Year,State,Race,Gender,State.code,Race.code,Gender.code)]

#check for NAs
eo.estimates[is.na(eo)]
min(eo.estimates$eo)

save(eo.estimates,file = 'Results/Check_eo.RData')

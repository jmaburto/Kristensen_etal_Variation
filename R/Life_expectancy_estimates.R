############# Written by JMA
############# 11/04/2020
rm(list=ls(all=TRUE))
library(data.table)
library(reshape2)

# Loading data
load('Data/Single_mx.RData')

# source functions
source('R/Functionsl.R')

# check life expectancy estimates

eo.estimates <- Single.mx[,list(eo = LifeExpectancy(mx,sex = ifelse(Gender.code[1]== 'M','m','f'))),
                          by = list(Year,State,Race,Gender.code,State.code,Race.code)]

min(eo.estimates$eo)
max(eo.estimates$eo)

eo.estimates[eo <65]
eo.estimates[eo > 80]

#check for NAs
save(eo.estimates,file = 'Results/Check_eo.RData')

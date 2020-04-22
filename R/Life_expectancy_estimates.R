############# Written by JMA
############# 11/04/2020
rm(list=ls(all=TRUE))
library(data.table)
library(reshape2)

# Loading data
load('Data/Single_mx.RData')

# source functions
source('R/Functionsl.R')

#selection of states
select <- data.table(read.csv('Data/include_states.csv',sep = ';',stringsAsFactors = F))
selection <- unique(select[include == 1]$state)
# check life expectancy estimates

Single.mx <- Single.mx[State %in% selection]

Results.e0.edagger <- Single.mx[,list(eo = LifeExpectancy(mx,sex = ifelse(Gender.code[1]== 'M','m','f')),
                                e.dagger = edagger.frommx(mx,sex = ifelse(Gender.code[1]== 'M','m','f'))),
                          by = list(Year,State,Race,Gender.code,State.code,Race.code)]

min(eo.estimates$eo)
max(eo.estimates$eo)

eo.estimates[eo <65]
eo.estimates[eo > 80]

#check for NAs
save(Results.e0.edagger,file = 'Results/eo_edagger.RData')

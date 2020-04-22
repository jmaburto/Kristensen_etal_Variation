############# Written by JMA
############# 11/04/2020
rm(list=ls(all=TRUE))
library(data.table)
library(reshape2)
library(DemoDecomp)
library(ggplot2)

# Loading data
load('Data/Single_mx.RData')

# source functions
source('R/Functionsl.R')

#selection of states
select <- data.table(read.csv('Data/include_states.csv',sep = ';',stringsAsFactors = F))
selection <- unique(select[include == 1]$state)
# check life expectancy estimates

Single.mx <- Single.mx[State %in% selection]
Counts.data <- Counts.data[State %in% selection]
Counts.data[,HIV.prop := HIV/Deaths]


Single.mx   <- Single.mx[order(Year,State,Race,Age)]
Counts.data <- Counts.data[order(Year,State,Race,Age)]

Single.mx$HIV.prop <- Counts.data$HIV.prop
Single.mx$rest.prop <- 1-Single.mx$HIV.prop 

#DT <- Single.mx[Race.code == '2054-5' & State == 'Alabama' & Gender.code == 'M']

time.decomp.results.ed <-Single.mx[,decomp.function.time1(DT = .SD,sex = tolower(Gender.code[1])), 
                                   by = list(State,State.code,Gender.code,Race,Race.code)]
time.decomp.results.ed$Indicator <- 'e.dagger'

time.decomp.results.eo <-Single.mx[,decomp.function.time2(DT = .SD,sex = tolower(Gender.code[1])), 
                                   by = list(State,State.code,Gender.code,Race,Race.code)]

time.decomp.results.eo$Indicator <- 'e0'

decomp.results.time <- rbind(time.decomp.results.ed,time.decomp.results.eo)


############### for now take out NA

decomp.results.time <- decomp.results.time[!(is.na(rest))]
decomp.results.time[(is.infinite(hiv))] <- 0


results.time <- melt(decomp.results.time,id.vars = c('State','State.code','Gender.code','Race', 'Race.code','period','Indicator','age'),
                     variable.name = 'Cause',value.name = 'Contribution')

results.time[,Age.5:= cut(age,c(seq(0,100,5),Inf),include.lowest = T,right = F,labels = seq(0,100,5))]

results.time.5 <- results.time[, list(Contribution = sum(Contribution)), 
                               by = list(State,State.code,Gender.code,Race,Race.code,period,Indicator,Cause,Age.5)] 

results.time.5$Age.5 <- as.numeric(as.character(results.time.5$Age.5))
results.time.5$Cause <- as.character(results.time.5$Cause)

save(decomp.results.time,results.time.5, file = 'Results/Decomp_results_time.RData')


ggplot() +
  ggtitle( 'Age-cause-contribution to difference in e.dagger 1980-1990' , subtitle = 'Males')+
  facet_wrap(~State)+
  geom_bar(data = results.time.5[period == '1980-1990' & Gender.code == 'M' & Indicator == 'e.dagger' & Age.5 < 100 & Race.code == '2054-5'],
           aes(x = Age.5, y = Contribution, fill = Cause,group = Cause), stat = "identity",position = "stack")

ggplot() +
  ggtitle( 'Age-cause-contribution to difference in e0 1980-1990' , subtitle = 'Males')+
  facet_wrap(~State)+
  geom_bar(data = results.time.5[period == '1980-1990' & Gender.code == 'M' & Indicator == 'e0' & Age.5 < 100 & Race.code == '2054-5'],
           aes(x = Age.5, y = Contribution, fill = Cause,group = Cause), stat = "identity",position = "stack")

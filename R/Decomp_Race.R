############# Written by JMA
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

#DT <- Single.mx[Year == 1990 & State == 'Alabama' & Gender.code == 'M']

decomp.results.ed <-Single.mx[,decomp.function(DT = .SD,sex = tolower(Gender.code[1])), by = list(Year,State,Gender.code,State.code)]
decomp.results.ed$Indicator <- 'e.dagger'

decomp.results.eo <-Single.mx[,decomp.function2(DT = .SD,sex = tolower(Gender.code[1])), by = list(Year,State,Gender.code,State.code)]
decomp.results.eo$Indicator <- 'e0'

decomp.results <- rbind(decomp.results.ed,decomp.results.eo)

decomp.results[,Age:= 0:110, by = list(Year,State,Gender.code,Indicator)]

names(decomp.results)[5:6] <- c('Rest', 'HIV')

decomp.results <- melt(decomp.results,id.vars = c('Year','State','State.code','Gender.code','Indicator','Age'),variable.name = 'Cause',value.name = 'Contribution')

decomp.results[,Age.5:= cut(Age,c(seq(0,100,5),Inf),include.lowest = T,right = F,labels = seq(0,100,5))]

decomp.results.5 <- decomp.results[, list(Contribution = sum(Contribution)), by = list(Year,State,State.code,Gender.code,Indicator,Cause,Age.5)] 
decomp.results.5$Age.5 <- as.numeric(as.character(decomp.results.5$Age.5))

  ggplot() +
    ggtitle( 'Age-cause-contribution to difference in e.dagger (B - W)' , subtitle = 'Males')+
    facet_wrap(~State)+
    geom_bar(data = decomp.results.5[Year == '1990' & Gender.code == 'M' & Indicator == 'e.dagger' & Age.5 < 100],
             aes(x = Age.5, y = Contribution, fill = Cause,group = Cause), stat = "identity",position = "stack")
  
  ggplot() +
    ggtitle( 'Age-cause-contribution to difference in life expectancy (W - B)' , subtitle = 'Males')+
    facet_wrap(~State)+
    geom_bar(data = decomp.results.5[Year == '1990' & Gender.code == 'M' & Indicator == 'e0' & Age.5 < 100],
             aes(x = Age.5, y = -Contribution, fill = Cause,group = Cause), stat = "identity",position = "stack")


save(decomp.results,decomp.results.5, file = 'Results/Decomp_results.RData')


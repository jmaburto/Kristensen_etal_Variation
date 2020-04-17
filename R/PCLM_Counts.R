### written by JMA 11/04/2020
### program to have single age estimates
library(data.table)
library(reshape2)
library(ungroup)

rm(list=ls(all=TRUE))

load("Data/US_HIVmortality_update_16042020.RData")

source('R/Functionsl.R')

Data   <- data.table(HIV_mortality)
Data[, State := as.character(State)]
Data[, Race := as.character(Race)]
Data[, Gender := as.character(Gender)]
Data[, Race.code := as.character(Race.code)]
Data[, Gender.code := as.character(Gender.code)]
Data[, Age_group := as.character(Age_group)]
str(Data)

Data   <- Data[order(Year,State,Gender,Race,Age_group.code)]
Data   <- Data[,-c(4,9,10,11,15)]

Age.n  <- c(0,1,5,10,15,20,25,35,45,55,65,75,85)
n.last <- 26

#add note to all 
Data[,Age.lower:= Age.n,by = list(Year,State,Race,Gender.code,State.code,Race.code)]
Data[,Note.HIV:= ifelse(length(Note_HIV_update[Note_HIV_update != ""])>0,'Out','In'),by = list(Year,State,Race,Gender.code,State.code,Race.code)]
Data[,Note.total:= ifelse(length(Note_total_update[Note_total_update != ""])>0,'Out','In'),by = list(Year,State,Race,Gender.code,State.code,Race.code)]

#Deal with Pop, total and hiv separately to have consistency
Total.deaths <- Data
HIV.deaths   <- Data

#Handle HIV.deaths'
vec1         <- which(HIV.deaths$Note_HIV_update != "")-1
HIV.deaths   <- HIV.deaths[-vec1]
HIV.deaths[which(HIV.deaths$Note_HIV_update != "")]$Age.lower <- HIV.deaths[which(HIV.deaths$Note_HIV_update != "")]$Age.lower - 5
HIV.counts <- HIV.deaths[,my.ungroup.fun1(x = Age.lower,y = HIV_deaths_update,nlast = n.last),by = list(Year,State,Race,Gender.code,State.code,Race.code)]

#Handle total deaths
vec2         <- which(Total.deaths$Note_total_update != "")-1
Total.deaths <- Total.deaths[-vec2]
Total.deaths[which(Total.deaths$Note_total_update != "")]$Age.lower <- Total.deaths[which(Total.deaths$Note_total_update != "")]$Age.lower - 5
Total.counts <- Total.deaths[,my.ungroup.fun1(x = Age.lower,y = total_deaths_update,nlast = n.last),by = list(Year,State,Race,Gender.code,State.code,Race.code)]

#regroup in 5-year-age groups

#Population Update
Population.counts <- Data[,my.ungroup.fun1(x = Age.n,y = Population_update,nlast = n.last),by = list(Year,State,Race,Gender.code,State.code,Race.code)]

#save(HIV.counts,Total.counts,Population.counts, file = 'Data/Consistent_counts1.RData')
#load('Data/Consistent_counts.RData')

#put them all together now
Counts.data <- merge(Population.counts,Total.counts,by = c('Year','State','Race','Gender.code','State.code','Race.code','V1'))
Counts.data <- merge(Counts.data,HIV.counts,by = c('Year','State','Race','Gender.code','State.code','Race.code','V1'))
names(Counts.data)[7:10] <- c('Age','Population','Deaths','HIV')

#Do some checks
Counts.data[is.na(Age)]
Counts.data[is.na(Population)]
Counts.data[is.na(Deaths)]
Counts.data[is.na(HIV)]$HIV <- 0

x <- Counts.data[,sum(HIV,na.rm = T), by =  list(Year,State,Race,Gender.code,State.code,Race.code)]
x <- x[order(Year,State,Race,Gender.code)]
y <- Data[,sum(HIV_deaths_update,na.rm = T), by =  list(Year,State,Race,Gender.code,State.code,Race.code)]
y <- y[order(Year,State,Race,Gender.code)]
z <- x$V1 - y$V1


save(Counts.data,file = 'Data/Consistent_counts.RData')



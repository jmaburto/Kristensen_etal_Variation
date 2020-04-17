### written by JMA 11/04/2020
### program to have single age estimates
library(data.table)
library(reshape2)
library(ungroup)

rm(list=ls(all=TRUE))

load("Data/Consistent_counts.RData")

my.ungroup.fun2 <- function(x = Age.n,y = Test.data$total_deaths,nlast = n.last,offset = Test.data$Population){
  m       <- pclm(x = x,y = y,nlast = nlast,offset = offset, control = list(lambda = NA, opt.method = "AIC"))
  m.table <- list(0:110,m$fitted,m$ci$upper,m$ci$lower)
  return(m.table)
}

Age.n  <- c(0,1,5,10,15,20,25,35,45,55,65,75,85)

Counts.data <- Counts.data[order(Year,State,Race,Gender.code,State.code,Age)]

Counts.data[,Age.5:= cut(Age,c(Age.n,Inf),include.lowest = T,right = F,labels = Age.n)]
Counts.data$Age.5 <- as.numeric(as.character(Counts.data$Age.5)) 

Counts.data.5 <- Counts.data[,list(Population = sum(Population),
                                  Deaths = sum(Deaths),
                                  HIV = sum(Deaths)), by = list(Year,State,Race,Gender.code,State.code,Race.code,Age.5)]

#check <- Counts.data.5[Year == 1980 & State == 'Alabama']


Single.mx <- Counts.data.5[,my.ungroup.fun2(x = Age.5,y = Deaths,nlast = 26,offset = Population), 
                   by = list(Year,State,Race,Gender.code,State.code,Race.code)]


names(Single.mx)[7:10] <- c('Age','mx','mx.upper','mx.lower')

save(Counts.data,Single.mx, file = 'Data/Single_mx.RData')


### written by JMA 11/04/2020
### program to have single age estimates
library(data.table)
library(reshape2)
library(ungroup)

rm(list=ls(all=TRUE))

load("Data/US_HIVmortality.RData")

Data   <- data.table(HIV_mortality)
Data[, State := as.character(State)]
Data[, Race := as.character(Race)]
Data[, Gender := as.character(Gender)]
Data[, Race.code := as.character(Race.code)]
Data[, Gender.code := as.character(Gender.code)]
Data[, Age_group := as.character(Age_group)]
str(Data)

# Exclude Alaska and NAs for now
# Check Utah 2000
Data   <- Data[! (State %in% c('Hawaii','Idaho','Maine','Montana','New Hampshire','North Dakota','South Dakota','Vermont','Wyoming','Utah') )]


Data   <- Data[order(Year,State,Gender,Race,Age_group.code)]
Age.n  <- c(0,1,5,10,15,20,25,35,45,55,65,75,85)
n.last <- 26

## Try one first
# Test.data <- Data[Year == 1990 & State == 'Louisiana' & Gender == 'Male' & Race.code == '2106-3']
# m.total         <- pclm(x = Age.n,y = Test.data$total_deaths,nlast = n.last,offset = NULL,control = list(lambda = NA, opt.method = "AIC"))
# m.rest          <- pclm(x = Age.n,y = Test.data$total_deaths-Test.data$HIV_deaths,nlast = n.last,offset = NULL,control = list(lambda = NA, opt.method = "AIC"))
# plot(m.total)
# plot(m.rest)
# m.total$fitted
# m$ci$upper
# m$ci$lower
# length(m$fitted)

my.ungroup.fun <- function(x = Age.n,y = Test.data$total_deaths,nlast = n.last,offset = Test.data$Population){
  m       <- pclm(x = x,y = y,nlast = nlast,offset = offset, control = list(lambda = NA, opt.method = "AIC"))
  m.table <- list(0:110,m$fitted,m$ci$upper,m$ci$lower)
  return(m.table)
}

# check for NA
Data[is.na(total_deaths)]
Data[is.na(Population)]
Data[is.na(HIV_deaths)]

Single.mx          <-  Data[, my.ungroup.fun(x = Age.n,y = total_deaths,nlast = n.last,offset = Population), 
                            by = list(Year,State,Race,Gender,State.code,Race.code,Gender.code)]
names(Single.mx)[8:11] <- c('Age','Mx','Mx.upper','Mx.lower')

Ungroup.data.total <-  Data[, my.ungroup.fun(x = Age.n,y = total_deaths,nlast = n.last,offset = NULL), 
                      by = list(Year,State,Race,Gender,State.code,Race.code,Gender.code)]
names(Ungroup.data.total)[8:11] <- c('Age','Dx','Dx.upper','Dx.lower')

Ungroup.rest       <-  Data[, my.ungroup.fun(x = Age.n,y = total_deaths-HIV_deaths,nlast = n.last,offset = NULL), 
                            by = list(Year,State,Race,Gender,State.code,Race.code,Gender.code)]
names(Ungroup.rest)[8:11] <- c('Age','Dx.rest','Dx.rest.upper','Dx.rest.lower')

save(Single.mx,Ungroup.data.total,Ungroup.rest, file = 'Data/Ungrouped_estimates.RData')


#We should compare with USA.HMD values


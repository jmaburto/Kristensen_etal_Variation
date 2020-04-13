# Functions needed for 'Life expectancy and equality: a long run relationship'
# a nify del function:

#mx1<-HMDL[Year ==1921 & Sex == 'f' & PopName == 'AUS']$mx
#mx2<-HMDL[Year ==1930 & Sex == 'f' & PopName == 'AUS']$mx
# DT <- HMDL[Sex == 'f' & PopName == 'AUS']
# x <- 1921

Efficiency.fun <- compiler::cmpfun(function(mx,sex = "f"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]                   
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- lx[i.openage ] * ax[i.openage ]
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  ex[1]
})


Inequality.fun <- compiler::cmpfun(function(mx,epsilon = 0,age = 0:110, sex = "f"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]                   
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  
  H_index <- ifelse(epsilon != 1,(sum(dx*(age+ax)^(1-epsilon),na.rm = T))^(1/(1-epsilon)),
                    prod((age+ax)^(dx)))
  
  efficiency <- sum(dx*(age+ax),na.rm = T)
  
  inequality <- 1 - H_index/efficiency
  
  inequality
})



inner.decomp.function <- function(DT = .DF,epsilon = .5,age = 0:110,sex = 'f'){
  years <- unique(DT$Year)
  
  decomp.list <- do.call(rbind,lapply(years[-length(years)], function(x,DT,epsilon = epsilon, age = age,sex = sex){
    decomp.year <- H_decomp(mx1 = DT[DT$Year == x]$mx,
                            mx2 = DT[DT$Year == x + 1]$mx,
                            epsilon = epsilon,
                            age = age,
                            sex = sex)
    
    decomp.year$year.1 <- x
    
    data.table(decomp.year)
    
  }, DT = DT,epsilon = epsilon,age = age,sex = sex))
}

H_decomp <- compiler::cmpfun(function(mx1, mx2, epsilon = .5,age = 0:110, sex = "f"){
  H1 <- H_index(mx = mx1 ,epsilon = epsilon ,age = age, sex = sex)
  H2 <- H_index(mx = mx2 ,epsilon = epsilon ,age = age, sex = sex)
  
  Dif.original <- H2$H_index - H1$H_index
  
  Delta.mu <- .5*(f.function(H2$efficiency,H2$inequality)-f.function(H1$efficiency,H2$inequality)+
                    f.function(H2$efficiency,H1$inequality)-f.function(H1$efficiency,H1$inequality))
  
  Delta.I <- .5*(f.function(H2$efficiency,H2$inequality)-f.function(H2$efficiency,H1$inequality)+
                    f.function(H1$efficiency,H2$inequality)-f.function(H1$efficiency,H1$inequality))
  
  C.mu <- Delta.mu/Dif.original
  
  C.I <- Delta.I/Dif.original
  
  return(data.table(cbind(H2 = H2$H_index,H1 = H1$H_index ,mu2 = H2$efficiency,mu1= H1$efficiency,
                          I2 = H2$inequality,I1=H1$inequality,
                          Delta.mu,Delta.I,C.mu,C.I)))
})

f.function <- function(mu,I){
  mu*(1-I)
}

H_index <- compiler::cmpfun(function(mx,epsilon = 0,age = 0:110, sex = "f"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]                   
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  
  #I need to generalize this to each age
  
  H_index <- ifelse(epsilon != 1,(sum(dx*(age+ax)^(1-epsilon),na.rm = T))^(1/(1-epsilon)),
         prod((age+ax)^(dx)))
  
  efficiency <- sum(dx*(age+ax),na.rm = T)
  
  inequality <- 1 - H_index/efficiency
  
  return(data.table(cbind(H_index,efficiency,inequality)))
})

#H_index(mx = mx,epsilon = .5,age = 0:110,sex = 'f')

# Some useful fucntions: for ax and life table
AKm02a0        <- function(m0, sex = "m"){
  sex <- rep(sex, length(m0))
  ifelse(sex == "m", 
         ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
                ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
         # f
         ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
                ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}

LifeExpectancy <- compiler::cmpfun(function(mx,sex = "f"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]                   
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- lx[i.openage ] * ax[i.openage ]
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  ex[1]
})

LifeTable      <- function(mx,sex = "f"){
  mx <- as.matrix(mx)
  i.openage <- nrow(mx)
  ax        <- mx * 0 + .5
  ax[1, ]   <- AKm02a0(m0 = mx[1, ], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)        
  qx[i.openage, ]       <- ifelse(is.na(qx[i.openage, ]), NA, 1)
  ax[i.openage, ]       <- 1 / mx[i.openage, ]                   
  px 				      <- 1 - qx 																				
  px[is.nan(px)]  <- 0 
  lx 			        <- apply(px, 2, function(px., RADIX, OPENAGE){ 		
    if (all(is.na(px.))) {
      px.
    } else {
      c(RADIX, RADIX * cumprod(px.[1:OPENAGE]))
    }
  }, RADIX = 1, OPENAGE = i.openage - 1
  )
  rownames(lx)    <- 0:(i.openage - 1) 
  dx 				      <- lx * qx 																				
  Lx 				      <- lx - (1 - ax) * dx 														
  Lx[i.openage, ]	<- lx[i.openage, ] * ax[i.openage, ]
  Tx 				      <- apply(Lx, 2, function(Lx., i.openage, OPENAGE){
    c(rev(cumsum(rev(Lx.[1:OPENAGE]))),0) + Lx.[i.openage]	
  }, OPENAGE = i.openage - 1, i.openage = i.openage
  )
  rownames(Tx)    <- rownames(lx)
  ex 				      <- Tx / lx 	                              
  list(e0=ex[1,],ex=ex,lx=lx,mx=mx)
}


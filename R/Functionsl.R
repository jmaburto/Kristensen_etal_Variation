
decomp.function <- function(DT = .SD,sex = sex[1]){
  #sex <- 'm'
  white.mx <- c(DT[Race.code == '2106-3']$mx*DT[Race.code == '2106-3']$rest.prop,
                DT[Race.code == '2106-3']$mx*DT[Race.code == '2106-3']$HIV.prop)
  
  black.mx <- c(DT[Race.code == '2054-5']$mx*DT[Race.code == '2054-5']$rest.prop,
                DT[Race.code == '2054-5']$mx*DT[Race.code == '2054-5']$HIV.prop)
  
  decomp <- horiuchi(edfrommxc, pars1 = white.mx, pars2 = black.mx,N = 50, sex == sex)
  
  dim(decomp) <- c(111,2)
  
  data.table(decomp)
  
}

decomp.function2 <- function(DT = .SD,sex = sex[1]){
  #sex <- 'm'
  white.mx <- c(DT[Race.code == '2106-3']$mx*DT[Race.code == '2106-3']$rest.prop,
                DT[Race.code == '2106-3']$mx*DT[Race.code == '2106-3']$HIV.prop)
  
  black.mx <- c(DT[Race.code == '2054-5']$mx*DT[Race.code == '2054-5']$rest.prop,
                DT[Race.code == '2054-5']$mx*DT[Race.code == '2054-5']$HIV.prop)
  
  decomp <- horiuchi(e0frommxc, pars1 = white.mx, pars2 = black.mx,N = 50, sex == sex)
  
  dim(decomp) <- c(111,2)
  
  data.table(decomp)
  
}


e0frommxc <- function(mxcvec,sex){
  dim(mxcvec) <- c(111,length(mxcvec)/111)
  mx          <- rowSums(mxcvec)
  LifeExpectancy(mx,sex)
}

edfrommxc <- function(mxcvec,sex){
  dim(mxcvec) <- c(111,length(mxcvec)/111)
  mx          <- rowSums(mxcvec)
  edagger.frommx(mx,sex)
}



edagger.frommx <- function(mx,sex){
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
  l <- length(ex)
  v <- (sum(dx[-l]* (ex[-l] + ax[-l]*(ex[-1]-ex[-l]) )) + ex[l])
  return(v)
}



my.ungroup.fun1 <- function(x,y,nlast = n.last,offset = NULL){
  m       <- pclm(x = x,y = y,nlast = nlast,offset = offset, control = list(lambda = NA, opt.method = "AIC"))
  m.table <- list(0:110,m$fitted)
  return(m.table)
}

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


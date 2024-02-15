
# ---------------------------------------------------------------------------- #
# Title:   Population total variation - Functions
# Author:  Gomez-Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   1. Functions

# ---------------------------------------------------------------------------- #
#     1. Functions
# ---------------------------------------------------------------------------- #


# ------Life table from mortality rates------#
# Authors: Carlo-Giovanni Camarda and Ugofilippo Basellini
lifetable.mx <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}


# ------Gompertz density function------#
GompFX <- function(pars,ages){
  ## gompertz parameters
  a <- pars[1]
  b <- pars[2]
  ## density
  t1 <- (a/b)*(exp(b*ages)-1)
  fx <- a*exp(b*ages - t1)
  fx <- fx/sum(fx)  ## proper density
  return(fx)
}


# ------Gompertz force of mortality------#
GompMU <- function(pars,ages){
  ## gompertz parameters
  a <- pars[1]
  b <- pars[2]
  ## force of mortality
  mu <- a*exp(b*ages)
}


# ------Standard deivation of the ages at death------#
sdv_func <- function(LT){
  nages <- dim(LT)[1]
  ex.Age <- Age+LT$ex
  ax.Age <- Age+LT$ax
  
  V <- rev(cumsum(rev(LT$dx*(ax.Age-ex.Age)^2)))/LT$lx
  S <- sqrt(V)
  
  return(as.data.frame(S[1]))
}


# ------Total variation------#
tvd <- function(f1,f2){
  f.diff.abs <- abs(f1-f2)
  out <- (1/2)*sum(f.diff.abs)
  return(out)
}

# ------Standard deivation of the ages at death when age is not included in the LT------#
sdv_func_age <- function(LT){
  Age <- c(0:110)
  nages <- dim(LT)[1]
  exAge <- Age+LT$ex
  axAge <- Age+LT$ax
  
  V <- rev(cumsum(rev(LT$dx*(axAge-exAge)^2)))/LT$lx
  S <- sqrt(V)
  
  return(S)
}

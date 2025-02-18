
# ---------------------------------------------------------------------------- #
# Paper:   Reassessing socioeconomic inequalities in mortality via distributional similarities
# Title:   Functions
# ---------------------------------------------------------------------------- #
# Content:
#   1. Functions

# ---------------------------------------------------------------------------- #
#     1. Functions
# ---------------------------------------------------------------------------- #

# ------Life table from mortality rates------#
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

# ------Non-overlap between two distributions------#
noi_func <- function(a,b){
  FX <- cbind(a,b)
  fmin <- apply(FX,1,min)
  fmax <- apply(FX,1,max)
  
  out <- 1-(sum(fmin)/sum(fmax))
  return(out)
}

# ------Pairwise non-overlap index (Shi et al.)------# 
noi_pair_func_dx <- function(dx_w, n){
  w <- matrix(as.matrix(dx_w[(length(dx_w)/2+1):length(dx_w)]), ncol = n)
  w_dpi <- colSums(w)/sum(w)
  dx <- matrix(as.matrix(dx_w[1:(length(dx_w)/2)]), ncol = n)
  dx.d <- dx/colSums(dx)
  comb_dx <- combn(ncol(dx),2)
  
  Sij <- rep(NA,ncol(comb_dx))
  weights <- rep(NA,ncol(comb_dx))
  i <- 1
  for (i in 1:ncol(comb_dx)){
    Sij[i] <- noi_func(dx[,comb_dx[1,i]],dx[,comb_dx[2,i]])
    weights[i] <- w_dpi[comb_dx[1,i]]*w_dpi[comb_dx[2,i]]
  }
  
  noi_t <- sum(Sij*weights)/sum(weights)
  noi_t
}

# ------Pairwise outsurvival probability------#
ov_pair_func_dx <- function(dx, w, lx, ex, n){
  dx <- matrix(dx, ncol = n)
  w <- matrix(w, ncol = n)
  lx <- matrix(lx, ncol = n)
  ex <- matrix(ex, ncol = n)
  w_dpi <- colSums(w)/sum(w)
  # dx.d <- dx/colSums(dx)
  
  # comb_dx <- combn(ncol(dx),2)
  comb_dx <- t(expand.grid(rep(list(1:n), 2)))
  
  ex_0 <- ex[1,]
  
  m <- nrow(dx)
  OVij <- rep(NA,ncol(comb_dx))
  weights <- rep(NA,ncol(comb_dx))
  i <- 1
  for (i in 1:ncol(comb_dx)){
    if (ex_0[comb_dx[2,i]] > ex_0[comb_dx[1,i]]){
      OVij[i] <- sum(dx[1:(m - 1),comb_dx[1,i]] * lx[2:m,comb_dx[2,i]]) + sum(dx[,comb_dx[1,i]] * dx[,comb_dx[2,i]]) / 2
      weights[i] <- w_dpi[comb_dx[1,i]]*w_dpi[comb_dx[2,i]]      
    } else {next}
  }
  
  OV_t <- sum(OVij*weights, na.rm = TRUE)/sum(weights, na.rm = TRUE)
  OV_t
}

# ------Total variation (for working paper)------#
tvd <- function(f1,f2){
  f.diff.abs <- abs(f1-f2)
  out <- (1/2)*sum(f.diff.abs)
  return(out)
}

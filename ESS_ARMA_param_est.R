# estimate the autocorrelation and ESS
library(robustbase)
library(robcor)
library(tidyverse)
library(tsqn)

# yule Walker -------------------------------------------------------------

arp <- function(x,p){
  cov.para=acf(x,type="covariance",plot=FALSE,lag.max=p) 
  M = toeplitz(cov.para$acf[1:p])
  R = cbind(rbind(cov.para$acf[2:(p+1)],M),c(1,(rep(0,p)))) 
  C = cov.para$acf[1:(p+1)]
  para.YW=solve(R) %*%C 
  phi <-  para.YW[1:p]
  return(phi)
}



# classic with missing data -----------------------------------------------------------------

classic <- function(x){
  n = length(x)
  u <- x[1:(n-1)]
  v <- x[2:n] 
  # Empirical estimator
  # rho0.est <-cov(u,v)/var(x)
  rho0.est <- cor(u,v, use = "complete.obs") # with missing data
  
  return(rho0.est)
}


# pooled autocorrelation --------------------------------------------------

pooled_autocorrelation <- function(x, y){
  m1 = mean(x, na.rm = TRUE)
  m2 = mean(y, na.rm = TRUE)
  x_1 = x[-1]
  x0 = x[-(length(x))]
  y_1 = y[-1]
  y0 = y[-(length(y))]
  nume = sum((x0 - m1)*(x_1 - m1)) + sum((y0 - m2)*(y_1 - m2))
  deno = sum((x - m1)**2) + sum((y - m2)**2)
  r = nume/deno
  return(r)
}


# Ma&Genton autocorrelation ---------------------------------------------------------------
# pooled
ma.genton <- function(x.var){
  x <- x.var
  n = length(x)
  u <- x[1:(n-1)]
  v <- x[2:n] 
  # with missing data 
  Sum.uv <- u+v
  Minus.uv <- u-v
  rg.na <- which(is.na(Sum.uv))
  if (length(rg.na)>0){
    Sum.uv <- Sum.uv[-rg.na]
    Minus.uv <- Minus.uv[-rg.na]
  }
  rhoMaGG.est <- ((robustbase::Qn(Sum.uv))^2-(robustbase::Qn(Minus.uv))^2)/((robustbase::Qn(Sum.uv))^2+(robustbase::Qn(Minus.uv))^2)
  
  # Using the robust estimator of Ma and Genton withour missing data
  # rhoMaGG.est <- corQn(u,v)
  return(rhoMaGG.est)
}

# 2 sides
ar.ma.genton.2side <- function(signal){
  n = length(signal)
  n.2 = floor((length(signal))/2)
  
  x <- signal[1:(n/2)]
  y <- signal[(n/2+1):n]
  ux <- x[1:(n.2-1)]
  vx <- x[2:n.2] 
  uy <- y[1:(n.2-1)]
  vy <- y[2:n.2] 
  
  Sum.uvx <- ux+vx
  Minus.uvx <- ux-vx
  rg.nax <- which(is.na(Sum.uvx))
  if (length(rg.nax)>0){
    Sum.uvx <- Sum.uvx[-rg.nax]
    Minus.uvx <- Minus.uvx[-rg.nax]
  }
  
  Sum.uvy <- uy+vy
  Minus.uvy <- uy-vy
  rg.nay <- which(is.na(Sum.uvy))
  if (length(rg.nay)>0){
    Sum.uvy <- Sum.uvy[-rg.nay]
    Minus.uvy <- Minus.uvy[-rg.nay]
  }
  
  # Using the robust estimator of Ma and Genton
  # rhoMaGG.est <-  ((robustbase::Qn(ux+vx))^2-(robustbase::Qn(ux-vx))^2+(robustbase::Qn(uy+vy))^2-(robustbase::Qn(uy-vy))^2)/((robustbase::Qn(ux+vx))^2+(robustbase::Qn(ux-vx))^2+(robustbase::Qn(uy+vy))^2+(robustbase::Qn(uy-vy))^2)
  
  rhoMaGG.est <- ((robustbase::Qn(Sum.uvx))^2-(robustbase::Qn(Minus.uvx))^2+(robustbase::Qn(Sum.uvy))^2-(robustbase::Qn(Minus.uvy))^2)/((robustbase::Qn(Sum.uvx))^2+(robustbase::Qn(Minus.uvx))^2+(robustbase::Qn(Sum.uvy))^2+(robustbase::Qn(Minus.uvy))^2)
  return(rhoMaGG.est)
}

# Souhil autocorrelation --------------------------------------------------
souhil2 <- function(x.var){
  x <- diff(x.var)
  u <- x[1:(length(x)-1)]
  v <- x[2:(length(x))]
  # without missing data
  # MA.G <-  ((robustbase::Qn(u+v))^2-(robustbase::Qn(u-v))^2)/((robustbase::Qn(u+v))^2+(robustbase::Qn(u-v))^2)
  # with missing data
  Sum.uv <- u+v
  Minus.uv <- u-v
  rg.na <- which(is.na(Sum.uv))
  if (length(rg.na)>0){
    Sum.uv <- Sum.uv[-rg.na]
    Minus.uv <- Minus.uv[-rg.na]
  }
  MA.G <- ((robustbase::Qn(Sum.uv))^2-(robustbase::Qn(Minus.uv))^2)/((robustbase::Qn(Sum.uv))^2+(robustbase::Qn(Minus.uv))^2)
  rho.est2 <-1+2*MA.G #New (phd Souhil page 127)
  return(rho.est2)
}


# ESS ---------------------------------------------------------------------

ESS <- function(n, rho){
  ne = n*(1-rho)/(1+rho)
  return(ne)
}

ESS.ar <- function(n, rho){
  r = (n - 2*rho - n*rho**2 +2*rho**(n+1))/(n*(1-rho)**2)
  ne = n/r
  return(ne)
}

ESS.ma <- function(n, theta){
  rho1 = theta/(1+theta**2)
  r = 1+2*(1-1/n)*rho1
  ne = n/r
  return(ne)
}

ESS.arma <- function(n, rho, theta){
  rho1 = (1+rho*theta)*(rho+theta)/(1+theta**2 + 2*theta*rho)
  r = 1 + (2*rho1*(rho^n - n*rho +n-1))/(n*(rho-1)**2)
  ne = n/r
  return(ne)
}

# This function is used to simulate time series 
# input: auto, hetero, sigma, arma.model, length N, burn.in, 
# sigma of heteroskedatic case inclus N
# sigma in the monthlyvar = 1 case inclus 12 elements 
simulate.general <- function(N, arma.model, burn.in = 0, 
                             hetero, sigma, monthly.var = 1,
                             gaps = 0, 
                             outlier.model = 0, prob.outliers = 0, size.outliers = 0){
  if(burn.in > 0){
    N = N + burn.in
  }
  
  if(hetero == 1){
    sigma <- rep(sigma, times = round(N/length(sigma)))
  }
  
  if (monthly.var == 1){
    number.year = round(N/365.25)
    length.month = c(31,28,31,30,31,30,31,31,30,31,30,31)
    month.name = rep( rep(c(1:12), times = length.month), times = number.year)
    N = length(month.name)
  }
  # initiate noise and simulated series 
  Y.sim <- rep(NA, N)
  Y.sim[1] <- 0
  monthly.noise <- rep(NA, N)
  ini.noise =  rnorm(n = N, mean = 0, sd = 1)
  
  # considering the autocorrelation
  alpha =  arma.model[1]
  beta = arma.model[2]
  coef = (1+2*alpha*beta+beta**2)/(1-alpha**2)
  sigma = sqrt((sigma**2)/coef)

  if(monthly.var == 1){
    for (i in c(1:12)) {
      ind = which(month.name == i)
      monthly.noise[ind] <- ini.noise[ind] * sigma[i]
    }
  }else{
    monthly.noise <- ini.noise * sigma
  }
  
  for ( j in 2:N){
    Y.sim[j] <- Y.sim[j-1]*alpha + monthly.noise[j] + monthly.noise[j-1]*beta
  }
  sim.series = Y.sim[(burn.in+1):N] 
  
  n = length(sim.series)
  cluster.true <- rep(1, n)
  if (outlier.model == 1) {
    outliers = sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0)
    sim.series[pos.outliers] <-  outliers[which(outliers!=0)]
    
    cluster.true[pos.outliers]=2
  } else if(outlier.model == 2) {
    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0 )
    outliers.amp <- rnorm(length(Y[pos.outliers]), mean = size.outliers, sd = 1)
    sim.series[pos.outliers] <- outliers.amp
    
    cluster.true[outliers>0]=2
    cluster.true[outliers<0]=3
  }
  
  return(sim.series)
}


simulate.general1 <- function(N, arma.model, burn.in = 0, 
                             hetero, sigma,
                             gaps = 0, 
                             outlier.model = 0, prob.outliers = 0, size.outliers = 0){
  # if(hetero == 1){
  #   sigma <- rep(sigma, times = round(N/length(sigma)))
  # }
  
  # initiate noise and simulated series 
  # Y.sim <- rep(NA, N)
  # Y.sim[1] <- 0
  # monthly.noise <- rep(NA, N)
  # ini.noise =  rnorm(n = N, mean = 0, sd = 1)
  # 
  # considering the autocorrelation
  alpha =  arma.model[1]
  beta = arma.model[2]
  if(alpha!=0 | beta!=0){
    coef = (1+2*alpha*beta+beta**2)/(1-alpha**2)
    sigma = sqrt((sigma**2)/coef)
    x.ini = arima.sim(model = list(ar = alpha, ma = beta), n = N, n.start = burn.in, sd = 1)
   
  }else{
    x.ini = rnorm(n = N, mean =0, sd = 1)
  }
  
  sim.series = x.ini*sigma
  return(sim.series)
}


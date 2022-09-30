# This function is used to simulate time series 
# input: auto, hetero, sigma, arma.model,  
if (yearly ==1){
  
}
if(auto == 1){
  alpha =  arma.model[1]
  beta = arma.model[2]
}
monthly = length(sigma)

if (monthly == 1 & (alpha+beta != 0)){
  sim.series <- arima.sim(model = list(ar = alpha, ma = beta), n = N, sd = sigma.1)
} else if(monthly == 1 & (alpha+beta == 0)){
  sim.series <- rnorm(n = N, sd = sigma.1, mean = 0)
}else if(monthly != 1){
  number.year = N/(sum(length.month))
  Y.sim <- rep(NA, N)
  month.name = rep( rep(c(1:monthly), times = length.month), times = number.year)
  # generate monthly noise 
  ini.noise =  rnorm(n = N, mean = 0, sd = 1)
  monthly.noise <- rep(NA, N)
  coef = (1+2*alpha*beta+beta**2)/(1-alpha**2)
  sigma.2 = sqrt((sigma.1**2)/coef)
  for (i in 1:monthly) {
    ind = which(month.name == i)
    monthly.noise[ind] <- ini.noise[ind] * sigma.2[i]
  }
  Y.sim[1] <- mean.1
  for ( j in 2:N){
    Y.sim[j] <- Y.sim[j-1]*alpha + monthly.noise[j] + monthly.noise[j-1]*beta
  }
  sim.series = Y.sim
}
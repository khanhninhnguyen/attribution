# test the screening methods 
gaps.list <- c(5,10,15)
k = 1

n0 = 100
sim.ar <- arima.sim(model = list(ar = 0.32), n = n0, sd = 1)
ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
sim.ar <- sim.ar + ind.out*sign(sim.ar)
Q1 <- quantile(sim.ar, .25)
Q3 <- quantile(sim.ar, .75)
IQR <- IQR(sim.ar)
up = Q3+IQR*1.5
down = Q1-IQR*1.5
a = sim.ar[which(sim.ar>down & sim.ar<up)]
IQR/1.349
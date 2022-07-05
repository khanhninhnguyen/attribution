# draft -------------------------------------------------------------------

# test the scaling factor between the variance of the difference a --------

phi = 0.8
theta = -0.5
var0 = 1
var_x = (1 + 2*phi*theta+theta**2)*var0/(1-phi**2)
print(var_x)
scale.diff = 2*(1-phi)*(1+theta**2-theta+phi*theta)/((1 + 2*phi*theta+theta**2))
scale.ma = 2*(1+ theta**2-theta)/(1+theta**2)
nb.sim = 1000
var.est = rep(NA, nb.sim)
var.diff.est = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  sim.series <- arima.sim(model = list(ar = phi, ma = theta), n = 1000, sd = sqrt(var0))
  sim.diff <- diff(sim.series)
  var.est[i] <- var(sim.series)
  var.diff.est[i] <- var(sim.diff)/scale.diff
}
summary(var.est)
summary(var.diff.est)


# test the impact of estimate 2n and the mean -----------------------------
n = 200
nb.sim = 1000
var.2n = rep(NA, nb.sim)
var.mean = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  x = rnorm(n, mean = 0, sd = 1)
  var.2n[i] = var(x)
  var.mean[i] = (var(x[1:100]) + var(x[101:200]))/2
}
summary(var.2n)
summary(var.mean)
var(var.2n)
var(var.mean)


nearby_search_distance(coor.file, list.name, list.gnss, horizontal, vertical, version_name, nearby.ver)
  
a = read.table(file = validation.file.ref, header = TRUE)

list.all =list.files(paste0(path_homo,"34/"))
name.full = substr(list.all,start = 6, stop = 9)

name.full %in% list.nearby.station.homo





# check the parameters from ARMA ------------------------------------------

res <- data.frame(matrix(NA, ncol = 3, nrow = 100))
for (i in c(1:100)) {
  x = arima.sim(model = list(ar=0.5), n = 1000, sd=1) + 1
  fit = arima(x, order = c(1,0,0), include.mean = TRUE)
  res[i,] <- c(mean(x), fit$coef)
}

summary(res$X2)
summary(res$X3)

yt <- arima.sim(list(order=c(1,0,0), ar=.5), n=500)
xt <- yt + 10   


# impact of gaps on the param ---------------------------------------------

nb.sim = 1000
gaps.list = seq(10, 80, 20)
n0 = 360
sdv = 1
Res1 <- list()
for (k in 1:length(gaps.list)) {
  # ar.order <- data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
  # ma.order <- data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
  arma.order <- data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
  for (i in 1:nb.sim) {
    # sim.ar <- arima.sim(model = list(ar = 0.3), n = n0, sd = sdv)
    # sim.ma <- arima.sim(model = list(ma = 0.3), n = n0, sd = sdv)
    sim.arma <- arima.sim(model = list(ar = 0.5, ma = 0), n = n0, sd = sdv)
    gaps =  rbinom(n0, 1, gaps.list[k]/100)
    ind.gaps = which(gaps != 0)
    # sim.ar[ind.gaps] <- NA
    # sim.ma[ind.gaps] <- NA 
    sim.arma[ind.gaps] <- NA
    
    # arfit2 = fit.arima(sim.ar, select.sig = 1)
    # mafit2 = fit.arima(sim.ma, select.sig = 1)
    # armafit2 = fit.arima(sim.arma, select.sig = 1)
  
    # ar.order[i,] <- 
    # ma.order[i,] <- c(unlist(mafit2))
    arma.order[i,c(1,2)] <- arima(sim.arma, order = c(1,0,0), method = "ML")$coef
    arma.order[i,3] <- classic(sim.arma)
  }
  # ar.2 = ar.order[which(ar.order$X3==1 & ar.order$X4==0),]
  # ma.2 = ma.order[which(ma.order$X3==0 & ma.order$X4==1),]
  # arma.2 = arma.order[which(arma.order$X3==1 & arma.order$X4==1),]
  # 
  # res = c(nrow(ar.2), nrow(ma.2), nrow(arma.2))
  # names(res) <- c("ar1", "ma1", "arma11")
  Res1[[k]] <-arma.order
} 

Res1.df = bind_rows(Res1)
b <- data.frame( value = Res1.df$X3, gaps = rep(as.factor(gaps.list), each = nb.sim ))

jpeg(paste0(path_results,"figure/ATTRIBUTION-TEST/" , "per-atuoarima-gaps-phi0.5.jpeg"),
     width = 3000, height = 2000,res = 300) # change name
ggplot(b, aes(x = gaps, col = gaps, y =value))+
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)  +
# scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0,1))+
  geom_hline(yintercept =0.5) +
  labs(y = "phi", x = "Gaps percentage(%)")+ 
  theme_bw()
dev.off()

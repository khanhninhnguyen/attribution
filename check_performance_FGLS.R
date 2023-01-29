# check efficiency of the FGLS
list.param.ar = seq(0,0.9,0.15)
list.param.sig = seq(0.1, 0.45, 0.05)

hetero = 0
autocor = 1
y.axis = "FPR"
jump = ifelse(y.axis=="FPR", 0, 0.3)
noise.name = "ARMA(1,1)"
k = 1
x.axis = "rho"
n =100
# true.sim.sd = sqrt(0.5 - list.param.sig[k]*cos(2*pi*(c(1:730)/365) -pi))
# ar = 0.3

true.sim.sd = sqrt(0.5 - list.param.sig[4]*cos(2*pi*(c(1:730)/365) -pi))
ar = list.param.ar[k]
# if(hetero==0){
#   if(autocor==1){
#     ar = list.param.ar[k]
#   }
#   true.sim.sd = sqrt(0.5)
# }else{
#   if(autocor==1){
#     ar = 0.3
#   }else{
#     if(x.axis == "rho"){
#       true.sim.sd = sqrt(0.5 - list.param.sig[4]*cos(2*pi*(c(1:730)/365) -pi))
#     }else{}
#     ar =0
#     x.axis = "sig"
#     noise.name = "white"
#     true.sim.sd = sqrt(0.5 - list.param.sig[k]*cos(2*pi*(c(1:730)/365) -pi))
#     # 
#   }
# }
Res = get(load(file = paste0(path_results,"attribution/N100/",hetero,"auto",autocor, x.axis, y.axis,noise.name,"R.Data")))

true.sd.beta = sqrt(Res$total[[k]][[1]]$gls$vcov[1,1])

true.sd.beta

# Check the variance estimates 

if(hetero!=0){
  fgls.sd = rowMeans(sapply(c(1:n), function(x) sqrt(Res$total[[k]][[x]]$fgls$var)))
}else{
  fgls.sd = colMeans(sapply(c(1:n), function(x) sqrt(Res$total[[k]][[x]]$fgls$var)))
}

# Check the coefficients estimates, phi, theta
fgls.beta = sapply(c(1:n), function(x) Res$total[[k]][[x]]$fgls$t.table$Estimate[1])
# Check the variance of beta 
fgls.sd.beta = sapply(c(1:n), function(x) Res$total[[k]][[x]]$fgls$t.table$`Std. Error`[1])
# Check the t-statistic 
fgls.t = sapply(c(1:n), function(x) Res$total[[k]][[x]]$fgls$t.table$`t value`[1])
# plot sd (offset)
jpeg(paste0(path_results,"attribution/FGLS",hetero,autocor,x.axis, y.axis,k,noise.name="arma1", ".jpg"), width = 800, height = 600)

par(mfrow=c(2,2))

mean.est = mean(fgls.sd.beta)
sd.est = sd(fgls.sd.beta)
hist(fgls.sd.beta, breaks = 100, xlab = "SD(jump)",
     main = paste0("True = ",round(true.sd.beta, digits = 3), ", Est = ", round(mean.est, digits = 3), ", SD = ", round(sd.est, digits = 3)))
abline(v=true.sd.beta, col ="red")
abline(v=mean.est, col ="blue")

# plot var estimates
if(hetero==0){
  mean.est = mean(fgls.sd)
  sd.est = sd(fgls.sd)
  hist(fgls.sd, breaks = 100, xlab = "SD",
       main = paste0("True = ",round(true.sim.sd, digits = 3), ", Est = ", round(mean.est, digits = 3), ", SD = ", round(sd.est, digits = 3)))
  abline(v=true.sim.sd, col ="red")
  abline(v=mean.est, col ="blue")
}else{
  plot(true.sim.sd, col="blue", type = "l", main =" simulated SD")
  lines(fgls.sd, col="red", type = "l")
}

# plot jump estimate 
mean.est = mean(fgls.beta)
sd.est = sd(fgls.beta)
hist(fgls.beta, breaks = 100, xlab = "jump",
     main = paste0("True = ",round(jump, digits = 3), ", Est = ", round(mean.est, digits = 3), ", SD = ", round(sd.est, digits = 3)))
abline(v=jump, col ="red")
abline(v=mean.est, col ="blue")

# plot t statistics 

mean.est = mean(fgls.t)
sd.est = sd(fgls.t)
hist(fgls.t, breaks = 100, xlab = "t-value",
     main = paste0("True = ",round((jump/true.sd.beta), digits = 3), ", Est = ", round(mean.est, digits = 3), ", SD = ", round(sd.est, digits = 3)))
abline(v=(jump/true.sd.beta), col ="red")
abline(v=mean.est, col ="blue")
dev.off()

table(fgls.t>1.96)
# check value 
fgls.sd


# FOR PAPER ---------------------------------------------------------------
rate.ext <- function(Res, list.sample, off.set){
  tot.tpr = data.frame(matrix(NA, ncol = 4, nrow = length(list.sample)))
  for (j in c(1:length(list.sample))) {
    #FGLS
    fgls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$Estimate[1])
    mu.fgls = fgls.beta + off.set
    fgls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$`Std. Error`[1])
    SD.jump.fgls = fgls.sd.beta 
    t.values.fgls = mu.fgls/SD.jump.fgls
    tpr.fgls = length(which(abs(t.values.fgls) > 1.96))
    
    # GLS 
    gls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$Estimate[1])
    mu.gls = gls.beta + off.set
    gls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$`Std. Error`[1])
    SD.jump.gls = gls.sd.beta 
    t.values.gls = mu.gls/SD.jump.gls
    tpr.gls = length(which(abs(t.values.gls) > 1.96))
    
    #OLS
    ols.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.ols$Estimate[1])
    mu.ols = ols.beta + off.set
    ols.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.ols$`Std. Error`[1])
    SD.jump.ols = ols.sd.beta 
    t.values.ols = mu.ols/SD.jump.ols
    tpr.ols = length(which(abs(t.values.ols) > 1.96))
    
    #HAC
    hac.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.hac$`Std. Error`[1])
    SD.jump.hac = hac.sd.beta 
    t.values.hac = mu.ols/SD.jump.hac
    tpr.hac = length(which(abs(t.values.hac) > 1.96))
    
    tot.tpr[j,] = c(tpr.fgls, tpr.gls, tpr.ols, tpr.hac)
  }
  colnames(tot.tpr) = c("FGLS", "GLS", "OLS", "OLS-HAC")
  out = tot.tpr[,c("GLS","FGLS", "OLS-HAC", "OLS")]
  return(out)
  
}

list.param.ar = seq(0,0.9,0.15)
list.param.sig = seq(0, 0.8, 0.2)

N = 400
nb.sim=1000
hetero = 1
autocor = 1
y.axis = "FPR"
# x.axis = "rho"
extrem = 1
x.axis = "sig"
# list.sample = list.param.ar
list.sample = list.param.sig

noise.name = ifelse(autocor==0, "white", "ar1")
Res = get(load(file = paste0(path_results,"attribution/N1000n/",hetero,"auto",autocor, x.axis, y.axis,noise.name,extrem,"R.Data")))


rate.ext(Res, list.sample = list.sample, off.set = 0)

FPR = rate.ext(Res, list.sample = list.sample, off.set = 0)/nb.sim
TPR = rate.ext(Res, list.sample = list.sample, off.set = 0.4)/nb.sim
if(x.axis=="rho"){
  FPR$phi = list.param.ar
  TPR$phi = list.param.ar
  dat.fpr = reshape2::melt(FPR, id ="phi")
  dat.tpr = reshape2::melt(TPR, id ="phi")
  
}else{
  FPR$range = list.param.sig*100
  TPR$range = list.param.sig*100
  dat.fpr = reshape2::melt(FPR, id ="range")
  dat.tpr = reshape2::melt(TPR, id ="range")
}

# plot FPR 
unit1 = "cm"
unit2 = "pt"
unit3 = "null"
x.axis.name = ifelse(x.axis=="rho", "phi","range")

y.axis.name = "FPR"
thres = 0.05
if(x.axis.name =="phi"){
  break1 = list.sample
}else{
  break1 = list.sample*100
}
p2 <- eval(parse(
  text=paste0("ggplot(dat.fpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + 
  xlab(x.axis.name) + geom_hline(yintercept =", thres,", lwd = 0.2)+
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=4),
        axis.title = element_text(size = 6), legend.key.size = unit(0.2, unit1), plot.tag = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0,unit2 ), plot.margin = rep(unit(0,unit3),4))")))

ggsave(paste0(path_results,"attribution/FGLS_HAC/h",hetero,"auto",autocor, x.axis.name, y.axis.name,noise.name,N,extrem,".jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)


# plot TPR


y.axis.name = "TPR"
thres = 0.95
p2 <- eval(parse(
  text=paste0("ggplot(dat.tpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + 
  xlab(x.axis.name) + geom_hline(yintercept =", thres,", lwd = 0.2)+
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=4),
        axis.title = element_text(size = 6), legend.key.size = unit(0.2, unit1), plot.tag = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0,unit2 ), plot.margin = rep(unit(0,unit3),4))")))

ggsave(paste0(path_results,"attribution/FGLS_HAC/h",hetero,"auto",autocor, x.axis.name, y.axis.name,noise.name,N,extrem,".jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)



# test -------------------- -----------------------------------------------


one.year = 50
nb.sim = 1000
n = 100

periodic = cos(2*pi*(c(1:n)/(n/2)) )
sigma.0 = (0.7 - periodic*0.5)
sigma.0 = rep(1,n)
var0 = sigma.0^2

testr = rep(NA, nb.sim)
testw = rep(NA, nb.sim)
sd = rep(NA, nb.sim)
for (i in c(1:nb.sim)) {
  set.seed(i)
  y = simulate.general1(N = n, arma.model = c(0,0), burn.in = burn.in, hetero = 1, sigma = (sigma.0))
  sd[i] = sd(y)
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:n])
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = n/2, one.year = one.year)
  Data.mod = Data.mod[,c(1,10,11)]
  fit.wls = lm(signal~right, data = Data.mod, weights=(sigma.0))
  wls.test = lmtest::coeftest(fit.wls,df=(n-3)) [, ] %>% as.data.frame()
  a = GLS(0,0, var0, Data.mod)
  testr[i] = a$t.table$`t value`[1]
  testw[i] = wls.test$`t value`[2]
  print(i)
}

X= as.matrix(Data.mod[,c(2,3)])
solve(t(X) %*% solve(diag(sigma.0)) %*% X)

table(abs(testr)<1.962)

table(abs(testw)<1.96)
sd(testr)



# power of test  ----------------------------------------------------------



compute_power <- function(Res, name.test, k, off.set, df1){
  nb.sim = length(Res$total[[1]])
  beta =  sapply(c(1:nb.sim), function(x) Res$total[[k]][[x]][[name.test]]$t.table$Estimate[1])
  sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[k]][[x]][[name.test]]$t.table$`Std. Error`[1])
  out = rep(NA, length(off.set))
  cv = qt(c(.975), df=df1) 
  for (i in c(1:length(off.set))) {
    t.values = (beta + off.set[i])/sd.beta
    out[i] = length(which(abs(t.values)>cv))
  }
  return(out)
}


n = c(200,400,600,800)
list.offset =  seq(0, 1, 0.02)
tpr = data.frame(matrix(NA, ncol = 4, nrow = length(list.offset)))
h = 0
nb.sim=100
for (i in c(1:4)) {
  n1 = n[i]
  b = get(load(
    file = paste0(path_results,"attribution/",hetero=h,"auto",autocor=1, x.axis="rho", y.axis = "FPR",
                  noise.name = "ar1",n=n1,"1R.Data")))
  tpri = compute_power(Res = b, name.test = "fgls", k = 1, off.set = list.offset, df1 = n-3)
  tpr[,i] = tpri 
}

colnames(tpr) = paste0("N=",n)
tpr$jump = list.offset 

dat.tpr = reshape2::melt(tpr, id ="jump")

unit1 = "cm"
unit2 = "pt"
unit3 = "null"

thres = 0.95
p2 <- ggplot(dat.tpr, aes(x = jump,y = value/nb.sim, col = variable))+
  geom_line(lwd = 0.3) +theme_bw() +
  ylab("Power") + 
  xlab("Jump") + geom_hline(yintercept =0.95, lwd = 0.2)+
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  # scale_x_continuous(breaks)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=4),
        axis.title = element_text(size = 6), legend.key.size = unit(0.2, unit1), plot.tag = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0,unit2 ), plot.margin = rep(unit(0,unit3),4))

ggsave(paste0(path_results,"attribution/power/h",h,"auto",autocor=1, "jump", "power",".jpg" ), 
       plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)




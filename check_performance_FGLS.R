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
list.param.ar = seq(0,0.9,0.15)
list.param.sig = seq(0, 0.8, 0.2)

nb.sim=1000
hetero = 1
autocor = 1
y.axis = "FPR"
x.axis = "rho"
# x.axis = "sig"
noise.name = ifelse(autocor==0, "white", "ar1")

Res = get(load(file = paste0(path_results,"attribution/N1000n/",hetero,"auto",autocor, x.axis, y.axis,noise.name,"1R.Data")))

rate.ext <- function(Res, list.sample, off.set){
  tot.tpr = data.frame(matrix(NA, ncol = 4, nrow = length(list.sample)))
  for (j in c(1:length(list.sample))) {
    #FGLS
    fgls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$Estimate[1])
    mu.fgls = fgls.beta + off.set
    fgls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$`Std. Error`[1])
    SD.jump.fgls = fgls.sd.beta 
    t.values.fgls = mu.fgls/SD.jump.fgls
    tpr.fgsl = length(which(abs(t.values.fgls) > 1.96))
    
    # GLS 
    gls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$Estimate[1])
    mu.gls = gls.beta + off.set
    gls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$`Std. Error`[1])
    SD.jump.gls = fgls.sd.beta 
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
    
    tot.tpr[j,] = c(tpr.fgsl, tpr.gls, tpr.ols, tpr.hac)
  }
  colnames(tot.tpr) = c("FGLS", "GLS", "OLS", "OLS-HAC")
  out = tot.tpr[,c("GLS","FGLS", "OLS-HAC", "OLS")]
  return(out)
  
}

list.sample = list.param.ar
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

ggsave(paste0(path_results,"attribution/h",hetero,"auto",autocor, x.axis.name, y.axis.name,noise.name,"1.jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)


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

ggsave(paste0(path_results,"attribution/h",hetero,"auto",autocor, x.axis.name, y.axis.name,noise.name,"1.jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)









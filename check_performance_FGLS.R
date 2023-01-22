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
Res = get(load(file = paste0(path_results,"attribution/N1000/",hetero,"auto",autocor, x.axis, y.axis,noise.name,"R.Data")))

true.sd.beta = sqrt(Res$total[[k]][[1]]$gls$vcov[1,1])

true.sd.beta

# Check the variance estimates 

if(hetero!=0){
  fgls.sd = rowMeans(sapply(c(1:1000), function(x) sqrt(Res$total[[k]][[x]]$fgls$var)))
}else{
  fgls.sd = colMeans(sapply(c(1:1000), function(x) sqrt(Res$total[[k]][[x]]$fgls$var)))
}

# Check the coefficients estimates, phi, theta
fgls.beta = sapply(c(1:1000), function(x) Res$total[[k]][[x]]$fgls$t.table$Estimate[1])
# Check the variance of beta 
fgls.sd.beta = sapply(c(1:1000), function(x) Res$total[[k]][[x]]$fgls$t.table$`Std. Error`[1])
# Check the t-statistic 
fgls.t = sapply(c(1:1000), function(x) Res$total[[k]][[x]]$fgls$t.table$`t value`[1])
# plot sd (offset)
jpeg(paste0(path_results,"attribution/FGLS",hetero,autocor,ar,x.axis, y.axis,k, ".jpg"), width = 800, height = 600)

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

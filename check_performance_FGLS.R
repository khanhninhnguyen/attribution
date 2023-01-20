# check efficiency of the FGLS
hetero = 0
autocor = 1
x.axis = "rho"
y.axis = "TPR"
noise.name = "ar1"
Res = get(load(file = paste0(path_results,"attribution/N1000/",hetero,"auto",autocor, x.axis, y.axis,noise.name,"R.Data")))
k = 4
fgls.sd = sapply(c(1:1000), function(x) Res$total[[k]][[x]]$fgls$t.table$`Std. Error`[1])
hac.sd = sapply(c(1:1000), function(x) Res$total[[k]][[x]]$fit.hac$`Std. Error`[1])
b = Res$total[[k]][[1]]$gls$vcov
sqrt(b)
hist(fgls.sd, breaks = 100)
abline(v=sqrt(b)[1,1], col ="red")
abline(v=mean(fgls.sd), col ="blue")

# 
hist(hac.sd, breaks = 100)
abline(v=sqrt(b)[1,1], col ="red")
abline(v=mean(hac.sd), col ="blue")
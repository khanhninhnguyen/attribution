# this prog used to apply the ancova/fgls on the real data 
library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att, "newUsed_functions.R"))
win.thres = 1
dat  = get(load(file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]

tot.res <- list()
# list.ind = c(1:length(dat)) # for the gps'-era', cqse 295 has error 
# list.ind = list.ind[-295]
data.test = gps.era.dat
list.ind = c(1:length(data.test)) 
for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  # if(which(Y.with.NA$date==date.detected.break)!=one.year) next
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  ########################
  # Test on the OLS estimate with the HAC covariance and selection of the 
  # significant parameters
  res.hac <- Test_OLS_vcovhac(Data.mod)
  ####################################
  # Test on the FGLS estimate with a ARMA(1,1) model
  # res.fgls <- Test_FGLS(Data.mod)
  # 
  tot.res[[name.dataset]] <- res.hac
  # tot.res[[name.dataset]] <- list(hac = round(res.fgls$res.gls, digits = 5),
  #                                 predicted = res.fgls$predicted)
  #                                 fgls = round(res.fgls$res.gls, digits = 5))
  print(k)
}

save(tot.res, file = paste0(path_results, "attribution/",name.series,".mean_test.RData"))

r <- rep(NA, length(tot.res))
for (i in c(1:length(tot.res))) {
  resi = tot.res[[i]]$hac
  r[i] <- resi$`Pr(>|z|)`[ which(rownames(resi) == "Jump")]

}

# plot results to investigate 

plot_HAC <- function(case.name, tot.res, data.in, name.var, ver ){
  Y = data.in[[case.name]]
  res.i = tot.res[[case.name]]
  Y$predict = rep(NA, 365*2)
  Y$predict[which(is.na(Y[name.var]) == FALSE)] <- res.i$predicted
  data.plot <- data.frame(date = Y$date, name.var = Y[name.var], predicted = Y$predict)
  a = reshape2::melt(data.plot, id = "date")
  
  jpeg(paste0(path_results,"attribution/OLS.HAC/", case.name, ver, ".jpeg" ),
       width = 3000, height = 1500,res = 300)
  p <- ggplot(data = a, aes(x = date, y = value, col = variable)) + 
    geom_line()+theme_bw()+
    geom_vline(xintercept = Y$date[365], size = 0.2)+
    ylab(name.var)+
    labs(title = case.name, 
      subtitle =  
        paste(rownames(res.i$fit.hac),  
              round(res.i$fit.hac$Estimate, digits = 4), 
              round(res.i$fit.hac$`Pr(>|z|)`, digits = 4), 
              collapse = ", ", sep = "; "))+ 
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size=15,face="bold"))
  print(p)
  dev.off()
  
}


# investigate the difference between selection and nonselection -----------


selec = get(load(file = paste0(path_results, "attribution/", name.series, ".mean_test.RData")))
full = get(load(file = paste0(path_results, "attribution/GPS1.ERA1.mean_test.RData")))
# full = get(load(file = paste0(path_results, "attribution/GPS.ERA.1step.mean_test.RData")))

hac.sel <- sapply(c(1:length(selec)), function(x){
  res.i = selec[[x]]$fit.hac
  ind = which(rownames(res.i) == "Xt")
  if(length(ind) == 0){
    NA
  }else{
    unlist(res.i$`Pr(>|z|)`[ind])
  }
})

hac.sel <- unlist(hac.sel)
hac.full <- sapply(c(1:length(full)), function(x){
  res.i = full[[x]]$hac
  ind = which(rownames(res.i) == "Jump")
  if(length(ind) == 0){
    NA
  }else{
    unlist(res.i$`Pr(>|z|)`[ind])
  }
})

r = c()
for (j in c(1:(length(unique.ind)-1))) {
  beg=unique.ind[j]
  end=unique.ind[j+1]
  r <- c(r, length(which(hac.full[beg:end]<0.05))/(end-beg))
}

ind.list = which(hac.sel<0.05)

for (j in ind.list) {
  plot_HAC(case.name = names(selec)[j], tot.res = selec, data.in = dat, name.var = name.series, ver = "s")
}


a = Data.mod[c(1:365),(-2)]
summary(lm( signal~., data = a))

a = Data.mod



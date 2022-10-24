# this prog used to apply the ancova/fgls on the real data 
library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att, "newUsed_functions.R"))
dat  = get(load(file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"screened.RData")))
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
    mutate(Jump=c(rep("Left",one.year),rep("Right",one.year))) %>% 
    mutate(complete.time=1:(2*one.year)) %>% 
    mutate(Xt=complete.time-one.year/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  ########################
  # Test on the OLS estimate with the HAC covariance and selection of the 
  # significant parameters
  # res.hac <- Test_OLS_vcovhac_1step(Data.mod)
  ####################################
  # Test on the FGLS estimate with a ARMA(1,1) model
  res.fgls <- Test_FGLS(Data.mod)

  tot.res[[name.dataset]] <- list(hac = round(res.fgls$res.gls, digits = 5),
                                  predicted = res.fgls$predicted)
  #                                 fgls = round(res.fgls$res.gls, digits = 5))
  print(k)
}

save(tot.res, file = paste0(path_results, "attribution/GPS.ERA.gls.mean_test.RData"))

r <- rep(NA, length(tot.res))
for (i in c(1:length(tot.res))) {
  resi = tot.res[[i]]$hac
  r[i] <- resi$`Pr(>|z|)`[ which(rownames(a) == "Jump")]

}

# plot results to investigate 

plot_HAC <- function(case.name, tot.res, data.in, name.var ){
  Y = data.in[[case.name]]
  res.i = tot.res[[case.name]]
  Y$predict = rep(NA, 365*2)
  Y$predict[which(is.na(Y[name.var]) == FALSE)] <- res.i$predicted
  data.plot <- data.frame(date = Y$date, name.var = Y[name.var], predicted = Y$predict)
  a = reshape2::melt(data.plot, id = "date")
  
  jpeg(paste0(path_results,"attribution/OLS.HAC/", case.name, ".jpeg" ),
       width = 3000, height = 1500,res = 300)
  ggplot(data = a, aes(x = date, y = value, col = variable)) + 
    geom_line()+theme_bw()+
    geom_vline(xintercept = Y$date[365], size = 0.2)+
    ylab(name.var)+
    labs(title = case.name, 
      subtitle =  paste(rownames(res.i$hac),  res.i$hac$Estimate, res.i$hac$`Pr(>|z|)`, collapse = ", ", sep = "; "))+ 
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size=15,face="bold"))
  dev.off()
  
}





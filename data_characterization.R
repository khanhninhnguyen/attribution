#### This function is used for the data characterization. This is step after pairing and screening data ####

# Heteroskedastic ---------------------------------------------------------

# from OLS resdiual 
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]
n = length(data.test)
sd.all <- rep(NA, n)
data.test = gps.era.dat
for (k in c(1:n)) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
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
  res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  sd.all[k] = sd(res.hac.1step$fit.ols$residuals)
}
# autocorrelation  --------------------------------------------------------



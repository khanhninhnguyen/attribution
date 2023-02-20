# this program is to explore the test result 
source(paste0(path_code_att,"plot_sic_diff_series.R"))

Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
all.var = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
all.arima = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)

# catch the problematic cases  --------------------------------------------

## case of G-G' = 0 --------------------------------------------------------
data.GG = Total.res[which(Total.res$`G-G'`==0),]
all.case.p = data.GG$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}
### variance plot 
dat.p = data.frame(mean.msd = (all.var$`mean.G-G'`)/(all.var$`mean.G-E`), sig = as.factor(Total.res$`G-G'`))
ggplot(data = dat.p, aes(x = sig, y = mean.msd))+theme_bw()+
  geom_boxplot()
### length plot  
dat.p = data.frame(n = Total.res$n1+ Total.res$n2, sig = as.factor(Total.res$`G-G'`))
ggplot(data = dat.p, aes(x = sig, y = n))+theme_bw()+
  geom_boxplot()

## case of E-E' != 0 --------------------------------------------------------
data.EE = Total.res[which(Total.res$`E-E'`!=0),]
all.case.p = data.EE$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}

dat.p = data.frame(mean.msd = (all.var$`mean.E-E'`)/(all.var$`mean.G-E`), sig = as.factor(Total.res$`E-E'`))
ggplot(data = dat.p, aes(x = sig, y = mean.msd))+theme_bw()+
  geom_boxplot()

dat.p = data.frame(n = Total.res$n1+ Total.res$n2, sig = as.factor(Total.res$`E-E'`))
ggplot(data = dat.p, aes(x = sig, y = n))+theme_bw()+
  geom_boxplot()

## case of G'-E' != 0 --------------------------------------------------------
data.G1E1 = Total.res[which(Total.res$`G'-E'`!=0),]
all.case.p = data.G1E1$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}

dat.p = data.frame(mean.msd = (all.var$`mean.G'-E'`)/(all.var$`mean.G-E`), sig = as.factor(Total.res$`G'-E'`))
ggplot(data = dat.p, aes(x = sig, y = mean.msd))+theme_bw()+
  geom_boxplot()

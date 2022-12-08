# estimation for specific station 
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

# read data for station 
cas.ind = which(names(dat) == "hers.2012-03-15.bar3")
case.name = names(dat)[cas.ind]
# construct the model 

Y.with.NA = dat[[case.name]]
name.series = "gps.era"
Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
  rename(signal=name.series) %>% 
  mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>%
  mutate(complete.time=1:(2*one.year*win.thres)) %>%
  mutate(Xt=complete.time-one.year*win.thres/2) %>% 
  dplyr::select(-date)
for (i in 1:4){
  eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
}
Data.mod <- Data.mod %>% dplyr::select(-Xt,-complete.time)

# OLS 
ols.fit = lm(signal~., data = Data.mod)
fit.ols=lmtest::coeftest(ols.fit,df=ols.fit$df.residual)[, ] %>% as.data.frame()
print(fit.ols[2,])

# check the noise model 
residus = get(load(file = paste0(path_results,"attribution/norm.residual.ols.RData")))
dat.i = residus[[name.series]][[cas.ind]]
bef.residus = dat.i[1:(one.year*10)]
aft.residus = dat.i[-c(1:(one.year*10))]
fit.arima(bef.residus)
fit.arima(aft.residus)

# OLS-HAC 
vcov.para=sandwich::kernHAC(ols.fit, prewhite = TRUE, approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
fit.hac=lmtest::coeftest(ols.fit,df=(ols.fit$df.residual),vcov.=vcov.para)[, ] %>% as.data.frame()
print(fit.hac[2,])
# FGLS
# read the moving variance 
sd.all= get(load( file = paste0(path_results,"attribution/sd.all_",  win.thres,"years_", nearby_ver,"screened.RData")))
names(sd.all)[cas.ind]
d = sd.all[[cas.ind ]]
weight1 = c(unlist(d$bef[name.series]), unlist(d$aft[name.series]))
Data.mod$weight1 = weight1^2
gls.fit <- gls(signal~Jump+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,
               data=Data.mod,
               correlation =  corAR1(form = ~ 1),
               na.action=na.omit,
               weights=varFixed(value = ~weight1))
vcov.para = gls.fit$varBeta
fit.gls=lmtest::coeftest(gls.fit,df=(ols.fit$df.residual),vcov.=vcov.para)[, ] %>% as.data.frame()
print(fit.gls[2,])

Y.with.NA$var = weight1
data.wt.na = Y.with.NA[which(is.na(Y.with.NA$gps.era)==FALSE),]
min.dat = min(data.wt.na$date)
max.dat = max(data.wt.na$date)
data.lim = Y.with.NA[which(Y.with.NA$date >= min.dat & Y.with.NA$date <= max.dat),]

jpeg(paste0(path_results,"attribution/std1.jpg" ),width = 2600, height = 1800,res = 300)
p <- ggplot(data.lim, aes(x=date, y=var)) + 
  geom_line(col="black")+theme_bw()+ 
  geom_vline(xintercept = Y.with.NA$date[3650], linetype = "dashed")+
  xlab("") + ylab("Moving standard deviation of GPS-ERA")+
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=20))

print(p)
dev.off()








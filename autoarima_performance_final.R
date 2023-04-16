# Final prog to study the auto.arima to put into the thesis
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

# SIMULATION ---------------------------------------------------------------
## AR(1) model -------------------------------------------------------------

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_AR1.RData"))

## MA(1) model -------------------------------------------------------------
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=0,ma=ar0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_MA1.RData"))
## ARMA(1,1) model -------------------------------------------------------------

length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.3-ar0
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData"))




# PLOT
arma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData")))
ar = get(load(file = paste0(path_results, "attribution0/performance_autoarima_AR1.RData")))
ma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_MA1.RData")))

get_data <- function(list.ini, param.val, true.model, details){
  nb.sim = length(list.ini[[1]][[1]])
  # choose details case at specific length 
  if(details!=0){
    tot.df = data.frame(matrix(NA, ncol = length(list.ini[[1]]), nrow = nb.sim))
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[5]][[j]][[x]]$pq))
      tot.df[,j] = model.est
    }
  }else{
    
    if(param.val == 0){
      tot.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
      for (i in c(1:length(list.ini))) {
        for (j in c(1:length(list.ini[[1]]))) {
          model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[i]][[j]][[x]]$pq))
          y = length(which(model.est == true.model))
          tot.df[j,i] = y
        }
      }
    }else{ # could be filtered the case of true model 
      tot.df = data.frame(matrix(NA, ncol = 2*length(list.ini), nrow = length(list.ini[[1]])))
      for (i in c(1:length(list.ini))) {
        for (j in c(1:length(list.ini[[1]]))) {
          phi = sapply(c(1:nb.sim), function(x) list.ini[[i]][[j]][[x]]$coef[1])
          theta = sapply(c(1:nb.sim), function(x) list.ini[[i]][[j]][[x]]$coef[3])
          tot.df[j,(2*i-1)] = mean(phi, na.rm = TRUE)
          tot.df[j,(2*i)] = mean(theta, na.rm = TRUE)
        }
      }
    }
  }
  return(tot.df)
}
# Fig 1: TPR as a fc of length for different coefs
model.plot = "ARMA(1,1)"
model.data = arma
param.name = "phi"
selected.ind = c(2:7)
df = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = 0)[selected.ind,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = as.factor(rep(coef.list[selected.ind], each = 10))) %>% 
  mutate(TPR = value/nb.sim)

p = ggplot(data = df, aes(x = N, y = TPR, col = phi)) +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lwd = 0.25, col = "black")+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 4.5),
        legend.title=element_text(size = 4.5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(0, "pt"),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5))+
  guides(color=guide_legend(title = param.name)) 
# 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.N.", model.plot, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

# Fig 2: TPR as a fc of coefs for different length

model.plot = "ARMA(1,1)"
model.data = arma
param.name = "N"
selected.ind = c(2:7)
df = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = 0)[selected.ind,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = as.factor(length.list)) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = rep(coef.list[selected.ind], each = 10)) %>% 
  mutate(TPR = value/nb.sim)

p = ggplot(data = df, aes(x = phi, y = TPR, col = N)) +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lwd = 0.25, col = "black")+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  scale_x_continuous(breaks = coef.list[selected.ind], 
                     limits = c(0.1, 0.61))+ 
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 3),
        legend.title=element_text(size = 4.5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(0, "pt"),
        legend.title.align=0.5,
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.y = unit(-0.1, 'cm'),
        plot.subtitle = element_text(size = 5))+
  guides(color=guide_legend(title = param.name, nrow = 5, byrow = TRUE))
# 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.coef.", model.plot, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)




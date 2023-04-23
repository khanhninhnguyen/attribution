# Final prog to study the auto.arima to put into the thesis
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)
nb.sim = 1000
burn.in = 1000
hetero = 0
sigma.sim = 1
selected.ind = c(2:7)

# SIMULATION ---------------------------------------------------------------
## AR(1) model -------------------------------------------------------------

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
## ARMA(1,1)a model -------------------------------------------------------------

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




## ARMA(1,1)b model -------------------------------------------------------------

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.5-ar0
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

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_ARMA1b.RData"))

## Accuracy of the coeff estimation ----------------------------------------
### when detect the true model      ----------------------------------------
#### AR(1)                          ----------------------------------------
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
      fit.ar = arima( y.ar, order = c(1,0,0), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_AR1_true.RData"))
#### MA(1)                          ----------------------------------------
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
      fit.ar = arima( y.ar, order = c(0,0,1), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_MA1_true.RData"))
#### ARMA(1,1)a                      ----------------------------------------
set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.3 - ar0
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = arima( y.ar, order = c(1,0,1), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_ARMA1a_true.RData"))

#### ARMA(1,1)b                      ----------------------------------------
set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.5 - ar0
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = arima( y.ar, order = c(1,0,1), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_ARMA1b_true.RData"))


### when detect the wrong model     ----------------------------------------
#### AR(1) but detect MA(1)         ----------------------------------------
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
      fit.ar = arima( y.ar, order = c(0,0,1), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_AR1_wrong.RData"))

#### MA(1) but detect AR(1)         ----------------------------------------
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
      fit.ar = arima( y.ar, order = c(1,0,0), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_MA1_wrong.RData"))

#### ARMA(1,1) but detect MA(1)         ----------------------------------------
set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.3 - coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = arima( y.ar, order = c(0,0,1), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_ARMA1a_wrong_MA.RData"))

#### ARMA(1,1) but detect AR(1)         ----------------------------------------
set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.3 - coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = arima( y.ar, order = c(1,0,0), method="ML")
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_arima_ARMA1a_wrong_AR.RData"))


# PLOT ---------------------------------------------------------------
arma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData")))
armab = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1b.RData")))

ar = get(load(file = paste0(path_results, "attribution0/performance_autoarima_AR1.RData")))
ma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_MA1.RData")))

get_data <- function(list.ini, param.val, true.model, details){
  nb.sim = length(list.ini[[1]][[1]])
  # choose details case at specific length 
  if(details!=0){
    tot.df = data.frame(matrix(NA, ncol = length(list.ini[[1]]), nrow = nb.sim))
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[details]][[j]][[x]]$pq))
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
get_data1 <- function(list.ini, param.val, true.model, details){
  nb.sim = length(list.ini[[1]][[1]])
  # choose details case at specific length 
  if(details!=0){
    tot.df = data.frame(matrix(NA, ncol = length(list.ini[[1]]), nrow = nb.sim))
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[details]][[j]][[x]]$coef))
      tot.df[,j] = model.est
    }
  }else{
    
    if(param.val == 0){
      tot.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
      for (i in c(1:length(list.ini))) {
        for (j in c(1:length(list.ini[[1]]))) {
          model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[i]][[j]][[x]]$coef))
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
## Fig 1: TPR as a fc of length for different coefs-------------------
model.plot = "ARMA(1,1)"
model.data = armab
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
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.N.", model.plot, "b.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

## Fig 2: TPR as a fc of coefs for different length---------------------

model.plot = "ARMA(1,1)"
model.data = armab
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
 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.coef.", model.plot, "b.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


## Fig 3 detail of the model identification --------------------------
model.plot = "ARMA(1,1)"
model.data = armab
param.name = "phi"
selected.ind = c(2:7)
details = 5
all.model = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = details)[,selected.ind]
ar.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "AR(1)")))
ma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "MA(1)")))
arma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "ARMA(1,1)")))
white.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "White")))

df = data.frame(white = white.est, ar = ar.est, ma = ma.est, arma = arma.est, phi = coef.list[selected.ind]) %>%
  reshape2::melt(id = "phi") %>%
  mutate(predict = as.factor(rep(c( "White", "AR(1)", "MA(1)", "ARMA(1,1)"), each = 6))) %>%
  mutate(phi = as.factor(phi))

p = ggplot(data = df, aes(x = phi, y = value/nb.sim, fill = predict))+
  theme_bw()+
  geom_bar(position="stack", stat="identity", width = 0.5)+
  ylab("Percentage")+
  labs(subtitle = paste0("Model: " , model.plot, ", N = ", length.list[details]))+
  scale_y_continuous(labels = scales::percent) +
  # geom_text(aes(label = paste0(value*100/nb.sim,"%")), 
  #           position = position_stack(vjust = 0.5), size = 1)+
  ggrepel::geom_text_repel(aes(label = paste0(value*100/nb.sim,"%")), 
                           position = position_stack(vjust = 0.5), size = 1, direction = "y", 
                           box.padding = unit(0.01, "lines"))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=4),
        legend.title = element_text(size = 4.5),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        legend.title.align = 0.5,
        plot.subtitle = element_text(size = 4)) +
  guides(color = guide_legend(title = param.name)) 


ggsave(paste0(path_results,"attribution0/auto.arima.TPR.detail.",details, model.plot, "b.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

## Fig 4: accuracy of param estimates
arma = get(load(file = paste0(path_results, "attribution0/performance_arima_MA1_wrong.RData")))
## Fig 4  s.e of estimates as a fc of n for different coef-------------------

get_data_acc <- function(list.ini, phi, theta){
  phi.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
  theta.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
  param = c("", "")
  if(phi == 1){ param[1] = "ar1"}
  if(theta == 1){ param[2] = "ma1"}
  
  for (i in c(1:length(list.ini))) {
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) { list.ini[[i]][[j]][[x]]$coef[param]})
      phi.df[j,i] = sd( model.est[1,], na.rm = TRUE)
      theta.df[j,i] = sd( model.est[2,], na.rm = TRUE)
    }
  }
  tot.df = list(phi = phi.df, theta = theta.df)
  return(tot.df)
}
get_data_acc_mean <- function(list.ini, phi, theta){
  phi.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
  theta.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
  param = c("", "")
  if(phi == 1){ param[1] = "ar1"}
  if(theta == 1){ param[2] = "ma1"}
  
  for (i in c(1:length(list.ini))) {
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) { list.ini[[i]][[j]][[x]]$coef[param]})
      phi.df[j,i] = mean( model.est[1,], na.rm = TRUE)
      theta.df[j,i] = mean( model.est[2,], na.rm = TRUE)
    }
  }
  tot.df = list(phi = phi.df, theta = theta.df)
  return(tot.df)
}

model.plot = "MA1"
model.true = "MA(1)"
model.est = "MA(1)"

if(model.true == model.est){
  est.true = "true"
}else{est.true = "wrong"}

arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true,".RData")))
# arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true, "_", model.est,".RData")))

df.all = get_data_acc(list.ini = arma.acc, phi = 1, theta =1)
### plot for the true model -------------------------------------------------
param.name = "theta"
df = df.all[[param.name ]][selected.ind,]  %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = as.factor(rep(coef.list[selected.ind], each = 10))) 

p = ggplot(data = df, aes(x = N, y = value, col = phi))+
  theme_bw()+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  scale_y_continuous(breaks = seq(0, max(df$value)+0.02,0.02),
                     limits = c(0, max(df$value)))+ 
  ylab(paste0("Standard error of ", param.name))+
  labs(subtitle = paste0("Simulated model: " , model.true, ", Estimate model:", model.est))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=4),
        legend.title = element_text(size = 4.5),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        legend.title.align = 0.5,
        plot.subtitle = element_text(size = 5)) +
  guides(color = guide_legend(title = param.name)) 


ggsave(paste0(path_results,"attribution0/auto.arima.TPR.", model.plot, "_", est.true, param.name, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

### plot for the wrong model -------------------------------------------------
#### AR1 model estimated as MA1 
model.plot = "AR1"
model.true = "AR(1)"
model.est = "MA(1)"

if(model.true == model.est){
  est.true = "true"
}else{est.true = "wrong"}
arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true,".RData")))

df.mean = get_data_acc_mean(list.ini = arma.acc, phi = 1, theta =1)
df.sd = get_data_acc(list.ini = arma.acc, phi = 1, theta =1)
param.name = "theta"
param.name.true = "rho"

df.mean.m = df.mean[[param.name ]][selected.ind,]  %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = as.factor(rep(coef.list[selected.ind], each = 10))) 

df.sd.m = df.sd[[param.name ]][selected.ind,]  %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = as.factor(rep(coef.list[selected.ind], each = 10))) 

df = left_join(df.mean.m, df.sd.m, by =c("N","phi"))

p = ggplot(data = df, aes(x = N, y = value.x, col = phi))+
  theme_bw()+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  # geom_errorbar(aes( ymin = value.x - value.y, ymax = value.x + value.y)) +
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  scale_y_continuous(breaks = seq(0, 0.6, 0.1),
                     limits = c(0,0.6))+
  ylab(paste0("Mean of ", param.name))+
  labs(subtitle = paste0("Simulated model: " , model.true, ", Estimate model:", model.est))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=4),
        legend.title = element_text(size = 4.5),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        legend.title.align = 0.5,
        plot.subtitle = element_text(size = 5)) +
  guides(color = guide_legend(title = param.name.true )) 


ggsave(paste0(path_results,"attribution0/auto.arima.TPR.", model.plot, "_", est.true, param.name, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)







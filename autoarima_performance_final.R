# Final prog to study the auto.arima to put into the thesis
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))
path_arima = paste0(path_results, "attribution0/performance_autoarima/")

length.list = seq(200, 2000, 200)
coef.list = seq(-0.8, 0.8, 0.1)
coef.list1 = seq(0, 0.6, 0.1)

nb.sim = 1000
burn.in = 1000
hetero = 0
sigma.sim = 1
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
        print(i)
        for (j in c(1:length(list.ini[[1]]))) {
          print(j)
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
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      # fit.ar = fit.arima(y.ar)
      fit.ar = forecast::auto.arima( y.ar, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                    max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_AR1_test.RData"))

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
      # fit.ar = fit.arima(y.ar)
      fit.ar = forecast::auto.arima( y.ar, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                     max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_ARMA1a.RData"))




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

## general ARMA(1,1) -------------------------------------------------------------

sum.list = seq(0.1, 0.6, 0.1)
for (s in c(1:length(sum.list))) {
  sum.i = sum.list[s]
  set.seed(1)
  tot.res = list()
  for (l in c(1:length(length.list))) {
    n = length.list[l]
    tot.fit = list()
    for (i in c(1:length(coef.list))) {
      fit.i = list()
      ar0 = coef.list[i]
      ma0 = sum.i-ar0
      print(c(ar0, ma0))
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
  save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_ARMA1", sum.i ,".RData"))
}


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


## general ARMA(1,1) -------------------------------------------------------------

sum.list = seq(0.1, 0.6, 0.1)
for (s in c(1:length(sum.list))) {
  sum.i = sum.list[s]
  set.seed(1)
  tot.res = list()
  for (l in c(1:length(length.list))) {
    n = length.list[l]
    tot.fit = list()
    for (i in c(1:length(coef.list))) {
      fit.i = list()
      ar0 = coef.list[i]
      ma0 = sum.i-ar0
      print(c(ar0, ma0))
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
  save(tot.res, file = paste0(path_results, "attribution0/truemodel_autoarima_ARMA1", sum.i ,".RData"))
}

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
# arma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData")))
sum.param = 0.1
arma = get(load(file = paste0(path_arima, "performance_autoarima_ARMA1",sum.param,".RData")))

ar = get(load(file = paste0(path_arima, "performance_autoarima_AR1.RData")))
ma = get(load(file = paste0(path_arima, "performance_autoarima_MA1.RData")))

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
## Fig 1: TPR as a fc of length for different coefs for AR(1) and MA(1) -------------------
model.plot = "AR(1)"
model.data = ar
param.name = "phi"
coef.list = seq(0, 0.6, 0.1)
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
  scale_y_continuous(breaks = seq(0,1,0.1), 
                     limits = c(0,1))+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5))+
  guides(color=guide_legend(title = param.name)) 
# 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.N.", model.plot, ".jpg" ), plot = p, width = 8, height = 4.8, units = "cm", dpi = 600)
## Fig 1: TPR as a fc of length for different coefs for ARMA(1,1)-------------------
sum.param = 0.1
arma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1",sum.param,".RData")))

model.plot = "ARMA(1,1)"
model.data = arma
param.name = "phi"
theta = sum.param - coef.list
selected.ind = which(!near(coef.list, sum.param, tol = 0.01) & !near(abs(coef.list), 0, tol = 0.01) & !(abs(theta)>0.9) )

# selected.ind = which(!near(coef.list, sum.param, tol = 0.01) & !near(abs(coef.list), 0, tol = 0.01) )[-c(1,2,3)]
# selected.ind = c(10:17)[-3]
# selected.ind = c(2:8)[-5]
# selected.ind = c(2:17)[c(-5, -11)]

coef.list[selected.ind]

df = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = 0)[selected.ind,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = rep(coef.list[selected.ind], each = 10)) %>% 
  mutate(theta = rep(round((sum.param-coef.list[selected.ind]), digits = 2), each = 10)) %>% 
  mutate(lab = paste(phi, theta, sep = ", ")) %>%
  # mutate(shape.p = as.factor(ifelse(((phi*10)%%2)==0, 3, 4))) %>%
  # rowwise() %>% mutate(min.pa = min(abs(phi), abs(theta))) %>% 
  mutate(phi = as.factor(phi)) %>% 
  
  # mutate(theta = as.factor(rep(round((sum.param-coef.list[selected.ind]), digits = 2), each = 10))) %>% 
  # mutate(lab = paste(phi, theta, sep = ", ")) %>%
  mutate(TPR = value/nb.sim)

p = ggplot(data = df, aes(x = N, y = TPR, col = lab)) +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lwd = 0.25, col = "black")+
  geom_point(size = 0.3) +
  geom_line( lwd = 0.25)+
  # labs(color="phi, theta") +  
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  scale_y_continuous(breaks = seq(0,1,0.1), 
                     limits = c(0,1))+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5))+
  guides(color=guide_legend(title = "phi, theta"),
         linetype = FALSE) 
# 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.N.", model.plot,sum.param, "all.jpg" ), plot = p, width = 8, height = 4.8, units = "cm", dpi = 600)


## Fig 2: TPR as a fc of coefs for different length---------------------

model.plot = "AR(1)"
model.data = ar
param.name = "phi"
coef.list = seq(0, 0.6, 0.1)
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
  scale_y_continuous(breaks = seq(0,1,0.1), 
                     limits = c(0,1))+
  xlab(param.name) + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5))
  # guides(color=guide_legend(title = param.name, byrow = TRUE))
 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.coef.", model.plot, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)
### For ARMA(1,1)
sum.param = 0.6
arma = get(load(file = paste0(path_arima, "performance_autoarima_ARMA1",sum.param,".RData")))

model.plot = "ARMA(1,1)"
model.data = arma
param.name = "phi, theta"
theta = sum.param - coef.list
selected.ind = which(!near(coef.list, sum.param, tol = 0.01) & !near(abs(coef.list), 0, tol = 0.01) & !(abs(theta)>0.9) )

df = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = 0)[selected.ind,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = rep(coef.list[selected.ind], each = 10)) %>% 
  mutate(theta = rep(round((sum.param-coef.list[selected.ind]), digits = 2), each = 10)) %>% 
  mutate(lab = paste(phi, theta, sep = "\n ")) %>%
  mutate(N = as.factor(N)) %>% 
  mutate(TPR = value/nb.sim)

p = ggplot(data = df, aes(x = phi, y = TPR, col = N)) +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lwd = 0.25, col = "black") +
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25) +
  scale_x_continuous(breaks = coef.list[selected.ind], 
                     limits = c((min(df$phi)-0.01), max((df$phi)+0.01)),
                     labels= unique(df$lab)) + 
  scale_y_continuous(breaks = seq(0,1,0.1), 
                     limits = c(0,1)) +
  xlab(param.name) + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5))
  # guides(color=guide_legend(title = param.name, byrow = TRUE))

ggsave(paste0(path_results,"attribution0/auto.arima.TPR.coef.", model.plot,sum.param, "all.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


## Fig 3 detail of the model identification --------------------------
sum.param = 0.6
arma = get(load(file = paste0(path_arima, "performance_autoarima_ARMA1",sum.param,".RData")))

model.plot = "ARMA(1,1)"
model.data = arma
param.name = "phi,theta"
details = 5
# AR, MA
# selected.ind = c(2:7)
# all.model = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = details)[,selected.ind]
# ARMA
theta = sum.param - coef.list
selected.ind = which(!near(coef.list, sum.param, tol = 0.01) & !near(abs(coef.list), 0, tol = 0.01) & !(abs(theta)>0.9) )
all.model = get_data(list.ini = model.data, param.val = 0, true.model = model.plot, details = details)[,selected.ind]


ar.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "AR(1)")))
ma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "MA(1)")))
arma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "ARMA(1,1)")))
white.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "White")))

df <- data.frame(white = white.est, ar = ar.est, ma = ma.est, arma = arma.est, phi = coef.list[selected.ind]) %>%
  # df <- data.frame(white = white.est, ar = ar.est, ma = ma.est, arma = arma.est, phi = coef.list1[selected.ind]) %>% uncomment for AR, MA
  reshape2::melt(id = "phi") %>%
  mutate(predict = as.factor(rep(c( "White", "AR(1)", "MA(1)", "ARMA(1,1)"), each = length(selected.ind)))) %>%
  mutate(theta = round((sum.param-phi), digits = 2)) %>% 
  mutate(lab = paste(phi, theta, sep = "\n ")) %>% 
  mutate(phi = as.factor(phi))
  
p = ggplot(data = df, aes(x = phi, y = value/nb.sim, fill = predict))+
  theme_bw()+
  geom_bar(position="stack", stat="identity", width = 0.5)+
  ylab("Percentage")+
  labs(subtitle = paste0("Model: " , model.plot, ", N = ", length.list[details]))+
  scale_x_discrete(breaks = coef.list[selected.ind], 
                     # limits = c((min(df$phi)-0.01), max((df$phi)+0.01)),
                     labels= unique(df$lab)) + 
  scale_y_continuous(labels = scales::percent,
                     limits = c(-0.05,1)) +
  # geom_text(aes(label = paste0(value*100/nb.sim,"%")), 
  #           position = position_stack(vjust = 0.5), size = 1)+
  ggrepel::geom_text_repel(aes(label = paste0(value*100/nb.sim,"%")), 
                           position = position_stack(vjust = 0.5), size = 1.5, direction = "y", 
                           box.padding = unit(0.01, "lines")) +
  xlab(param.name) + 
  theme(axis.text.x = element_text(size = 6),
        # axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(color = guide_legend(title = param.name))

# AR, MA
# ggsave(paste0(path_results,"attribution0/auto.arima.TPR.detail.",details, model.plot, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)
# ARMA 
ggsave(paste0(path_results,"attribution0/auto.arima.TPR.detail.",details, model.plot,sum.param, ".jpg" ), plot = p, width = 12, height = 5, units = "cm", dpi = 600)

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

model.plot = "ARMA1"
model.true = "ARMA(1,1)"
model.est = "ARMA(1,1)"

if(model.true == model.est){
  est.true = "true"
}else{est.true = "wrong"}

# AR1, MA1
# coef.list = seq(0, 0.6, 0.1)
# selected.ind = c(2:7)
# arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true,".RData")))
# ARMA 
arma.acc = get(load(file = paste0(path_results, "attribution0/truemodel_autoarima_ARMA1", sum.param = 0.3 ,".RData")))
selected.ind = which(!near(coef.list, sum.param, tol = 0.01) & !near(abs(coef.list), 0, tol = 0.01) )[-c(1,2,3)]

df.all = get_data_acc(list.ini = arma.acc, phi = 1, theta =1)
### plot for the true model -------------------------------------------------
param.name = "theta"
leg.title = "phi, theta"
df = df.all[[param.name ]][selected.ind,]  %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(N = length.list) %>% 
  reshape2::melt(id = "N") %>% 
  mutate(phi = rep(coef.list[selected.ind], each = 10)) %>%  
  mutate(theta = round((sum.param-phi), digits = 2)) %>% 
  mutate(lab = paste(phi, theta, sep = ", ")) %>%
  mutate(phi = as.factor(phi))
    
p = ggplot(data = df, aes(x = N, y = value, col = lab))+
  theme_bw()+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000)) + 
  scale_y_continuous(breaks = seq(0, max(df$value)+0.02,0.03),
                     limits = c(0, max(df$value))) + 
  ylab(paste0("Standard error of ", param.name)) +
  labs(subtitle = paste0("Simulated model: " , model.true, ", Estimate model:", model.est))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        plot.tag = element_text(size = 5),
        legend.title.align = 0.5,
        plot.subtitle = element_text(size = 5)) +
  guides(color = guide_legend(title = leg.title)) 

ggsave(paste0(path_results,"attribution0/auto.arima.sdcoef.", model.plot, "_", est.true, param.name, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

### plot for the wrong model -------------------------------------------------
#### AR1 model estimated as MA1 ----------------------------------------------
model.plot = "MA1"
model.true = "MA(1)"
model.est = "AR(1)"

if(model.true == model.est){
  est.true = "true"
}else{est.true = "wrong"}
arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true,".RData")))

df.mean = get_data_acc_mean(list.ini = arma.acc, phi = 1, theta =1)
df.sd = get_data_acc(list.ini = arma.acc, phi = 1, theta =1)
param.name = "phi"
param.name.true = "theta"

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

p = ggplot(data = df, aes(x = N, y = value.y, col = phi))+
  theme_bw()+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  # geom_errorbar(aes( ymin = value.x - value.y, ymax = value.x + value.y)) +
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  # scale_y_continuous(breaks = seq(0, 0.6, 0.1),
  #                    limits = c(0,0.6))+
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


ggsave(paste0(path_results,"attribution0/auto.arima.sd.coef.", model.plot, "_", est.true, param.name, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

#### ARMA model estimated as MA1 ----------------------------------------------
model.plot = "ARMA1a"
model.true = "ARMA(1,1)"
model.est = "MA(1)"

if(model.true == model.est){
  est.true = "true"
}else{est.true = "wrong"}
arma.acc = get(load(file = paste0(path_results, "attribution0/performance_arima_", model.plot, "_", est.true, "_","MA.RData")))

df.mean = get_data_acc_mean(list.ini = arma.acc, phi = 1, theta =1)
df.sd = get_data_acc(list.ini = arma.acc, phi = 1, theta =1)
param.name = "theta"
param.name.true = "phi"

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

p = ggplot(data = df, aes(x = N, y = value.y, col = phi))+
  theme_bw()+
  geom_point(size = 0.3) +
  geom_line(lwd = 0.25)+
  # geom_errorbar(aes( ymin = value.x - value.y, ymax = value.x + value.y)) +
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+ 
  # scale_y_continuous(breaks = seq(0, 0.6, 0.1),
  #                    limits = c(0,0.6))+
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


ggsave(paste0(path_results,"attribution0/auto.arima.sdcoef.", model.plot, "_", est.true, param.name, ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)









# check if the significance test is necessary -----------------------------

phi = sapply(c(1:nb.sim), function(x) arimaorder(fit.i[[x]]))
est.mol = as.data.frame(t(phi))
est.mod = sapply(c(1:nb.sim), function(x) model.iden(as.numeric(est.mol[x,])))
table(est.mod)
ind.s = which(est.mod == "ARMA(1,1)")
phi = as.data.frame(t(sapply(ind.s, function(x) fit.i[[x]]$coef)))
phi.se = as.data.frame(t(sapply(ind.s, function(x) sqrt(diag(fit.i[[x]]$var.coef)))))
t = phi/phi.se




# GAPs --------------------------------------------------------------------

gap.list = seq(0, 0.5, 0.1)
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
model.plot = "AR(1)"
tot.res <- list()
tot.fit = list()

for (i in c(1:length(gap.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=0.6,ma=-0.3), burn.in = burn.in, hetero = 0, sigma = sqrt(sigma.sim), gaps = gap.list[i])
    # fit 
    fit.ar = fit.arima(y.ar)
    TPR.ar[j] = model.iden(fit.ar$pq)
    fit.i[[j]] = fit.ar
  }
  tot.res[[i]] = TPR.ar
  tot.fit[[i]] = fit.i
}


save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_gap_arma.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_gap_arma_all.RData"))

# plot the impact of gaps 
tot.res = get(load(file = paste0(path_results, "attribution0/performance_autoarima_gap_ar.RData")))
tot.fit = get(load(file = paste0(path_results, "attribution0/performance_autoarima_gap_ar_all.RData")))

# TPR ---------------------------------------------------------------------
all.model = as.data.frame(tot.res) 
colnames(all.model) <- NULL

ar.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "AR(1)")))
ma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "MA(1)")))
arma.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "ARMA(1,1)")))
white.est = sapply(c(1:ncol(all.model)), function(x) length(which(all.model[,x] == "White")))

df <- data.frame(white = white.est, ar = ar.est, ma = ma.est, arma = arma.est, gap.per = gap.list) %>%
  # df <- data.frame(white = white.est, ar = ar.est, ma = ma.est, arma = arma.est, phi = coef.list1[selected.ind]) %>% uncomment for AR, MA
  reshape2::melt(id = "gap.per") %>%
  mutate(predict = as.factor(rep(c( "White", "AR(1)", "MA(1)", "ARMA(1,1)"), each = length(gap.list)))) %>%
  # mutate(lab = )
  mutate(gap.per = as.factor(gap.per*100))
nb.sim = 10000
param.name = "Percentage of gap"
p = ggplot(data = df, aes(x = gap.per, y = value/nb.sim, fill = predict))+
  theme_bw()+
  geom_bar(position="stack", stat="identity", width = 0.5)+
  ylab("Percentage")+
  labs(subtitle = paste0("Model: " , model.plot, ", N = 1000"))+
  scale_x_discrete(breaks = c(gap.list*100)) + 
  scale_y_continuous(labels = scales::percent,
                     limits = c(-0.05,1)) +
  # geom_text(aes(label = paste0(value*100/nb.sim,"%")), 
  #           position = position_stack(vjust = 0.5), size = 1)+
  ggrepel::geom_text_repel(aes(label = paste0(value*100/nb.sim,"%")), 
                           position = position_stack(vjust = 0.5), size = 1.5, direction = "y", 
                           box.padding = unit(0.01, "lines")) +
  xlab(param.name) + 
  theme(axis.text.x = element_text(size = 6),
        # axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(color = guide_legend(title = param.name))

ggsave(paste0(path_results,"attribution0/auto.arima.TPR.gap.", model.plot, ".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)

# parameter ---------------------------------------------------------------

coef.dat = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
for (x in c(1:6)) {
  for (y in c(1:nb.sim)) {
    z = NA
    if(as.character(all.model[y,x]) == "AR(1)"){
      z = tot.fit[[x]][[y]]$coef[1]
    } else if(as.character(all.model[y,x]) == "MA(1)"){
      z = tot.fit[[x]][[y]]$coef[3]
    }
    coef.dat[y,x] = z
  }
}

colnames(all.model) = paste0("c", gap.list)
colnames(coef.dat) = paste0("c", gap.list)

df = data.frame(matrix(NA, ncol = 3, nrow = nb.sim*6))
for (i in c(1:6)) {
  df[c((10000*i-9999):(10000*i)), 1] = gap.list[i] 
  df[c((10000*i-9999):(10000*i)), 2] = all.model[,i]
  df[c((10000*i-9999):(10000*i)), 3] = coef.dat[,i] 
}
colnames(df) = c("gap", "model", "coef")
df$gap = as.factor(df$gap)

p = ggplot(data = df, aes(x = gap, y = coef, fill = model))+
  theme_bw()+
  geom_boxplot(size = 0.3, outlier.size = 0.3)+
  ylab("Coefficient")+
  xlab("Percentage of gap")+
  theme(axis.text.x = element_text(size = 6),
        # axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5.5),
        plot.tag = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 5)) 

ggsave(paste0(path_results,"attribution0/auto.arima.TPR.gap.coef.", model.plot, ".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)


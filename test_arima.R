# test arima for missing data coef = data.frame(matrix(NA, nrow = 10000, ncol = 5))
set.seed(1)
l = 1000
classic.rho <- function(x){
  l1 = length(x)
  u = x[1:(l1-1)]
  v = x[2:l1]
  t1 = u - mean(u)
  t2 = v- mean(v)
  nor = sum(t1*t2)
  denor = sqrt((sum(t1**2)) *(sum(t2**2)))
  rho = nor/denor
  return(rho)
}
for (i in c(1:10000)) {
  gap = sample(c(1:l), size = 200, replace = F)
  y = arima.sim(model = list(ar=0.5), n = l, sd = 1 )
  y.na = y
  y.na[gap] = NA
  y.wtna = y[-gap]
  
  fit1 = arima(y.na, order = c(1,0,0))
  fit2 = arima(y.wtna, order = c(1,0,0))
  fit3 = arima(y, order = c(1,0,0))
  
  coef[i,] = c(classic.rho(y), fit3$coef[1], classic.rho(y.wtna), fit2$coef[1],  fit1$coef[1])
  
}

colnames(coef) <- c("F.classic", "F.arima", "M.classic", "M.arima.m", "M.arima.f")
df <- coef[,-4] %>% 
  reshape2::melt() %>%
  mutate(simulate = c(rep("Without NA",20000), rep("With NA", 20000))) %>%
  mutate(estimator = as.factor(rep(c("classical", "arima", "classical", "arima"), each = 10000)))

p = ggplot(data = df, aes(x = estimator, y = value, col = simulate))+theme_bw()+
  geom_boxplot(lwd = 0.3, outlier.size = 0.3)+
  geom_hline(yintercept = 0.5, lwd = 0.3)+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=4),
        legend.title = element_text(size = 4.5),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        legend.title.align = 0.5,
        plot.subtitle = element_text(size = 4))

ggsave(paste0(path_results,"attribution0/arima.treat.NA.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)







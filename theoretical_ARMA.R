library("ggplot2")
# make plot of theoretical ACF and PACF of 3 models  ----------------------
lagmax = 14
phi = 0.3
theta = 0.3
phi1 = 0.6
theta1 = -0.3
phi2 = 0.3
theta2 = 0.2
Lag = c(0:lagmax)
# compute ACF and PACF
## AR(1)
acf.ar = ARMAacf(ar = phi, ma = 0, lag.max = lagmax, pacf = FALSE)
pacf.ar = c(0,ARMAacf(ar = phi, ma = 0, lag.max = lagmax, pacf = TRUE))
## MA(1)
acf.ma = ARMAacf(ar = 0, ma = theta, lag.max = lagmax, pacf = FALSE)
pacf.ma = c(0,ARMAacf(ar = 0, ma = theta, lag.max = lagmax, pacf = TRUE))
## ARMA(1,1)a
acf.arma.a = ARMAacf(ar = phi1, ma = theta1, lag.max = lagmax, pacf = FALSE)
pacf.arma.a = c(0,ARMAacf(ar = phi1, ma = theta1, lag.max = lagmax, pacf = TRUE))
## ARMA(1,1)b
acf.arma.b = ARMAacf(ar = phi2, ma = theta2, lag.max = lagmax, pacf = FALSE)
pacf.arma.b = c(0,ARMAacf(ar = phi2, ma = theta2, lag.max = lagmax, pacf = TRUE))
# PLOT
df.plot = data.frame(Lag = Lag, acf.ar = acf.ar, acf.ma = acf.ma, acf.arma.a = acf.arma.a, acf.arma.b = acf.arma.b,
                     pacf.ar = pacf.ar, pacf.ma = pacf.ma, pacf.arma.a = pacf.arma.a, pacf.arma.b = pacf.arma.b) %>% 
  reshape2::melt(id = "Lag") %>% 
  mutate(type = rep(c("ACF", "PACF"), each = 4*(lagmax+1))) %>%
  mutate(model = rep(rep(c("AR(1)", "MA(1)", "ARMA(1,1)a",  "ARMA(1,1)b"), each = (lagmax+1)), 2)) %>%
  mutate(model = factor(model, levels = c("AR(1)", "MA(1)", "ARMA(1,1)a",  "ARMA(1,1)b")))

p <- ggplot(df.plot, aes(x = Lag, y = value)) + 
  theme_bw() +
  geom_col(color = "steelblue", width = 0.2) + 
  facet_grid(type ~ model) +
  scale_x_continuous(breaks = seq(0,15,3)) +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6),
        strip.background = element_rect(fill="white"),
        plot.tag = element_text(size = 6),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 6))
ggsave(paste0(path_results,"attribution/ACF_PACF_theretical.jpg" ), plot = p, width = 16, height = 8, units = "cm", dpi = 1200)

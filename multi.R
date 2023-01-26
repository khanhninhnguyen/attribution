# study the multicolinearity 

n = 730
y = rnorm(n, 0,1)
y[(n/2):n] <- y[(n/2):n] 
df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:n])
Data.mod = construct.design(data.df = df, name.series = "y", break.ind = n/2)
Data.mod = Data.mod[,-11]
list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
mod.X <-  list.para %>% stringr::str_c(collapse = "+")
mod.expression <- c("signal","~",mod.X) %>% stringr::str_c(collapse = "")
# ols
ols.fit = lm(mod.expression, data = Data.mod)
car::vif(ols.fit)
# adjust the p values due to the dependence 
# this program is to explore the test result 

Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
all.var = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
all.arima = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)

# convert to p value 

t.table = Total.res[,paste0("t", list.name.test)]
ttop = function(x){
  reut2*pnorm(-abs(t.table[x,]), mean = 0, sd = 1, lower.tail = TRUE)
}
p.table = apply(t.table, 2, function(x) )
